

#!/bin/tcsh -xef



# ==================== { fmri_prep_hj.sh } ========================================
# ||
# ||  [Purpose]: For a given subject, do the preprocessing of MVPA blocks from DeMo :
# ||
# ||             -----[   1.  Despiking  ]
# ||             -----[   2.  Slicetiming  ]
# ||             -----[   3.  Motion Correction  ]
# ||             -----[   4.  EPI unwarping  ]
# ||             -----[   5.  Coregistration  ]
# ||             -----[   6.  Normalization  ]
# ||             -----[   7.  Smoothing  ]
# ||             -----[   8.  Scaling  ]
# ||
# ||  [Usage]: - to execute via tcsh: 
# ||              >   tcsh -xef fmri_prep_hj.sh 0xx |& tee output.fmri_prep_hj.0xx.log.sh
# ||           - to execute via bash: 
# ||              >   tcsh -xef fmri_prep_hj.sh 0xx 2>&1 | tee output.fmri_prep_hj.oxx.log.sh
# ||
# ||           - modify some parameters and the preprocessing blocks in the script if needed 
# ||
# ||  [Notice]: a)
# ||
# ===========================================================================================
# Author: Shuo, Chen [shuo.chen2@mail.mcgill.ca]
#   Date: 2023.0710 [tcsh shell]
#       @ JD lab, BIC of MNI, McGill University, Canada


echo "(version 7.10, February 18, 2019)"
echo "execution started: `date`"



# =========================== auto block: setup ============================
# script setup

# take note of the AFNI version
afni -ver

# check that the current AFNI version is recent enough
afni_history -check_date 27 Jun 2019
if ( $status ) then
    echo "** this script requires newer AFNI binaries (than 27 Jun 2019)"
    echo "   (consider: @update.afni.binaries -defaults)"
    exit
endif

# the user may specify a single subject to run with
if ( $#argv > 0 ) then
    set subj = sub-$argv[1]
else
    echo "-- < ! > -- please provide a subject ID"
    exit 1
endif

# assign output directory name

set afni_dir = "/Users/shuo/abin"
set script_dir = "/Users/shuo/OneDrive - McGill University/Git_Codes/DeMo" 
set root_dir = /Volumes/JDlab_Shuo/DeMo/Hannah
set ses_dir = MVPA

set data_dir = $root_dir/$subj/ses-$ses_dir/func
set T1_dir = $root_dir/$subj/ses-D1/anat
set fmap_dir = $root_dir/$subj/ses-$ses_dir/fmap


# set list of task runs
set tasks = (block1 block2) # block2 block3 block4 block5
set grp_base = block2

set output_dir = $root_dir/derivative/$subj/ses-$ses_dir/func
# set output_t1 = $root_dir/derivative/$subj/ses-D1/anat

set block2do = (9 10) #(0 1 2 3 4 5 6 7 8 9 10 cs)

set versions2do = (DySy DySn) #(DySy DySn DnSy DnSn)

# some parameters used in the following steps
set CMthreshold = 6 # mm ; if center distance > threshold, then do the center matching before the coregistration between EPI sessions
set Motion_Outlier = 0.3 # threshold for outliers, calculated by afni enorm
set Smth_kernel = 5 # mm


# verify that the results directory does not yet exist
if ( -d $output_dir ) then
    echo output dir "$output_dir" already exists
    # exit
else
	mkdir -p $output_dir
	# mkdir -p $output_t1

endif

# -------------------------------------------------------
# enter the results directory (can begin processing data)
cd $data_dir

# get the length of scans of each task
set tr_counts = ()
foreach run ( $tasks )
    set tr_counts = ( $tr_counts `3dinfo -nv *$run*.nii.gz` )
end

# echo $tr_counts

# ============================ block 00 : radcor ============================


# # apply 3dTcat to copy input dsets to results dir,
# # while removing the first 0 TRs
# 3dTcat -prefix $output_dir/pb00.$subj.r01.tcat                   \
#     sub-027_ses-D2_task-rest1_part-mag_bold.nii.gz'[0..$]'
# 3dTcat -prefix $output_dir/pb00.$subj.r02.tcat                   \
#     sub-027_ses-D2_task-testA_part-mag_bold.nii.gz'[0..$]'

# # and make note of repetitions (TRs) per run
# set tr_counts = ( 280 135 )


if (" $block2do " =~ *" 0 "*) then
    echo "Proprocessing block 00 is ongoing."

    # ---------------------------------------------------------
    # data check: compute correlations with spherical ~averages
    rm -rf $output_dir/pb00.radcor

    @radial_correlate -nfirst 0 -do_clean yes -rdir $output_dir/pb00.radcor \
    	                  *.nii.gz

else
    echo "Proprocessing block 00 is skipped!"
endif




# ========================== block 01 : outcount ==========================

if (" $block2do " =~ *" 1 "*) then
    echo "Proprocessing block 01 is ongoing."

    rm -rf $output_dir/pb01.outcount
    mkdir $output_dir/pb01.outcount

    # data check: compute outlier fraction for each volume
    foreach run ( $tasks )
        touch $output_dir/pb01.outcount/out.task-$run.pre_ss_warn.txt

        3dToutcount -automask -fraction -polort 3 -legendre                     \
                    *$run*.nii.gz > $output_dir/pb01.outcount/outcount.task-$run.1D

        # outliers at TR 0 might suggest pre-steady state TRs
        if ( `1deval -a $output_dir/pb01.outcount/outcount.task-$run.1D"{0}" -expr "step(a-0.4)"` ) then
            echo "** TR #0 outliers: possible pre-steady state TRs in task-$run" \
                >> $output_dir/pb01.outcount/out.task-$run.pre_ss_warn.txt
        endif
    end

    # catenate outlier counts into a single time series
    cat $output_dir/pb01.outcount/outcount.task-$run.1D > $output_dir/pb01.outcount/outcount.taskALL.1D

    # get run number and TR index for minimum outlier volume
    set minindex = `3dTstat -argmin -prefix - $output_dir/pb01.outcount/outcount.taskALL.1D\'`
    echo $minindex
    set ovals = ( `1d_tool.py -set_run_lengths $tr_counts                       \
                              -index_to_run_tr $minindex` )


    # save run and TR indices for extraction of vr_base_min_outlier
    set minoutrun = $ovals[1]
    set minouttr  = $ovals[2]
    echo "min outlier: run $minoutrun, TR $minouttr" | tee $output_dir/pb01.outcount/out.min_outlier.txt


else
    echo "Proprocessing block 01 is skipped!"
endif


# ================================ block 02 : despike =================================
if (" $block2do " =~ *" 2 "*) then
    echo "Proprocessing block 02 is ongoing."

    rm -rf $output_dir/pb02.despike
    mkdir $output_dir/pb02.despike

    # apply 3dDespike to each run
    foreach run ( $tasks )
        3dDespike -NEW -prefix $output_dir/pb02.despike/pb02.$subj.task-$run.despike.nii.gz \
            *$run*.nii.gz
    end
else
    echo "Proprocessing block 02 is skipped!"
endif


# ================================= block 03 : tshift =================================
if (" $block2do " =~ *" 3 "*) then
    echo "Proprocessing block 03 is ongoing."

    rm -rf $output_dir/pb03.tshift
    mkdir $output_dir/pb03.tshift
    cd $output_dir/pb03.tshift

    # time shift data so all slice timing is the same 
    foreach run ( $tasks )

        # slice timing correction
        # jsonfile="`ls -S *$run*.json | head -1`"
        abids_json_tool.py -json2txt -overwrite \
             -input $data_dir/*$run*.json -prefix pb03.junk.txt

        grep SliceTiming pb03.junk.txt | sed -e 's/^SliceTiming *://' > SliceTimes.1D
        rm pb03.junk.txt

        # The default shifted target time is the average of the 'tpattern' values
        3dTshift -TR 1.5 -wsinc9 -tpattern @SliceTimes.1D \
                 -prefix pb03.$subj.task-$run.tshift.DySy.nii.gz \
                 $output_dir/pb02.despike/*$run*nii.gz

        # Prepare the version without despiking
        3dTshift -TR 1.5 -wsinc9 -tpattern @SliceTimes.1D \
                 -prefix pb03.$subj.task-$run.tshift.DnSy.nii.gz \
                 $data_dir/*$run*nii.gz

        # # ref:
        # https://reproducibility.stanford.edu/slice-timing-correction-in-fmriprep-and-linear-modeling/
        # https://wiki.humanconnectome.org/display/PublicData/HCP+fMRI+slice-timing+acquisition+parameters
        # https://en.wikibooks.org/wiki/SPM/Slice_Timing#cite_note-sladky-2011-2

        3dcopy $output_dir/pb02.despike/*$run*nii.gz pb03.$subj.task-$run.tshift.DySn.nii.gz
        rm $output_dir/pb02.despike/*$run*nii.gz
        3dcopy $data_dir/*$run*nii.gz pb03.$subj.task-$run.tshift.DnSn.nii.gz

    end

else
    echo "Proprocessing block 03 is skipped!"
endif

# ================================= block 04: motion =================================
if (" $block2do " =~ *" 4 "*) then
    echo "Proprocessing block 04 is ongoing."

    rm -rf $output_dir/pb04.motion
    mkdir $output_dir/pb04.motion
    cd $output_dir/pb04.motion

    foreach run ( $tasks )

        cd $output_dir/pb01.outcount

        # # extract MIN_OUTLIER index for current run
        # set min_outlier_index = `3dTstat -argmin -prefix - outcount.task-$run.1D\'`

        # cd $output_dir/pb04.motion

        # # extract volreg base for this run
        # 3dbucket -prefix intra_base_task-$run.DySy.nii.gz                         \
        #     "$output_dir/pb03.tshift/pb03.$subj.task-$run.tshift.DySy.nii[$min_outlier_index]"
              

        cd $output_dir/pb04.motion
                   

        # get the index of this run in a given session
        set r_index = 1
        foreach run_search ($tasks)
            if ($run_search == $run) then
                echo "-- motion correcting task-$tasks[$r_index] now --"
                break
            endif
            @ r_index++
        end

        3dcopy $fmap_dir/*sbref_dir-PA_run-0$r_index*.nii.gz sref-for.task-$run.nii.gz

        foreach ver ( $versions2do )
            # intra-EPI: register each volume to the base image
            3dvolreg -zpad 4 -twopass -heptic                                \
                     -base sref-for.task-$run.nii.gz                           \
                     -1Dfile motion.task-$run.$ver.1D \
                     -prefix rm.volreg.task-$run.$ver.nii.gz   \
                     -1Dmatrix_save mat.task-$run.vr.intra.$ver.1D                         \
                     $output_dir/pb03.tshift/*$run*$ver.nii.gz

            unset r_index run_search

            # create a temporal mask to mark down huge instant (delta/relative) movement volumes; motion censoring/scrubbing
            1d_tool.py -infile mat.task-$run.vr.intra.$ver.1D \
                    -set_nruns 1 -show_censor_count \
                    -censor_prev_TR \
                    -censor_motion $Motion_Outlier task-$run.$ver # -censor_next_TR -censor_prev_TR

            # try stricter threshold as 0.2 for both sleep and task sessions. More scrubbing no harm, as long we have enough samples for both training and prediction
            # https://afni.nimh.nih.gov/afni/community/board/read.php?1,161451,161463#msg-161463
            # https://afni.nimh.nih.gov/afni/community/board/read.php?1,143809,144115
            # <if scrubbed out beyond 50% volumes, discard the session>

            # in case other indice for calculating to-be-scrubbed volumes are required
              # fsl_motion_outliers -i "${data}.nii.gz" -o Mo_DV_ConfM_${subj}_${session}.txt --dvars -s Mo_DV_raw_${subj}_${session}.txt -p Mo_DV_Plot_${subj}_${session}  #  based on intensity differences within the realigned timeseries; In the case that the motion correction is not accurate, then using motion correction parameters (rotation angles and translations in mm) is a poor way to estimate the outliers.
              # fsl_motion_outliers -i "${data}.nii.gz" -o Mo_FD_ConfM_${subj}_${session}.txt --fd -s Mo_FD_raw_${subj}_${session}.txt -p Mo_FD_Plot_${subj}_${session} # Power
              # fsl_motion_outliers -i "${data}.nii.gz" -o Mo_FDrms_ConfM_${subj}_${session}.txt --fdrms -s Mo_FDrms_raw_${subj}_${session}.txt -p Mo_FDrms_Plot_${subj}_${session} # Jenkinson
              # <compare fsl --refrms --fdrms and afni enorm .3>
        end
    end
else
    echo "Proprocessing block 04 is skipped!"
endif


# ================================== block 05 : unwarp ==================================
# compute blip up/down non-linear distortion correction for EPI

if (" $block2do " =~ *" 5 "*) then
    echo "Proprocessing block 05 is ongoing."

    rm -rf $output_dir/pb05.unwarp
    mkdir $output_dir/pb05.unwarp
    cd $output_dir/pb05.unwarp


    foreach run ( $tasks )

        # get the index of this run in a given session
        set r_index = 1
        foreach run_search ($tasks)
            if ($run_search == $run) then
                echo "-- unwarping $tasks[$r_index] now --"
                break
            endif
            @ r_index++
        end

        3dcopy $fmap_dir/*sbref_dir-PA_run-0$r_index*.nii.gz sref-for.task-$run
        3dcopy $fmap_dir/*sbref_dir-AP_run-0$r_index*.nii.gz sref-rev.task-$run


        # if the centers are far, roughly align the center first.
        set CMDist = `@Center_Distance -dset sref-for.task-$run+orig sref-rev.task-$run+orig`

        echo
        echo "|| ----- AP-PA mass-center distance = $CMDist mm-----"
        echo

        if (`echo "${CMDist} > ${CMthreshold}" | bc` == 1) then
            echo
            echo "|| ----- align (rigid body) the AP PA image first, because they are far at the beginning -----"
            echo

            3dcopy sref-rev.task-$run+orig sref-rev.far.task-$run+orig
            rm *sref-rev.task-$run+orig

            # intra-EPI: register each volume to the base image
            3dvolreg -zpad 4 -twopass -Fourier                                \
                     -base sref-for.task-$run+orig                           \
                     -1Dfile rm.task-$run.blip.vr.1D -prefix sref-rev.task-$run+orig   \
                     sref-rev.far.task-$run+orig


        endif

        # automask the single band reference datasets 
        3dAutomask -apply_prefix rm.sref-for.masked.task-$run.nii.gz sref-for.task-$run+orig
        3dAutomask -apply_prefix rm.sref-rev.masked.task-$run.nii.gz sref-rev.task-$run+orig

        # compute the midpoint warp between the reference datasets
        3dQwarp -plusminus -pmNAMES rev for                           \
                -pblur 0.05 0.05 -blur -1 -1                          \
                -noweight -minpatch 9                                 \
                -source rm.sref-rev.masked.task-$run.nii.gz                     \
                -base   rm.sref-for.masked.task-$run.nii.gz                     \
                -prefix blip.task-$run

        # warp reference datasets (forward and each masked) for QC checks
        # (and preserve obliquity)
        3dNwarpApply -quintic -nwarp blip.task-"$run"_for_WARP+orig          \
                     -source sref-for.task-$run+orig                                 \
                     -prefix unwarp.sref-for.task-$run.nii.gz
        3drefit -atrcopy sref-for.task-$run+orig IJK_TO_DICOM_REAL                   \
                         unwarp.sref-for.task-$run.nii.gz

        # 3dNwarpApply -quintic -nwarp blip.task-"$run"_for_WARP+orig          \
        #              -source rm.sref-for.masked.task-$run.nii.gz                \
        #              -prefix unwarp.sref-for.masked.task-$run.nii.gz
        # 3drefit -atrcopy sref-for.task-$run+orig IJK_TO_DICOM_REAL                   \
        #                  unwarp.sref-for.masked.task-$run.nii.gz

        # 3dNwarpApply -quintic -nwarp blip.task-"$run"_rev_WARP+orig          \
        #              -source rm.sref-rev.masked.task-$run.nii.gz                \
        #              -prefix unwarp.sref-rev.masked.task-$run.nii.gz
        # 3drefit -atrcopy sref-rev.task-$run+orig IJK_TO_DICOM_REAL                   \
        #                  unwarp.sref-rev.masked.task-$run.nii.gz

        # warp EPI time series data
        foreach ver ( $versions2do )

            3dNwarpApply -quintic -nwarp blip.task-"$run"_for_WARP+orig          \
                         -source $output_dir/pb04.motion/rm.volreg.task-$run.$ver.nii.gz         \
                         -prefix rm.unwarp.task-$run.$ver.nii.gz
            3drefit -atrcopy sref-for.task-$run+orig IJK_TO_DICOM_REAL                   \
                             rm.unwarp.task-$run.$ver.nii.gz
        
        end
        unset r_index run_search

    end
else
    echo "Proprocessing block 05 is skipped!"
endif



# ================================= block 06 : coregistration =================================
if (" $block2do " =~ *" 6 "*) then
    echo "Proprocessing block 06 is ongoing."

    rm -rf $output_dir/pb06.coreg
    mkdir $output_dir/pb06.coreg
    cd $output_dir/pb06.coreg

    foreach run ( $tasks )

        if (" $run " =~ *$grp_base*) then

            echo "Run $grp_base was chosen as the group reference, no need to coregistrate to itself"
            3dcopy $output_dir/pb05.unwarp/unwarp.sref-for.task-$grp_base.nii.gz \
                sref-for.task-$grp_base.grp_base.nii.gz


            # Just do once for e2a: compute anat alignment transformation to EPI registration base
            # (new anat will be intermediate, stripped, sub-0xx_ses-D1_T1w_ns+orig)

            # run uniformity correction on EPI base
            3dUnifize -T2 -input sref-for.task-$grp_base.grp_base.nii.gz \
                -prefix sref-for.task-$grp_base.grp_base.unif.nii.gz

            3dcopy $T1_dir/*T1w.nii.gz T1w.nii.gz

            align_epi_anat.py -anat2epi -anat $T1_dir/*T1w.nii.gz     \
                   -save_skullstrip -suffix _cr                  \
                   -epi sref-for.task-$grp_base.grp_base.unif.nii.gz  \
                   -epi_base 0     \
                   -epi_strip 3dAutomask                              \
                   -cost lpc+ZZ -giant_move                           \
                   -volreg off -tshift off -ex_mode quiet

            # optional tunig parameters
            # align_epi_anat.py -epi2anat -ginormous_move -anat_has_skull no  \
            #        -save_resample -align_centers on -Allineate_opts '-weight_frac 1.0 -maxrot 6 -maxshf 10 -VERB -warp shr'

        else

            # and compute xforms for cross-run allin to EPIcr_base
            3dAllineate -base $output_dir/pb05.unwarp/unwarp.sref-for.task-$grp_base.nii.gz                                \
                        -source $output_dir/pb05.unwarp/unwarp.sref-for.task-$run.nii.gz                            \
                        -1Dfile EPIcr_allin_dfile.task-$run.1D                      \
                        -1Dmatrix_save mat.EPIcr_allin.task-$run.1D               \
                        -autoweight -source_automask                                  \
                        -lpa -cubic -prefix rm.coreg.sref-for.task-$run.nii.gz  
        endif

    end

else
    echo "Proprocessing block 06 is skipped!"
endif


# ================================= block 07 : normalization + segmentation =================================
if (" $block2do " =~ *" 7 "*) then
    echo "Proprocessing block 07 is ongoing."

    rm -rf $output_dir/pb07.norm+seg
    mkdir $output_dir/pb07.norm+seg
    cd $output_dir/pb07.norm+seg

    3dcopy $afni_dir/MNI152_T1_2009c+tlrc rm.mnit1_09c_1mm.nii.gz

    # warp anatomy to standard space (non-linear warp)
    auto_warp.py -base /Users/shuo/abin/MNI152_T1_2009c+tlrc \
                 -input $output_dir/pb06.coreg/*ns+orig.HEAD                       \
                 -skull_strip_input yes

    # move results up out of the awpy directory
    # - NL-warped anat, affine warp, NL warp
    # - use typical standard space name for anat
    # - be sure NIFTI sform_code=2 means standard space
    3dbucket -DAFNI_NIFTI_VIEW=tlrc                                 \
             -prefix T1w@MNI.nii.gz awpy/*T1w*.aw.nii*

    mv awpy/anat.un.aff.Xat.1D .
    mv awpy/anat.un.aff.qw_WARP.nii .

    # segmentation
    3dSeg -anat $output_dir/pb06.coreg/*ns+orig.HEAD -mask AUTO -classes 'CSF ; GM ; WM'


else
    echo "Proprocessing block 07 is skipped!"
endif


# ================================= block 08 : EPI one-step transformation =================================
if (" $block2do " =~ *" 8 "*) then
    echo "Proprocessing block 08 is ongoing."

    rm -rf $output_dir/pb08.EPItrans
    mkdir $output_dir/pb08.EPItrans
    cd $output_dir/pb08.EPItrans

    foreach run ( $tasks )

        if (" $run " =~ *$grp_base*) then

            # catenate epi2anat/tlrc xforms, because there is no EPI coregistration for this run
            cat_matvec -ONELINE                                                       \
                       $output_dir/pb07.norm+seg/anat.un.aff.Xat.1D                                             \
                       $output_dir/pb06.coreg/*T1w_cr_mat.aff12.1D -I > mat.task-$run.warp.aff12.1D

        else

            # catenate EPIcoreg/epi2anat/tlrc xforms
            cat_matvec -ONELINE                                                       \
                       $output_dir/pb07.norm+seg/anat.un.aff.Xat.1D                                             \
                       $output_dir/pb06.coreg/*T1w_cr_mat.aff12.1D -I                     \
                       $output_dir/pb06.coreg/mat.EPIcr_allin.task-$run.1D > mat.task-$run.warp.aff12.1D

        endif


        foreach ver ( $versions2do )

            # apply catenated xform: volreg/unwarp/EPIcoreg/epi2anat/tlrc/NLtlrc
            # then apply non-linear standard-space warp
            3dNwarpApply -master $output_dir/pb07.norm+seg/T1w@MNI.nii.gz -dxyz 2.5                 \
                         -source $output_dir/pb03.tshift/*.task-$run.*DySy.nii.gz                     \
                         -nwarp "$output_dir/pb07.norm+seg/anat.un.aff.qw_WARP.nii \
                         mat.task-$run.warp.aff12.1D \
                         $output_dir/pb05.unwarp/blip.task-"$run"_for_WARP+orig \
                         $output_dir/pb04.motion/mat.task-$run.vr.intra.$ver.1D" \
                         -ainterp wsinc5                                              \
                         -prefix rm.epi.nomask.task-$run.$ver.nii.gz

            # create an all-1 dataset to mask the extents of the warp
            3dcalc -overwrite -a $output_dir/pb05.unwarp/rm.unwarp.task-$run.$ver.nii.gz -expr 1                   \
                   -prefix rm.epi.all1.task-$run.$ver.nii.gz

            # warp the all-1 dataset for extents masking 
            3dNwarpApply -master $output_dir/pb07.norm+seg/T1w@MNI.nii.gz -dxyz 2.5                 \
                         -source rm.epi.all1.task-$run.$ver.nii.gz                                     \
                         -nwarp "$output_dir/pb07.norm+seg/anat.un.aff.qw_WARP.nii \
                         mat.task-$run.warp.aff12.1D" \
                         -interp cubic                                                \
                         -ainterp NN -quiet                                           \
                         -prefix rm.epi.1.task-$run.$ver

            # make an extents intersection mask of this run
            3dTstat -min -prefix rm.epi.min.task-$run.$ver rm.epi.1.task-"$run".$ver+tlrc
        end

    end

    # warp the volreg base EPI dataset to make a final version
    3dNwarpApply -master $output_dir/pb07.norm+seg/T1w@MNI.nii.gz -dxyz 2.5                     \
                 -source $output_dir/pb05.unwarp/sref-for.task-$grp_base+orig                                 \
                 -nwarp "$output_dir/pb07.norm+seg/anat.un.aff.qw_WARP.nii \
                 mat.task-$grp_base.warp.aff12.1D \
                 $output_dir/pb05.unwarp/blip.task-"$grp_base"_for_WARP+orig" \
                 -ainterp wsinc5                                                  \
                 -prefix sref-for.task-$grp_base.grp_base@MNI.nii.gz


    foreach ver ( $versions2do )
        # ----------------------------------------
        # create the extents mask: mask_epi_extents+tlrc
        # (this is a mask of voxels that have valid data at every TR)
        3dMean -datum short -prefix rm.epi.mean.$ver rm.epi.min.task-*.$ver*.HEAD 
        3dcalc -a rm.epi.mean.$ver+tlrc -expr 'step(a-0.999)' -prefix mask_epi_extents.$ver

        # and apply the extents mask to the EPI data 
        # (delete any time series with missing data)
        foreach run ( $tasks )
            3dcalc -a rm.epi.nomask.task-$run.$ver.nii.gz -b mask_epi_extents.$ver+tlrc               \
                   -expr 'a*b' -prefix pb08.task-$run.EPItrans.$ver.nii.gz
        end


        # ---------------------------------------------------------
        # data check: compute correlations with spherical ~averages
        @radial_correlate -nfirst 0 -do_clean yes -rdir radcor.post.$ver           \
                          pb08.task-*.EPItrans.$ver.nii.gz
    end

else
    echo "Proprocessing block 08 is skipped!"
endif


# ================================= block 09 : Smoothing =================================

if (" $block2do " =~ *" 9 "*) then
    echo "Proprocessing 09 is ongoing."

    rm -rf $output_dir/pb09.smooth
    mkdir $output_dir/pb09.smooth
    cd $output_dir/pb09.smooth


    # # prepare a mask for below smoothing

    3dmask_tool -input $afni_dir/MNI152_T1_2009c+tlrc \
            -prefix "rm.anat_mask.nii.gz" \
            -fill_holes -fill_dirs xy -dilate_input 1 -1

    foreach ver ( $versions2do )    
        3dresample -master $output_dir/pb08.EPItrans/*$grp_base.EPItrans.$ver.nii.gz \
                -prefix "rm.anat_mask.rs.$ver.nii.gz" \
                -inset "rm.anat_mask.nii.gz"

        foreach run ( $tasks )
            # Do the smoothing
            3dBlurInMask -input $output_dir/pb08.EPItrans/*$run.EPItrans.$ver.nii.gz \
                    -FWHM $Smth_kernel \
                    -mask rm.anat_mask.rs.$ver.nii.gz \
                    -float -quiet -prefix pb09.task-$run.smooth.$ver.nii.gz
        end
    end
else
    echo "Proprocessing block 09 is skipped!"
endif

# ================================= block 10 : Scaling =================================
# scale each voxel time series to have a mean of 100
# (be sure no negatives creep in)
# (subject to a range of [0,200])

if (" $block2do " =~ *" 10 "*) then
    echo "Proprocessing 10 is ongoing."

    rm -rf $output_dir/pb10.scale
    mkdir $output_dir/pb10.scale
    cd $output_dir/pb10.scale

    foreach ver ( $versions2do )
        foreach run ( $tasks )
            3dTstat -prefix rm.mean_task-$run.$ver.nii.gz $output_dir/pb09.smooth/*$run.smooth.$ver.nii.gz

            3dcalc -a $output_dir/pb09.smooth/*$run.smooth.$ver.nii.gz \
                   -b rm.mean_task-$run.$ver.nii.gz \
                   -c $output_dir/pb08.EPItrans/mask_epi_extents.$ver+tlrc                            \
                   -expr 'c * min(200, a/b*100)*step(a)*step(b)'       \
                   -prefix pb10.task-$run.scale.$ver.nii.gz

           echo "-- < ! > -- The maximum value in data : `ls pb10.task-$run.scale.$ver.nii.gz` \
                is     `3dinfo -dmaxus pb10.task-$run.scale.$ver.nii.gz`"

        end
    end

else
    echo "Proprocessing block 10 is skipped!"
endif

# ================================= block: cleaning the space =================================
if (" $block2do " =~ *" cs "*) then
    echo "The intermediate outcome will be cleaned after the quality-check is done." 

    read -p "-- < ! > -- If the quality-check is done, press [Enter] key, \
            then space will be cleaned out."

    cd $output_dir/
    find . -type f -name 'rm.*' -delete

    # take the advantage of the naming short-cut, as the final outcome are all named after in pb10.*
    find . -type f -name 'pb0*.*.nii.gz' -delete
else
    echo "No space-cleaning was executed"
endif



# ================================= block: for developing =================================

if (" $block2do " =~ *" 99 "*) then
    echo "Proprocessing TEST is ongoing."


    # cd $output_dir/pb08.EPItrans

    # # create 'full_mask' dataset (union mask)
    # foreach run ( $tasks )
    #     3dAutomask -prefix rm.mask_task-$run pb08.task-$run.EPItrans.nii.gz
    # end

    # # create union of inputs, output type is byte
    # 3dmask_tool -inputs rm.mask_task-*+tlrc.HEAD -union -prefix full_mask.$subj

    # # ---- create subject anatomy mask, mask_anat.$subj+tlrc ----
    # #      (resampled from tlrc anat)
    # 3dresample -master full_mask.$subj+tlrc -input $output_dir/pb07.norm+seg/T1w@MNI.nii.gz \
    #            -prefix rm.resam.anat

    # # convert to binary anat mask; fill gaps and holes
    # 3dmask_tool -dilate_input 5 -5 -fill_holes -input rm.resam.anat+tlrc      \
    #             -prefix mask_anat.$subj

    # # compute tighter EPI mask by intersecting with anat mask
    # 3dmask_tool -input full_mask.$subj+tlrc mask_anat.$subj+tlrc              \
    #             -inter -prefix mask_epi_anat.$subj



endif




# ================================= the end =================================


cd "$script_dir"
mv output.*.sh $output_dir/`ls output.*.sh`
echo "execution finished: `date`"














# notes for Shuo
# 1. avoid multiple temporal projection, or try to remove the effect from the previous step in a sequential order.
# 3. exact length of trs, extract scans before black screen ending switch 
# 4. think about the naming of files across sessions
# 5. 027 mvpa scans 45 degree tilt
# 8. test 3dvolreg for sleep data, with heptic+2pass vs Fourier
# 9. change the motion correction reference, if single band ref is blurring
# 10. masks erode? resample?

# quality-control check ordes:
# 1. steady, time series, 3dToutcount
# 2. film play the data, notice big motion, and data quality
# 3. check the all the single references, choose which to use as the group ref
