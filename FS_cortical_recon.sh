#!/bin/bash

##########################################################################################
# Prepare_for_Electrode_Reconstruction                                                   #
# This script runs a cortical reconstruction through Freesurfer and converts the files to#
# run the Hermes method                                                                  #
#                                                                                        #
# Michael Randazzo 8/6/14                                                                #
#                                                                                        #
##########################################################################################

#Step 0
#Setup freesurfer
export FREESURFER_HOME=/Applications/freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh
#reset FS_LOAD_DWI to avoid error
export FS_LOAD_DWI=0

# 1st Step
# Acquire CT and Pre-Op MRI images and place them in the Subjects folder (Use OsiriX to
# separate and export the desired scans from the remaining DICOM images)
# Assume CT file structure SUBJECT_ID/Anatomy/CT/*/(images)
# Assume MRI file structure SUBJECT_ID/Anatomy/MRI/preop/(images)

# Set the subject ID and the subject directory below:
export SUBJECT_ID=DBS40712
export SUBJECTS_DIR=/Users/Dengyu/Documents/Subjects/${SUBJECT_ID}
# export SUBJECTS_DIR=/Volumes/Nexus/DBS/${SUBJECT_ID}/Anatomy

#####

# 2nd Step
# Reconstruct the cortical surface using Freesurfer

# Set Subject ID
#export SUBJECT_ID=preop

# Cortical Reconstruction (Step will take an estimated 12-18 hours)
recon-all -subject ${SUBJECT_ID}_FS -i /Users/Dengyu/Documents/Subjects/${SUBJECT_ID}/MRI/preop/preop_used/IM-0015-0001.dcm -autorecon-all -notal-check

#####

# 3rd Step
# Convert Freesurfer output to image files
# Also convert previous MRI and CT image from DICOM to nii

# Convert the T1 MRI image from .mgz to .nii
mri_convert -i ${SUBJECTS_DIR}/${SUBJECT_ID}_FS/mri/T1.mgz -o ${SUBJECTS_DIR}/${SUBJECT_ID}_FS/mri/T1.nii  -it mgz -ot nii

# Convert surface from mgz to nii
mri_convert -i ${SUBJECTS_DIR}/${SUBJECT_ID}_FS/mri/lh.ribbon.mgz -o ${SUBJECTS_DIR}/${SUBJECT_ID}_FS/mri/gray_left.nii -it mgz -ot nii
mri_convert -i ${SUBJECTS_DIR}/${SUBJECT_ID}_FS/mri/rh.ribbon.mgz -o ${SUBJECTS_DIR}/${SUBJECT_ID}_FS/mri/gray_right.nii -it mgz -ot nii
mri_convert -i ${SUBJECTS_DIR}/${SUBJECT_ID}_FS/mri/ribbon.mgz -o ${SUBJECTS_DIR}/${SUBJECT_ID}_FS/mri/t1_class.nii -it mgz -ot nii

# Convert the MRI from DICOM to nii
mri_convert -i ${SUBJECTS_DIR}/MRI/preop/preop_used/IM-0003-0001.dcm -o ${SUBJECTS_DIR}/Strip_${SUBJECT_ID}/anat_t1.nii -it dicom -ot nii

## mri_convert -i ${SUBJECTS_DIR}/leaddbs_mricon/DICOM/t2/IM-0019-0001-0001.dcm -o ${SUBJECTS_DIR}/leaddbs_mricon/anat_t2.nii -it dicom -ot nii

## mri_convert -i ${SUBJECTS_DIR}/${SUBJECT_ID}/Anatomy/MRI/postop/postop/IM-0002-0001-0001.dcm -o ${SUBJECTS_DIR}/${SUBJECT_ID}/Anatomy/postopmri.nii -it dicom -ot nii


# Convert the CT from DICOM to nii
mri_convert -i ${SUBJECTS_DIR}/CT/preop/preop_used/IM-0002-0001.dcm -o ${SUBJECTS_DIR}/Strip_${SUBJECT_ID}/preop_ct.nii -it dicom -ot nii

mri_convert -i ${SUBJECTS_DIR}/CT/postop/postop_used/IM-0002-0001.dcm -o ${SUBJECTS_DIR}/Strip_${SUBJECT_ID}/postop_ct.nii -it dicom -ot nii

## mri_convert -i ${SUBJECTS_DIR}/CT/intraop/intraop_used/IM-0001-0001.dcm -o ${SUBJECTS_DIR}/Strip_${SUBJECT_ID}/intraopct.nii -it dicom -ot nii

#If ct conversion runs into error, then do: 
#export FS_LOAD_DWI=0


# Check if nii files are in the subject folders
