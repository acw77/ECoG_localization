% define path variables
SUBJECT='DBS3021';
PATH_DATA='\\136.142.16.9\Nexus\DBS';
PATH_FREESURFER = [PATH_DATA filesep SUBJECT '\Anatomy\FreeSurfer\preop'];
PATH_SYNC = [PATH_DATA filesep SUBJECT '\Preprocessed Data\Sync'];

% load cortical surface and electrode locations
load([PATH_FREESURFER filesep 'cortex_indiv'], 'cortex');
load([PATH_FREESURFER filesep 'Electrode_Locations\CortElecLocL_eq.mat'], 'CortElecLoc');

% load annot tables
session    = bml_annot_read([PATH_SYNC filesep 'annot' filesep SUBJECT '_session.txt']);
electrode    = bml_annot_read([PATH_SYNC filesep 'annot' filesep SUBJECT '_electrode.txt']);
mercoords    = ...
    readtable([PATH_DATA filesep ...
               SUBJECT filesep ...
               'Anatomy' filesep ...
               'leaddbs_' SUBJECT  filesep ...
               'annot' filesep ...
               SUBJECT '_MER_coords.txt']);

           
           
% add columns to electrode table
electrode.mni_x = ones(size(electrode.id))*NaN;
electrode.mni_y = ones(size(electrode.id))*NaN;
electrode.mni_z = ones(size(electrode.id))*NaN;
           

% write new electrode file
writetable(electrode, [PATH_SYNC filesep 'annot' filesep SUBJECT '_electrode.txt'], 'Delimiter', '\t'); 
           
