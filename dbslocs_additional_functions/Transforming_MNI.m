% This script transforms electrode coordinates and brain surface
% Template script. Keep one for each subject

% coordinates from individual space to MNI space
% Export brainstorm matlab file to matlab as sMRI before running this script

% add brainstorm to path
addpath(genpath('E:\MATLAB\brainstorm3'));
% define path variables
SUBJECT='DBS3021';
PATH_DATA='\\136.142.16.9\Nexus\DBS';
PATH_FREESURFER = [PATH_DATA filesep SUBJECT '\Anatomy\FreeSurfer\preop'];

% load cortical surface and electrode locations
load([PATH_FREESURFER filesep 'cortex_indiv'], 'cortex');
load([PATH_FREESURFER filesep 'Electrode_Locations\CortElecLocL_eq.mat'], 'CortElecLoc');

% load T1 image, freesurfer function MRIread.m required
f = MRIread([PATH_FREESURFER filesep 'mri\T1.nii']);

% load in sMRI if already there 
load([PATH_FREESURFER filesep 'sMRI_' SUBJECT])

% cell2mat of CortElecLoc
elec = reshape(cell2mat(CortElecLoc),3,length(CortElecLoc))';
elec0 = elec;

cortex0 = cortex;

% vox2ras transformation matrix
vox2ras = [1,0,0,-128;0,1,0,-128;0,0,1,-128;0,0,0,1];

%% Convert  electrode from freesurfer coordinate (tkr-ras) to voxel
clear elec_vox
for k=1:size(elec,1)
    a=inv(vox2ras) * f.tkrvox2ras/f.vox2ras*[elec(k,:) 1]';
    elec_vox(k,:)=a(1:3)';
end

%transform from voxel to mni, brainstorm function cs_convert.m required
elec_mni_m = cs_convert(sMRI,'voxel','mni',elec_vox);
elec_mni_mm = elec_mni_m * 1000; %  in milimeters
% at this stage, add mni coords to the electrode table (Processed Data/Sync/annot)

% store in CortElecLoc_MNI
CortElecLoc_MNI = mat2cell(elec_mni_mm, ones(length(CortElecLoc),1), 3)';
save([PATH_FREESURFER filesep 'Electrode_Locations\CortElecLoc_MNI.mat'],'CortElecLoc_MNI')

%% convert individual brain to mni space
for k=1:size(cortex.vert,1)
    a=inv(vox2ras) * f.tkrvox2ras/f.vox2ras*[cortex.vert(k,:) 1]';
    cortex.vert(k,:)=a(1:3)';
end

% transform from voxel to mni, brainstorm function cs_convert.m required
cortex.vert = cs_convert(sMRI,'voxel','mni',cortex.vert);
cortex.vert = cortex.vert * 1000; %  in milimeters

% store in cortex_indov_mni
save([PATH_FREESURFER filesep 'cortex_indiv_mni.mat'],'cortex');

%% save sMRI as sMRI_DBSXXXX
save([PATH_FREESURFER filesep 'sMRI_' SUBJECT '.mat'], 'sMRI');
%% Plot

% load in MNI brain
load('E:\MATLAB\lead_v2.1.5\templates\space\MNI_ICBM_2009b_NLIN_ASYM\cortex\CortexHiRes.mat', 'Vertices', 'Faces');

figure; % MNI electrode plotted with MNI brain
Hp = patch('vertices', Vertices,'faces', Faces,...
    'facecolor',[.50 .50 .50],'edgecolor','none',...
    'facelighting', 'gouraud', 'specularstrength', .50);
camlight('headlight','infinite');
axis off; axis equal
alpha 0.5
hold on; plot3(elec_mni_mm(:,1), elec_mni_mm(:,2), elec_mni_mm(:,3), 'b.', 'MarkerSize', 15)

figure; % electrode and brain in individual space
Hp = patch('vertices',cortex0.vert,'faces',cortex0.tri(:,[1 3 2]),...
    'facecolor',[.85 .50 .50],'edgecolor','none',...
    'facelighting', 'gouraud', 'specularstrength', .50);
camlight('headlight','infinite');
axis off; axis equal
alpha 0.5
hold on; plot3(elec0(:,1), elec0(:,2), elec0(:,3), 'b.', 'MarkerSize', 15)

figure; % transformed individual brain vs. MNI brain

Hp1 = patch('vertices', Vertices,'faces', Faces,...
    'facecolor',[.50 .50 .50],'edgecolor','none',...
    'facelighting', 'gouraud', 'specularstrength', .50);
camlight('headlight','infinite');
axis off; axis equal
alpha 0.5

hold on;
Hp2 = patch('vertices',cortex.vert,'faces',cortex.tri(:,[1 3 2]),...
    'facecolor',[.85 .50 .50],'edgecolor','none',...
    'facelighting', 'gouraud', 'specularstrength', .50);
camlight('headlight','infinite');