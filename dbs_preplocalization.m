function dbs_preplocalization(options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare_for_Localization                                                %
% Michael Randazzo 06/29/2015                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ptdirectory = char(options.uipatdirs);
fsdirectory = char(options.uifsdir);
cd(fsdirectory)
% Check for Electrode_Locations Folder
if ~exist([char(options.uifsdir) '/Electrode_Locations'])
    mkdir([char(options.uifsdir),'/Electrode_Locations']);
% if ~exist([directory '\Electrode_Locations'])  % Windows
%     mkdir([directory,'\Electrode_Locations']); % Windows
end

%%%%%%%%%%%% Cortex %%%%%%%%%%%%%%
[cortex.vert_lh,cortex.tri_lh]= read_surf(fullfile(fsdirectory,'surf/lh.pial')); % Reading left side pial surface
[cortex.vert_rh,cortex.tri_rh]= read_surf(fullfile(fsdirectory,'surf/rh.pial')); % Reading right side pial surface

% Generating entire cortex
cortex.vert = [cortex.vert_lh; cortex.vert_rh]; % Combining both hemispheres
cortex.tri = [cortex.tri_lh; (cortex.tri_rh + length(cortex.vert_lh))]; % Combining faces (Have to add to number of faces)

cortex.tri=cortex.tri+1; % freesurfer starts at 0 for indexing

% Reading in MRI parameters
% WJL edit 06/12/2017 
% f=MRIread(fullfile(fsdirectory,'mri/T1.mgz'));
f=MRIread(fullfile(fsdirectory,'mri/T1.nii'));

% Translating into the appropriate space
for k=1:size(cortex.vert,1)
    a=f.vox2ras/f.tkrvox2ras*[cortex.vert(k,:) 1]';
    cortex.vert(k,:)=a(1:3)';
end

save(fullfile(fsdirectory,'cortex_indiv.mat'),'cortex');

% Skull
% Get Mesh
[~,fsfiles] = dbs_subdir(fsdirectory);
% [patsub , ptfiles] = dbs_subdir(ptdirectory);
% ptfiles(~cellfun(@isempty,strfind(ptfiles,'.stl')))
try
    fname = char(fsfiles(~cellfun(@isempty,strfind(fsfiles,'.obj'))));
    disp('running    dbs_displayobj.......')
    [skull.vert,skull.tri] = dbs_displayobj(fname);
catch
    [fname,~,extension] = dbs_uigetfile(ptdirectory,'Choose Skull Mesh .obj');
    [skull.vert,skull.tri] = dbs_displayobj([fname extension]);
end

if ~exist(fullfile(fsdirectory,'ct_reg'))
    mkdir(fullfile(fsdirectory,'ct_reg'))
else
end
disp(['saving    ' fullfile(fsdirectory,'ct_reg','skull.mat') ' .......'])
save(fullfile(fsdirectory,'ct_reg','skull.mat'),'skull')

% Option to prevent overwriting:
% if ~exist(fullfile(directory,fname))
%     try
%         [skull.vert,skull.tri] = dbs_displayobj(fname);
%         save(fullfile(directory,'skull.mat'),'skull')
%     catch
%         cd(directory)
%         [skull.vert,skull.tri] = dbs_displayobj(fname);
%         save(fullfile(directory,'skull.mat'),'skull')
%     end
% else
% end

% Create cortical hull
grayfilename = [char(options.uifsdir) '/mri/t1_class.nii'];
outputdir='hull';
disp('running dbs_gethull.......')
[~,~] = dbs_gethull(grayfilename,outputdir,3,21,.3);
disp('** Process Done')
end

