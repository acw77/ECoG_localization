function dbs_initrecentfsfolders(handles,fssub)

if ~exist('fssub','var')
    fssub='freesurfer';
end

% make the user matlab path the root dirtectory
fprintf('Making %s the root directory for path variable storage.\n', userpath);
dbsroot = userpath;
%dbsroot=dbs_getroot;
if exist([dbsroot filesep 'dbs_recentfsfolders.mat']) == 2
    load([dbsroot filesep 'dbs_recentfsfolders.mat']);
else
    fullrfs={['No recent ',fssub,' found']};
end
save([dbsroot filesep 'dbs_recentfsfolders.mat'],'fullrfs');
dbs_updaterecentfsfolder(handles,fssub);