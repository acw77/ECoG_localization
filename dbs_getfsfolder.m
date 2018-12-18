function dbs_getfsfolder(handles)


p='/'; % default use root
try
    p=pwd; % if possible use pwd instead (could not work if deployed)
end
try % finally use last patient parent dir if set.
    % make the user matlab path the root dirtectory
    fprintf('Making %s the root directory for path variable storage.\n', userpath);
    dbsroot = userpath;
    %dbsroot=dbs_getroot;
    if exist([dbsroot filesep 'dbs_recentfsfolders.mat'])
        load([dbsroot filesep 'dbs_recentfsfolders.mat']);
        p=fileparts(fullrpts{1});
    else
    end
end

uifsdir=dbs_uigetdir(p,'Please choose freesurfer folder(s)...');

if isempty(uifsdir)
    return
end

dbs_load_fs(handles,uifsdir);