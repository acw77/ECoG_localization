function dbs_getpatients(handles)


p='/'; % default use root
try
    p=pwd; % if possible use pwd instead (could not work if deployed)
end
try % finally use last patient parent dir if set.
    % make the user matlab path the root dirtectory
    fprintf('Making %s the root directory for path variable storage.\n', userpath);
    dbsroot = userpath;
    %dbsroot=dbs_getroot;
    if exist([dbsroot filesep 'dbs_recentpatients.mat'])
        load([dbsroot filesep 'dbs_recentpatients.mat']);
        p=fileparts(fullrpts{1});
    else
    end
end

uipatdir=dbs_uigetdir(p,'Please choose patient folder(s)...');

if isempty(uipatdir)
    return
end

dbs_load_pts(handles,uipatdir);