function dbs_rcfscallback(handles)

if get(handles.recentfs,'Value')==1
    return
end
% make the user matlab path the root dirtectory
fprintf('Making %s the root directory for path variable storage.\n', userpath);
dbsroot = userpath;
%dbsroot=dbs_getroot;
load([dbsroot filesep 'dbs_recentfsfolders.mat']);
if iscell(fullrfs)
    fullrfs=fullrfs(get(handles.recentfs,'Value')-1);
end

if strcmp('No recent patients found',fullrfs)
   return 
end


dbs_load_fs(handles,fullrfs);