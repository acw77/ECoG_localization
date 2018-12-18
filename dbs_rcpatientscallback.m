function dbs_rcpatientscallback(handles)

if get(handles.recentpts,'Value')==1
    return
end
% make the user matlab path the root dirtectory
fprintf('Making %s the root directory for path variable storage.\n', userpath);
dbsroot = userpath;
%dbsroot=dbs_getroot;
load([dbsroot filesep 'dbs_recentpatients.mat']);
if iscell(fullrpts)
    fullrpts=fullrpts(get(handles.recentpts,'Value')-1);
end

if strcmp('No recent patients found',fullrpts)
    return
end


dbs_load_pts(handles,fullrpts);