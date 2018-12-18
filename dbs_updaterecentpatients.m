function dbs_updaterecentpatients(handles,patsub,nuchosenix)
if ~exist('patsub','var')
    patsub='patients';
end
% make the user matlab path the root dirtectory
fprintf('Making %s the root directory for path variable storage.\n', userpath);
dbsroot = userpath;
%dbsroot=dbs_getroot;
load([dbsroot filesep 'dbs_recentpatients.mat']);
for i=1:length(fullrpts)
    [~,fullrpts{i}]=fileparts(fullrpts{i});
end
fullrpts=[{['Recent ',patsub,':']};fullrpts];
try
set(handles.recentpts,'String',fullrpts);
catch
    return
end
if exist('nuchosenix','var')
   set(handles.recentpts,'Value',nuchosenix+1); 
end