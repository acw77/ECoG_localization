function dbs_updaterecentfsfolder(handles,fssub,nuchosenix)
if ~exist('fssub','var')
    fssub='freesurfer';
end
% make the user matlab path the root dirtectory
fprintf('Making %s the root directory for path variable storage.\n', userpath);
dbsroot = userpath;
%dbsroot=dbs_getroot;
load([dbsroot filesep 'dbs_recentfsfolders.mat']);
for i=1:length(fullrfs)
    [~,fullrfs{i}]=fileparts(fullrfs{i});
end
fullrfs=[{['Recent ',fssub,':']};fullrfs];
try
    set(handles.recentfs,'String',fullrfs);
catch
    return
end
if exist('nuchosenix','var')
   set(handles.recentfs,'Value',nuchosenix+1); 
end