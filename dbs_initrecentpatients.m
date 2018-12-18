function dbs_initrecentpatients(handles,patsub)

if ~exist('patsub','var')
    patsub='patients';
end

% make the user matlab path the root dirtectory
fprintf('Making %s the root directory for path variable storage.\n', userpath);
dbsroot = userpath;
%dbsroot=dbs_getroot;
if exist([dbsroot filesep 'dbs_recentpatients.mat']) == 2
    load([dbsroot filesep 'dbs_recentpatients.mat']);
else
    fullrpts={['No recent ',patsub,' found']};
end
save([dbsroot filesep 'dbs_recentpatients.mat'],'fullrpts');
dbs_updaterecentpatients(handles,patsub);
