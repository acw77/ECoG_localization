function [ptroot] = dbs_getrecentptfolder

% small function determining the location of the fs root directory.
% make the user matlab path the root dirtectory
fprintf('Making %s the root directory for path variable storage.\n', userpath);
root = userpath;
%dbsroot=dbs_getroot;
load(fullfile(root,'dbs_recentpatients.mat'))
ptroot = fullrpts{1};

end
