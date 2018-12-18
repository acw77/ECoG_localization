function [fsroot] = dbs_getrecentfsfolder

% small function determining the location of the fs root directory.
% make the user matlab path the root dirtectory
fprintf('Making %s the root directory for path variable storage.\n', userpath);
root = userpath;
%dbsroot=dbs_getroot;
load(fullfile(root,'dbs_recentfsfolders.mat'))
fsroot = fullrfs{1};

end