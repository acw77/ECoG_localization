function [Elec] = Elec_equalization_Fcn(cfg)
% Elec_equialization_Fcn fills in the electrodes of a strip based on the
% poistions of the corner electrodes.
%
% cfg.Nrow - number of rows in the strip
% cfg.NperR - total number of electrodes per row 
% cfg.CorteElecLoc - locations of the electrodes in the corners of the
%   strip. These electrode should be marked following the electrode numbering
%   scheme. Allowed number of points are 2, 4 or 6. Provide 2 points to equalize 
%   a row, 6 points to equalize 3 rows, or 4 points to equalize 3 rows
%   assigning midpoint to middle row. 
% cfg.mask_indices - the hull as structure as loaded from hull.ma (Freesurfer)
%
% Optional arguments
% cfg.preElec - number of electrodes in strip before the first marked 
%     electrode. Defaults to zero. Note that these electrodes should be
%     included in NperR. 
% cfg.postElec - number of electrodes in strip after the last marked 
%     electrode. Defaults to zero. Note that these electrodes should be
%     included in NperR.
% cfg.hull_buffer - buffer to select cortical points. Default to 5. 
%     consider increasing in preElec or postElec are used. 
% 
% Returns a (Nrow*NperR)x3 matrix with electrode coordinates

%% specify electrode characteristics
% row number
Nrow = cfg.Nrow;
% number of contacts per row
NperR = cfg.NperR;

CortElecLoc = cfg.CortElecLoc;
mask_indices = cfg.mask_indices;

if isfield(cfg,'preElec')
  preElec = cfg.preElec;
else
  preElec = 0;
end

if isfield(cfg,'postElec')
  postElec = cfg.postElec;
else
  postElec = 0;
end

if isfield(cfg,'hull_buffer')
  hull_buffer = cfg.hull_buffer;
else
  hull_buffer = 5;
end


% If there are multiple electrodes, you want to deal one at a time
assert(ismember(numel(CortElecLoc),[2,4,6]),'Invalid number of coordinates to equalize.')
CortElecLoc0 = CortElecLoc; %CortElecLoc0 will save as original electrode coordinates
elec0 = reshape(cell2mat(CortElecLoc0),3,length(CortElecLoc0))';   % N X 3


%% Get the region of hull mask of interest
mask_indices0 = mask_indices; % mask_indices0 stores the entire hull mask
Max_Coordinates = max(elec0);
Min_Coordinates = min(elec0);
mask_indices(find(...
    mask_indices(:,1)>Max_Coordinates(1)+hull_buffer | mask_indices(:,1)<Min_Coordinates(1)-hull_buffer | ...
    mask_indices(:,2)>Max_Coordinates(2)+hull_buffer | mask_indices(:,2)<Min_Coordinates(2)-hull_buffer | ...
    mask_indices(:,3)>Max_Coordinates(3)+hull_buffer | mask_indices(:,3)<Min_Coordinates(3)-hull_buffer),:) = [];

%% define normal vector

% For one-row electrode which have only two marked dots, expand the line to
% a plane: find a line that passes the mid point and roughly perpendicular
% to the row
if Nrow == 1
    mid_pt = mean(elec0);
    distances_oi = pdist2(mid_pt,mask_indices);
    idx_oi = find(distances_oi < prctile(distances_oi,10) & distances_oi > prctile(distances_oi,5));
    Combinations = nchoosek(idx_oi,2);
    dist_of_comb = diag(pdist2(mask_indices(Combinations(:,1),:), mask_indices(Combinations(:,2),:)));
    far_pairs_idx = find(dist_of_comb > prctile(dist_of_comb,80));
    far_pairs = Combinations(far_pairs_idx,:);
    vectors_to_test = mask_indices(far_pairs(:,1),:) - mask_indices(far_pairs(:,2),:);
    my_vector = elec0(1,:) - elec0(2,:);
    my_vector_dup = repmat(my_vector, [size(vectors_to_test,1),1]);
    [~,pair_wanted_idx] = min(abs(dot(vectors_to_test, my_vector_dup, 2)));
    pair_wanted = far_pairs(pair_wanted_idx,:);
    dots_used = [mask_indices(pair_wanted,:); elec0];
else % for multiple-row electrode, dots_used = elec0
    dots_used = elec0;
end
% get dots_used
% % clearvars -except CortElecLoc CortElecLoc0 cortex elec0 mask_indices mask_indices0 NperR Nrow dots_used

% get 4 normal vectors based on the four points in dots_used

dot_combinations = nchoosek(1:4,3);

norm_vect_repo = [];
for i = 1:size(dot_combinations,1)
    dots_i = dots_used(dot_combinations(i,:),:);
    norm_vect_repo(i,:) = cross(dots_i(2,:) - dots_i(1,:), dots_i(3,:) - dots_i(1,:));
end

% make sure the four vectors are in the same direction
opp_idx = find(dot(norm_vect_repo, repmat(norm_vect_repo(1,:),[4,1]),2)<0);
norm_vect_repo(opp_idx,:) = -norm_vect_repo(opp_idx,:);

% then adjust molds of the vector to be the same (100)
molds = sqrt(dot(norm_vect_repo,norm_vect_repo,2 ));
norm_vect_repo_adj = norm_vect_repo * 100 ./ molds;

% average to get final norm vector
norm_final = mean(norm_vect_repo_adj);
mid_pt = mean(elec0,1);
another_dot = mid_pt + norm_final;

% % figure (1); 
% % plot3(mask_indices(:,1), mask_indices(:,2), mask_indices(:,3), 'g.')
% % axis equal
% % camlight('headlight','infinite');
% % hold on; plot3(elec0(:,1), elec0(:,2), elec0(:,3), 'b.', 'MarkerSize', 15)
% % plot3([another_dot(:,1), mid_pt(:,1)], [another_dot(:,2), mid_pt(:,2)], [another_dot(:,3), mid_pt(:,3)], 'm', 'MarkerSize', 15);

k = mask_indices;
k0 = k;

ElecEq_total = [];
for M = 1:(size(elec0,1)/2) % loop through existing rows
    norm2this = cross(another_dot - elec0(2*M - 1 , :), another_dot - elec0(2*M, :));
    
    clear v
    for i= 1:length(k)
        v(i) = dot(k(i,:)-another_dot,norm2this);
    end    
    
    hull_res = 3e2;
    %k(find(abs(v)>hull_res),:)=[];
    k(abs(v)>hull_res,:)=[];
    
    % find the line path of k
    % reordering the points for the spline fitting
    dist_mat = pdist2(k,k);
    [rr,~] = find(dist_mat == max(dist_mat(:)));
    temp = k(rr(1),:);
    k(rr(1),:) = []; 
    k = [temp;k];
    
    path=1;
    idx = [];
    for i=1:length(k)-1
        nodes = setdiff(1:length(k),path);
        [~,idx] = pdist2(k(nodes,:), k(path(end),:), 'euclidean', 'smallest', 1);
        path = cat(1, path, nodes(idx));
    end
    k=k(path,:);  
    
    %performing spline fitting
    s=cscvn(k');
    
    %looking for points in spline that correspond to marked 
    %electrode of the row 
    % %     hold on; fnplt(s, 'r', 2)

%     brk=1;
%     while ~(isequal(fnval(s, s.breaks(brk)), elec0(2*M - 1,:)') || ...
%       isequal(fnval(s, s.breaks(brk)), elec0(2*M,:)'))
%     brk = brk+1;
%     end
%     brk1 = brk; 
%     brk = brk+1;
%     while ~(isequal(fnval(s, s.breaks(brk)), elec0(2*M - 1,:)') || ...
%       isequal(fnval(s, s.breaks(brk)), elec0(2*M,:)'))
%     brk = brk+1;
%     end
%     brk2 = brk;
    brk1=find(vecnorm(fnval(s,s.breaks)' -  elec0(2*M - 1,:),2,2) < eps);
    brk2=find(vecnorm(fnval(s,s.breaks)' -  elec0(2*M,:),2,2) < eps); 
    
    %determining start and end points in spline according to preElec and
    %postElec, ie if there are electrodes beyond those marked
    
    %calculating arc along spline
    s_vec = linspace(s.breaks(1),s.breaks(end),500);
    darcds = vecnorm(diff(fnval(s,s_vec)'),2,2)' ./ diff(s_vec); %derivative of arc with respect to s
    arc = cumtrapz(darcds)*mean(diff(s_vec)); 
    arc=[arc,arc(end)+mean(diff(arc))];

    %trasforming from s coord to arc
    arc1 = interp1(s_vec,arc,s.breaks(brk1));
    arc2 = interp1(s_vec,arc,s.breaks(brk2));
    
    %distance between elec along spline in arc units
    arc_elec_dist = (arc2 - arc1)/(NperR - postElec - preElec);
    arc_start = arc1 - arc_elec_dist .* preElec;
    arc_end = arc2 + arc_elec_dist .* postElec;
    
    %calculating elec position in s coords  
    s_elec_loc = interp1(arc,s_vec,linspace(arc_start,arc_end,NperR));
    ElecEq = fnval(s,s_elec_loc)';
    
    assert(~any(any(ismissing(ElecEq))),'Extrapolation failed. Increase hull_buffer');
    
    %updating marked electrode position for middle row
    if preElec > 0
      elec0(2*M - 1 , :) = ElecEq(1,:);
    end
    if postElec > 0
      elec0(2*M, :) = ElecEq(end,:);
    end
    
    
%   elec_line=fnval(s, linspace(s_start,s_end,400))';
%   elec_line=fnval(s, linspace(s.breaks(brk1),s.breaks(brk2),400))';
%   [ElecEq, Temp_distances] = SubCurve( elec_line, NperR);
    
%     figure();
%     hold on; 
%     plot3(ElecEq(:,1), ElecEq(:,2), ElecEq(:,3), 'r.', 'MarkerSize', 15);
%     plot3(ElecEq(1,1), ElecEq(1,2), ElecEq(1,3), 'y.', 'MarkerSize', 15);
%     disp('First contact of each row is shown in yellow; check the order');
    
    ElecEq_total = [ElecEq_total; ElecEq];
    k = k0;
end

if Nrow == 3 % for 3-row electrode, needs to add the middle row
    
    middle_start = mean([elec0(1,:);elec0(3,:)]);
    middle_end = mean([elec0(2,:);elec0(4,:)]);
    
    norm2this = cross(another_dot - middle_start, another_dot - middle_end);
    
    clear v
    for i= 1:length(k)
        v(i) = dot(k(i,:)-another_dot,norm2this);
    end  
    
    hull_res_here = hull_res;
    
    qualified_s = 0;
    safety_counter = 0;
    while qualified_s == 0
        k = k0;
        hull_res_here = hull_res_here + 100;
        k(find( abs(v)>hull_res_here),:)=[]; 

        % find the line path of k

        dist_mat = pdist2(k,k);
        [rr,cc] = find(dist_mat == max(dist_mat(:)));
        temp = k(rr(1),:);
        k(rr(1),:) = []; k = [temp;k];

        path=1;
        idx = [];
        for i=1:length(k)-1
            nodes = setdiff(1:length(k),path);
            [~,idx] = pdist2(k(nodes,:), k(path(end),:), 'euclidean', 'smallest', 1);
            path = cat(1, path, nodes(idx));
        end
        k=k(path,:);  

        s=cscvn(k');
        
        % test if s is qualified
        s_start = fnval(s,s.breaks(1));
        s_end = fnval(s,s.breaks(end));
        
        side1 = elec0(1:2,:);
        side2 = elec0(3:4,:);
  
% %         clearvars determine
        [~,idx] = pdist2(side1, s_start','euclidean', 'smallest', 1);
        determine(1) = dot(s_start' - side1(idx,:), side1(setdiff(1:2,idx),:) - side1(idx,:));
        
        [~,idx] = pdist2(side2, s_start','euclidean', 'smallest', 1);
        determine(2) = dot(s_start' - side2(idx,:), side2(setdiff(1:2,idx),:) - side2(idx,:));
        
        [~,idx] = pdist2(side1, s_end','euclidean', 'smallest', 1);
        determine(3) = dot(s_end' - side1(idx,:), side1(setdiff(1:2,idx),:) - side1(idx,:));
        
        [~,idx] = pdist2(side2, s_end','euclidean', 'smallest', 1);
        determine(4) = dot(s_end' - side2(idx,:), side2(setdiff(1:2,idx),:) - side2(idx,:));
        
        if any(determine>=0)
            qualified_s = 0;
        else
            qualified_s = 1;
        end
        
        if safety_counter > 100
          error('max iterations reached');
        else
          safety_counter = safety_counter + 1;
        end
        
    end  % get qualified s
    
% %     hold on; fnplt(s, 'r', 2)
    
    elec_line=fnval(s, linspace(s.breaks(1),s.breaks(end),400))';
    
    [~,idx_end1] = min(pdist2(elec_line, elec0(1,:)) + pdist2(elec_line, elec0(3,:)));
    
    [~,idx_end2] = min(pdist2(elec_line, elec0(2,:)) + pdist2(elec_line, elec0(4,:)));
    
    [ElecEq, Temp_distances] = SubCurve( elec_line(min(idx_end1,idx_end2):max(idx_end1, idx_end2),:), NperR);
    if idx_end1 > idx_end2
      ElecEq = flipud(ElecEq);
    end
    
    
% %     hold on; plot3(ElecEq(:,1), ElecEq(:,2), ElecEq(:,3), 'r.', 'MarkerSize', 15);
% %     plot3(ElecEq(1,1), ElecEq(1,2), ElecEq(1,3), 'y.', 'MarkerSize', 15);
% %     disp('First contact of each row is shown in yellow; check the order');
    
    ElecEq_total = [ElecEq_total(1:NperR,:); ElecEq; ElecEq_total((end - NperR + 1):end,:)];
end

Elec = mat2cell(ElecEq_total, ones(NperR * Nrow,1), 3)';

