% load in hull.mat, cortex_indiv.mat and CortElecLoc.mat

%% specify electrode characteristics
% row number
Nrow = 3;
% number of contacts per row
NperR = 21;

% If there are multiple electrodes, you want to deal one at a time

CortElecLoc0 = CortElecLoc; %CortElecLoc0 will save as original electrode coordinates
elec0 = reshape(cell2mat(CortElecLoc0),3,length(CortElecLoc0))';   % N X 3


%% Get the region of hull mask of interest
mask_indices0 = mask_indices; % mask_indices0 stores the entire hull mask
Max_Coordinates = max(elec0);
Min_Coordinates = min(elec0);
mask_indices(find(mask_indices(:,1)>Max_Coordinates(1) +5 | mask_indices(:,1)<Min_Coordinates(1) -5 | ...
    mask_indices(:,2)>Max_Coordinates(2)+5 | mask_indices(:,2) < Min_Coordinates(2)-5 | ...
    mask_indices(:,3)>Max_Coordinates(3)+5 | mask_indices(:,3) < Min_Coordinates(3)-5),:) = [];
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
clearvars -except CortElecLoc CortElecLoc0 cortex elec0 mask_indices mask_indices0 NperR Nrow dots_used

% get 4 normal vectors based on the four points in dots_used

dot_combinations = nchoosek(1:4,3);

norm_vect_repo = [];
for i = 1:size(dot_combinations,1)
    
    dots_i = dots_used(dot_combinations(i,:),:);
    
    norm_vect_repo(i,:) = cross(dots_i(2,:) - dots_i(1,:), dots_i(3,:) - dots_i(1,:));
end

% make sure the four vectors are in the same direction
opp_idx = find(dot(norm_vect_repo, repmat(norm_vect_repo(1,:),[4,1]),2)<0);
norm_vect_repo(opp_idx,:) = -norm_vect_repo(opp_idx,:)

% then adjust molds of the vector to be the same (100)
molds = sqrt(dot(norm_vect_repo,norm_vect_repo,2 ));

norm_vect_repo_adj = norm_vect_repo * 100 ./ molds;

% average to get final norm vector
norm_final = mean(norm_vect_repo_adj);

mid_pt = mean(elec0,1);

another_dot = mid_pt + norm_final;

figure (1); 
plot3(mask_indices(:,1), mask_indices(:,2), mask_indices(:,3), 'g.')
axis equal
camlight('headlight','infinite');
hold on; plot3(elec0(:,1), elec0(:,2), elec0(:,3), 'b.', 'MarkerSize', 15)
plot3([another_dot(:,1), mid_pt(:,1)], [another_dot(:,2), mid_pt(:,2)], [another_dot(:,3), mid_pt(:,3)], 'm', 'MarkerSize', 15);

k = mask_indices;
k0 = k;

ElecEq_total = [];
for M = 1: size(elec0,1)/2 % loop through existing rows
    norm2this = cross(another_dot - elec0(2*M - 1 , :), another_dot - elec0(2*M, :));
    
    clear v
    for i= 1:length(k)
        v(i) = dot(k(i,:)-another_dot,norm2this);
    end    
    
%     hull_res = 1.0e2;
%     hull_res = 2.5e2;
    hull_res = 3e2;
    k(find( abs(v)>hull_res),:)=[];
    
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
    hold on; fnplt(s, 'r', 2)
    brk=1;
    while ~(isequal(fnval(s, s.breaks(brk)), elec0(2*M - 1,:)') || ...
    isequal(fnval(s, s.breaks(brk)), elec0(2*M,:)'))
    brk = brk+1;
    end
    brk1 = brk; 
    brk = brk+1;
    while ~(isequal(fnval(s, s.breaks(brk)), elec0(2*M - 1,:)') || ...
    isequal(fnval(s, s.breaks(brk)), elec0(2*M,:)'))
    brk = brk+1;
    end
    brk2 = brk;

    elec_line=fnval(s, linspace(s.breaks(brk1),s.breaks(brk2),400))';


    [ElecEq, Temp_distances] = SubCurve( elec_line, NperR);
    figure (1);
    hold on; plot3(ElecEq(:,1), ElecEq(:,2), ElecEq(:,3), 'r.', 'MarkerSize', 15);
    plot3(ElecEq(1,1), ElecEq(1,2), ElecEq(1,3), 'y.', 'MarkerSize', 15);
    disp('First contact of each row is shown in yellow; check the order');
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
        
        [~,idx] = pdist2(side1, s_start','euclidean', 'smallest', 1);
        
        clearvars determine
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
    end  % get qualified s
    
    hold on; fnplt(s, 'r', 2)
    
    elec_line=fnval(s, linspace(s.breaks(1),s.breaks(end),400))';
    
    [~,idx_end1] = min(pdist2(elec_line, elec0(1,:)) + pdist2(elec_line, elec0(3,:)))
    
    [~,idx_end2] = min(pdist2(elec_line, elec0(2,:)) + pdist2(elec_line, elec0(4,:)))
    
    [ElecEq, Temp_distances] = SubCurve( elec_line(min(idx_end1,idx_end2):max(idx_end1, idx_end2),:), NperR);
    
    hold on; plot3(ElecEq(:,1), ElecEq(:,2), ElecEq(:,3), 'r.', 'MarkerSize', 15);
    plot3(ElecEq(1,1), ElecEq(1,2), ElecEq(1,3), 'y.', 'MarkerSize', 15);
    disp('First contact of each row is shown in yellow; check the order');
    
    ElecEq_total = [ElecEq_total(1:NperR,:); ElecEq; ElecEq_total(end - NperR + 1:end,:)];
end
CortElecLoc = mat2cell(ElecEq_total, ones(NperR * Nrow,1), 3)';

disp('Equalized electrode location is saved as CortElecLoc_eq.mat');
save('./Electrode_Locations/CortElecLocL_eq','CortElecLoc');