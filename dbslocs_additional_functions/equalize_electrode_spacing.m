%Equalize Cortical Electrode Spacing
%Last Updated 08/04/2017 Steven Lo

%Load Electrodes

%%
% Check contact spacing
for i=1:(length(CortElecLoc)-1)
    fprintf('Dist %d-%d: %4.2f\n', i, i+1, pdist2(CortElecLoc{i}, CortElecLoc{i+1}));
end



% Set electrode numbers for row
firstelec = 3;
lastelec = 4;
fprintf('Dist %d-%d: %4.2f\n', 1, 22, pdist2(CortElecLoc{1}, CortElecLoc{14}));
fprintf('Dist %d-%d: %4.2f\n', 15, 28, pdist2(CortElecLoc{15}, CortElecLoc{28}));

%%
% Load hull.mat
figure; plot3(mask_indices(:,1), mask_indices(:,2), mask_indices(:,3), 'g.')
axis equal
camlight('headlight','infinite');
% use brush tool to select a segment of the hull --> k

%%
% Plot segmented hull
l = k;
l0 = l;
figure; plot3(l(:,1), l(:,2), l(:,3), 'g.')
axis equal
camlight('headlight','infinite');

% Plots electrodes
CortElecLoc0 = CortElecLoc;
elec = reshape(cell2mat(CortElecLoc0),3,length(CortElecLoc0))';
hold on; plot3(elec(:,1), elec(:,2), elec(:,3), 'r.', 'MarkerSize', 25)

%%
% Using select cursor, click on central point from a perpendicular view
pos1=get(gca,'CurrentPoint');
hold on; plot3(pos1(:,1), pos1(:,2), pos1(:,3), 'm');

%%
% Plots normal vector and line
norm = cross(pos1(2,:)-elec(firstelec,:), pos1(2,:)-elec(lastelec,:));

clear v
for i= 1:length(l)
    v(i) = dot(l(i,:)-pos1(2,:),norm);
end

%hull_res = 1.0e2;
hull_res = 2.5e2;
%hull_res = 3e2;
l(find( abs(v)>hull_res),:)=[];


l = [elec(firstelec,:); l; elec(lastelec,:)];

%l(arrayfun(@(x) max(arrayfun(@(y) max(isequal(l(x,:), temp(y,:))), 1:length(temp))), 1:length(l)),:)=[];

%l = l(end:-1:1,:);

path=1;
for i=1:length(l)-1
    nodes = setdiff(1:length(l),path);
    [~,idx] = pdist2(l(nodes,:), l(path(end),:), 'euclidean', 'smallest', 1);
    path = cat(1, path, nodes(idx));
end
l=l(path,:);

epsilon = 1e-10;

s=cscvn(l');
hold on; fnplt(s, 'b', 2)
brk=1;
while ~(pdist2(fnval(s, s.breaks(brk))' - elec(firstelec,:), [0,0,0])<epsilon || ...
        pdist2(fnval(s, s.breaks(brk))' - elec(lastelec,:), [0,0,0])<epsilon)
brk = brk+1;
end
brk1 = brk; 
brk = brk+1;
while ~(pdist2(fnval(s, s.breaks(brk))' - elec(firstelec,:), [0,0,0])<epsilon || ...
        pdist2(fnval(s, s.breaks(brk))' - elec(lastelec,:), [0,0,0])<epsilon)
brk = brk+1;
end
brk2 = brk;

elec_line=fnval(s, linspace(s.breaks(brk1),s.breaks(brk2),400))';

% elec_num = size(elec,1);
elec_num = 21;
[ElecEq, distances] = SubCurve( elec_line, elec_num );
hold on; plot3(ElecEq(:,1), ElecEq(:,2), ElecEq(:,3), 'b.', 'MarkerSize', 25)
hold on; plot3(ElecEq(1,1), ElecEq(1,2), ElecEq(1,3), 'y.', 'MarkerSize', 25)

% continue electrodes
elec_num = 21;
mean_dist = mean(distances);
total_dist = (elec_num-1)*mean_dist;
while (pdist2(fnval(s, s.breaks(brk))',elec(firstelec,:))<total_dist)
brk = brk+1;
end
brk3 = brk;
elec_line=fnval(s, linspace(s.breaks(brk1),s.breaks(brk3),400))';
last_pt = 1;
while arclength3(elec_line(1:(last_pt+1),:)) < total_dist
    last_pt = last_pt+1;
    %fprintf('%4.2f\n', arclength3(elec_line(1:last_pt,:)))
end
elec_line = elec_line(1:last_pt,:);
[ElecEq, distances] = SubCurve( elec_line, elec_num );
hold on; plot3(ElecEq(:,1), ElecEq(:,2), ElecEq(:,3), 'r.', 'MarkerSize', 25)
hold on; plot3(ElecEq(1,1), ElecEq(1,2), ElecEq(1,3), 'y.', 'MarkerSize', 25)

elec = [ElecEq1; ElecEq2; ElecEq3];
CortElecLoc = mat2cell(elec, ones(63,1), 3)';
save('CortElecLocR_eq', 'CortElecLoc', 'CortElecLoc0')

CortElecLoc = mat2cell(ElecEq, ones(elec_num,1), 3)';
% check if electrode contact order got switched; if so, reverse order 
if ~isequal(CortElecLoc0{firstelec}, CortElecLoc{1})
    fprintf('Electrode order switched.\n')
    CortElecLoc = CortElecLoc(elec_num:-1:1);
end

%%
% Plot first electrode, check if is correct
hold on; plot3(elec(1,1), elec(1,2), elec(1,3), 'y.', 'MarkerSize', 25)

% Renumber if reversed
CortElecLoc = CortElecLoc(elec_num:-1:1);

% Replot to check
hold on; plot3(elec(1,1), elec(1,2), elec(1,3), 'm.', 'MarkerSize', 25)
%%
%FOR MULTI ROW STRIPS
CortElecLoc1=CortElecLoc;
distances1=distances;

CortElecLoc2=CortElecLoc;
distances2=distances;

CortElecLoc3=CortElecLoc;
disttances3=distances;

CortElecLoc6=CortElecLoc;
distances6=distances;

%%
% Plot all electrodes, check if correct
hold on; plot3(ElecEq(:,1), ElecEq(:,2), ElecEq(:,3), 'b.', 'MarkerSize', 25)
hold on; plot3(elec(1,1), elec(1,2), elec(1,3), 'y.', 'MarkerSize', 25)
hold on; plot3(elec(22,1), elec(22,2), elec(22,3), 'c.', 'MarkerSize', 25)
hold on; plot3(elec(43,1), elec(43,2), elec(43,3), 'm.', 'MarkerSize', 25)

%%
% Concatenates separate rows into one variable
%28
CortElecLoc=cat(2, CortElecLoc1, CortElecLoc2);

%54/63
CortElecLoc=cat(2, cat(2, CortElecLoc1, CortElecLoc2), CortElecLoc3);

%54/63 + 6
CortElecLoc=cat(2, cat(2, CortElecLoc1, CortElecLoc2), cat(2, CortElecLoc3, CortElecLoc6));

%%
% Saves equalized electrode locations with electrod spacing
save('CortElecLocLF485_eq', 'CortElecLoc0', 'CortElecLoc', 'distances');

save('CortElecLocLF615_eq', 'CortElecLoc0', 'CortElecLoc', 'distances', 'distances6');