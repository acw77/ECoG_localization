
load('E:\MATLAB\brainstorm_db\DBS2000_series\anat\@default_subject\tess_cortex_pial_high.mat', 'Atlas')
load('Z:\Commits\DBS_speech\CortElecLoc_MNI.mat', 'CortElecLoc_MNI')
elec = reshape(cell2mat(CortElecLoc_MNI),3,length(CortElecLoc_MNI))';


% Using select cursor, click on central point from a perpendicular view
pos1=get(gca,'CurrentPoint');
hold on; plot3(pos1(:,1), pos1(:,2), pos1(:,3), 'm');
hold on; plot3(pos1(:,1), pos1(:,2), pos1(:,3), 'm.', 'markersize', 25);


(pos1(1,:)-pos1(2,:))/norm(pos1(1,:)-pos1(2,:))

figure;
hold on; plot3(elec(:,1), elec(:,2), elec(:,3), 'r.', 'MarkerSize', 25)

elec1 = elec + 5*(pos1(1,:)-pos1(2,:))/norm(pos1(1,:)-pos1(2,:));

hold on; plot3(elec1(:,1), elec1(:,2), elec1(:,3), 'r.', 'MarkerSize', 25)

figure;
Hp = patch('Vertices',BS1.Vertices,'Faces',BS1.Faces,...
    'facecolor',[1 .5 .5],'edgecolor','none',...
    'facelighting', 'gouraud', 'specularstrength', .50);
camlight('headlight','infinite');
axis off; axis equal
alpha 1


idxAtlas = 3;
idxScouts = AnatomyCatIdx;

VColor = repmat([1 1 1], [length(BS1.Vertices), 1]);
for i=1:length(idxScouts)
VColor(Atlas(idxAtlas).Scouts(idxScouts(i)).Vertices,:) = ...
    repmat(Atlas(idxAtlas).Scouts(idxScouts(i)).Color, [length(Atlas(idxAtlas).Scouts(idxScouts(i)).Vertices), 1]);
end

figure;
Hp = patch('Vertices',BS1.Vertices,'Faces',BS1.Faces,'facecolor', 'interp',...
    'FaceVertexCData',VColor,'edgecolor','none',...
    'facelighting', 'gouraud', 'specularstrength', .50);
camlight('headlight','infinite');
axis off; axis equal
alpha 1


AnatomyCat = {};
AnatomyCatIdx = [];
for el=1:length(CortElecLoc_MNI)
    [~, vmin] = min(pdist2(CortElecLoc_MNI{el}, BS1.Vertices));
    AnatomyCatIdx(el) = find(arrayfun(@(x) ismember( vmin, x.Vertices), Atlas(idxAtlas).Scouts));
    AnatomyCat{el} = Atlas(idxAtlas).Scouts( AnatomyCatIdx(el) ).Label;
end

for e=1:length(elec)
   hold on; plot3(elec1(e,1), elec1(e,2), elec1(e,3), 'o', 'color', 'g', 'MarkerSize', 10)
   hold on; text(elec1(e,1), elec1(e,2), elec1(e,3), num2str(e))
end
