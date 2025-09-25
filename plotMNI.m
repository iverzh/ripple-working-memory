


addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/spm12'))

% clear
close all
matExportFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/matFiles';

dataDirectory = '/space/seh10/6/halgdev/projects/iverzh/data/Sternberg/preprocess/OrigUpload';
flst = dir(fullfile(dataDirectory, '*LFP_micro*'));
flst = {flst.name}';
flst = strrep(flst, '_LFP_micro.mat', '');
subj_list_full = flst;






% Using SPM
% template = spm_vol('/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/MNI152_T1_1mm.nii');
template = spm_vol('/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/CIT168_Labeling_Templates/CV_639317_T2w_avg.nii');
[volData, ~] = spm_read_vols(template);

templateMask = spm_vol('/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/icbm_avg_152_t1_tal_nlin_symmetric_VI_mask.nii');
[volMask, ~] = spm_read_vols(templateMask);
mask = (volMask(6:end-6, 6:end-6, 6:end-6));

scale_factor = 1; % reduces resolution to 25% of original (4mm voxels)
volLowRes = imresize3(volData, scale_factor, 'linear');
maskLowRes = imresize3(mask, scale_factor, 'linear');



%%
% Threshold the data to isolate the brain (adjust threshold as necessary)
thresh = 6000; % Adjust to isolate the cortex clearly
binaryCortex = volLowRes > 2.0 & volLowRes < 3.5; %& logical(maskLowRes);
% binaryCortex = volLowRes > 1.0 & volLowRes < 2; %& logical(maskLowRes);

% Example: Remove objects smaller than 50 pixels

% Smooth slightly to create a more natural cortical surface
binaryCortex = smooth3(binaryCortex, 'gaussian', 5);


scroll3DSlice(binaryCortex)


binaryCortex = double(binaryCortex > 1-5e-2);
minPixels = 1000;  % Set your threshold here
binaryCortex = bwareaopen(binaryCortex, minPixels);
scroll3DSlice(binaryCortex)

%%
regions = {'OFC', 'ACC', 'SMA', 'AMY', 'HIP'} ;  
regionColors =  brewermap(12, 'Dark2');

% Generate isosurface
fv = isosurface(binaryCortex, 0.5);

% Simple direct affine transformation (recommended)
vertices_vox = (fv.vertices - 1) / scale_factor; % undo scaling and indexing
vertices_mni = (template.mat * [vertices_vox ones(size(vertices_vox,1),1)]')';
vertices_mni = vertices_mni(:,1:3);
vertices_mni(:,1:3) = vertices_mni(:,1:3) - mean(vertices_mni(:,1:3));
% Define the angle in degrees
theta = -90;

% Create the rotation matrix using cosd and sind (which take degrees)
Rz = [cosd(theta) -sind(theta) 0;
      sind(theta)  cosd(theta) 0;
      0            0           1];

% Assuming points is an Nx3 matrix of coordinates:
vertices_mni = (Rz * vertices_mni')';
vertices_mni = vertices_mni(:,1:3);
vertices_mni(:,2) = vertices_mni(:,2) - 20;
vertices_mni(:,3) = vertices_mni(:,3) + 10;



% smoothed_vertices = vertices_mni;

% Set the target reduction ratio (e.g., 0.5 means keeping 50% of the original faces)
targetReduction = 0.0500;

% Use reducepatch to decimate the mesh
[reducedFaces, reducedVertices] = reducepatch(fv.faces, vertices_mni, targetReduction);

% Create a new structure for the reduced surface
fv_reduced.faces = reducedFaces;
fv_reduced.vertices = reducedVertices;

% iterations = 5;   % number of smoothing iterations
% lambda = 0.1;     % smoothing strength (between 0 and 1)
% smoothed_vertices = laplacian_smooth(fv_reduced.vertices, fv_reduced.faces, iterations, lambda);

exportToVTK(fv_reduced.vertices, fv_reduced.faces ,'/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/MNIsurf.vtk')

figure;
p = patch('Faces', fv_reduced.faces, 'Vertices', fv_reduced.vertices);
p.FaceColor = [0.75, 0.75, 0.75]; hold on;
p.EdgeColor = 'none';
p.FaceAlpha = 0.2;
hems = {'L','R'};
regionCoordsAll = {}; c = 1;
regionAll = {}; 
for hmi = 1:2
    for iR = 1:length(regions)
        regionCoords = [];
        for iS = 1:length(subj_list_full)
            subject = subj_list_full{iS};
            filename = sprintf('/space/seh10/6/halgdev/projects/iverzh/data/Sternberg/preprocess/NativeSamp/%s_LFP_micro.mat', subject);
            load(filename, 'chan_coords', 'chan_locations')
            chan_coordsScale = chan_coords;


            chan_labels = regexprep(chan_locations, '[^a-zA-Z_]', '');



            splitPattern = '(_right|_left)';
            locations = regexp(chan_labels, splitPattern, 'split')';
            locations(cellfun(@(X) isempty(X), locations)) = [];
            locations = strrep(locations, 'ventral_medial_prefrontal_cortex', 'OFC');
            locations = strrep(locations, 'dorsal_anterior_cingulate_cortex', 'ACC');
            locations = strrep(locations, 'pre_supplementary_motor_area', 'SMA');
            locations = strrep(locations, 'amygdala', 'AMY');
            locations = strrep(locations, 'hippocampus', 'HIP');
            
            hem = regexp(chan_labels, splitPattern, 'match')';
            hemi = contains(hem, 'left');
            locations(hemi) = cellfun(@(X) ['L' X], locations(hemi),  'UniformOutput', false);
            hemi = contains(hem, 'right');
            locations(hemi) = cellfun(@(X) ['R' X], locations(hemi),  'UniformOutput', false);




             iiR = find(contains(locations, sprintf('%s%s', hems{hmi}, regions{iR})));
             if isempty(iiR); continue; end
             regionCoords = [regionCoords; chan_coordsScale(iiR,:)];


        %     for ii = 1:size(chan_coordsScale, 1)
        %         p3  = plot3(chan_coordsScale(ii,1), chan_coordsScale(ii,2), chan_coordsScale(ii,3), '.'); hold on;
        %         p3.MarkerSize = 18;
        %         iiR = find(contains(regions, locations{ii}));
        %         p3.Color = regionColors(iiR, :);
        %     end
        end

        filename = sprintf('/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/elec_%s%s.vtk', hems{hmi}, regions{iR});
        exportToVTK(regionCoords, [] ,filename, 'PointsOnly', true)
        
        regionCoordsAll{c} = regionCoords;
        regionsAll{c} = sprintf('%s%s', hems{hmi}, regions{iR});
        c = c + 1;
        
    end
end
% lighting for 3D visualization
lighting gouraud
camlight headlight
material dull

axis equal
xlabel('X (MNI)');
ylabel('Y (MNI)');
zlabel('Z (MNI)');
title('Transparent Cortical Surface in MNI Space');
view([90, 0]);
 axis vis3d;
rotate3d on;
fig = gcf;
fig.Color = 'w';

%%

addpath('/usr/pubsw/packages/freesurfer/RH4-x86_64-R600/matlab')%read_curv for sulc
addpath('/usr/pubsw/packages/MMPS/MMPS_235/matlab/fsurf_tools')
addpath('/home/pubsw/packages/MMPS/MMPS_235/matlab/mmil_utils')


close all



viewAng = [90 0;...
    90 270;...
    270 0];
h={'lh','rh'};

parcType = 'fusSplit7';
ih = 1;
iView=1;

[srf.(h{ih}).vertices, srf.(h{ih}).faces] = read_surf(sprintf('/space/seh8/5/halgdev/projects/jgarret/FW_data/electrode_location_files/fsaverage/surf/%s.pial','lh'));
srf.(h{ih}).faces=srf.(h{ih}).faces+1;
srf.(h{ih}).vertNorm = patchnormals(srf.(h{ih}));
figTmp=figure;
set(figTmp,'Position',[1 1 1080 1080])
set(figTmp,'InvertHardCopy','off')
set(figTmp,'Units','inches')
set(figTmp,'PaperOrientation','landscape')
subtightplot(1,1,1,[0 0],[0 0],[0,0])

load(sprintf('/space/seh8/5/halgdev/projects/jgarret/FW_data/brainView/FWparcAnnot_%s_%s.mat',parcType,'lh'),'val')
load(sprintf('/space/seh8/5/halgdev/projects/jgarret/FW_data/brainView/parcClrs_%s.mat',parcType))

pH = trisurf(srf.(h{ih}).faces,srf.(h{ih}).vertices(:,1),srf.(h{ih}).vertices(:,2),srf.(h{ih}).vertices(:,3),...
    'LineStyle','none','FaceVertexCData',val/255,'CDataMapping','direct');
axis image off
set(figTmp,'color',[0 0 0])
camva(8)
pH.FaceAlpha = 1;
hold on
material dull;lighting gouraud
view(viewAng(iView,1)+180*ih,viewAng(iView,2))
delete(findobj(gcf,'type','light'));
camlight headlight;
drawnow
ax=gca;

ih = 1;
exportToVTK(srf.(h{ih}).vertices, srf.(h{ih}).faces ,sprintf('/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/fssurf_%s.vtk', h{ih}))
ih = 2;
exportToVTK(srf.(h{ih}).vertices, srf.(h{ih}).faces ,sprintf('/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/fssurf_%s.vtk', h{ih}))




%%

function smoothed_pts = laplacian_smooth(vertices, faces, iterations, lambda)
    % vertices: Nx3 matrix of points
    % faces: face connectivity (from isosurface)
    % iterations: number of smoothing iterations
    % lambda: smoothing factor (0 to 1)

    smoothed_pts = vertices;

    % Build adjacency matrix
    TR = triangulation(faces, vertices);
    A = vertexAttachments(TR);
    
    for iter = 1:iterations
        new_pts = smoothed_pts;
        for v = 1:size(vertices,1)
            neighbors_idx = unique(TR.ConnectivityList(cell2mat(A(v)),:));
            neighbors_idx(neighbors_idx == v) = [];
            neighbors = smoothed_pts(neighbors_idx,:);
            avg_neighbor = mean(neighbors,1);
            
            % Update point towards neighbors
            new_pts(v,:) = smoothed_pts(v,:) + lambda*(avg_neighbor - smoothed_pts(v,:));
        end
        smoothed_pts = new_pts;
    end
end
