% allen_browser_test
cd('C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF')

tv = readNPY('template_volume_10um.npy'); % grey-scale "background signal intensity"
av = readNPY('annotation_volume_10um_by_index.npy'); % the number at each pixel labels the area, see note below
st = loadStructureTree('structure_tree_safe.csv'); % a table of what all the labels mean

% slice viewer
f = allenAtlasBrowser(tv, av, st);

% wire mesh
bregma = allenCCFbregma();
isBrain = av > 1; % >0 for original av, >1 for by_index
figure;
gridIn3D(double(isBrain),0.5,50,bregma);
axis vis3d
set(gca, 'ZDir', 'reverse')
axis equal
axis off
view([-30,25]);

% Outline the CP
target_structure = 'CP';
st_idx = find(strcmp(st.acronym,target_structure));
target_structure_color = hex2dec(reshape(st.color_hex_triplet{st_idx},3,[]))./255;
slice_spacing = 10;
structure_patch = isosurface(permute(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) == st_idx,[3,1,2]),0);
structure_wire = reducepatch(structure_outline.faces,structure_outline.vertices,0.05);

figure;
hold on;
structure_patch_h = patch('Vertices',structure_patch.vertices*slice_spacing, ...
    'Faces',structure_patch.faces, ...
    'FaceColor',target_structure_color,'EdgeColor','none');
structure_wire_h = patch('Vertices',structure_wire.vertices*slice_spacing, ...
    'Faces',structure_wire.faces, ...
    'FaceColor','none','EdgeColor',target_structure_color);
camlight
lighting gouraud

% Plot an arbitrary slice

bregma = allenCCFbregma();
isBrain = av > 1; % >0 for original av, >1 for by_index
figure;
gridIn3D(double(isBrain), 0.5, 50, bregma);
axis vis3d
set(gca, 'ZDir', 'reverse')
axis equal
axis off
view([-30,25]);

max_slice_size = ceil(sqrt(size(tv,1)^2+size(tv,2)^2));

hsp = surf(1:max_slice_size,1:max_slice_size,zeros(max_slice_size)+500,'EdgeColor','none');
rotate(hsp,[1,-1,1],30)
xd = round(hsp.XData);
yd = round(hsp.YData);
zd = round(hsp.ZData);

use_xd = xd > 0 & xd < size(tv,1);
use_yd = yd > 0 & yd < size(tv,3);
use_zd = zd > 0 & zd < size(tv,2);
use_slice_px = use_xd & use_yd & use_zd;
curr_slice_idx = sub2ind(size(tv),xd(use_slice_px),zd(use_slice_px),yd(use_slice_px));

curr_slice_isbrain = false(max_slice_size);
curr_slice_isbrain(use_slice_px) = isBrain(curr_slice_idx);

grab_pix_idx = sub2ind(size(tv),xd(curr_slice_isbrain),zd(curr_slice_isbrain),yd(curr_slice_isbrain));

curr_slice = nan(max_slice_size);
curr_slice(curr_slice_isbrain) = tv(grab_pix_idx);

set(hsp,'CData',curr_slice);





% TO DO: toggle rotate about probe


%% Different strategy from even grids: downsample polygon

figure; hold on;
% Plot outside of the brain

slice_spacing = 10;

target_volume = permute(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) > 1,[3,1,2]);

structure_patch = isosurface(target_volume,0);
structure_wire = reducepatch(structure_patch.faces,structure_patch.vertices,0.01);

target_structure_color = [0.7,0.7,0.7];

% structure_patch_h = patch('Vertices',structure_patch.vertices*slice_spacing, ...
%     'Faces',structure_patch.faces, ...
%     'FaceColor',target_structure_color,'EdgeColor','none');

structure_wire_h = patch('Vertices',structure_wire.vertices*slice_spacing, ...
    'Faces',structure_wire.faces, ...
    'FaceColor','none','EdgeColor',target_structure_color);


axis vis3d
set(gca, 'ZDir', 'reverse')
axis equal
axis off
view([-30,25]);

% material(structure_patch_h,'dull');
% camlight
% lighting gouraud

% Outline the CP
target_structure = 'CP';
st_idx = find(strcmp(st.acronym,target_structure));
target_structure_color = hex2dec(reshape(st.color_hex_triplet{st_idx},3,[]))./255;
slice_spacing = 10;
structure_patch = isosurface(permute(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) == st_idx,[3,1,2]),0);
structure_wire = reducepatch(structure_patch.faces,structure_patch.vertices,0.01);

% structure_patch_h = patch('Vertices',structure_patch.vertices*slice_spacing, ...
%     'Faces',structure_patch.faces, ...
%     'FaceColor',target_structure_color,'EdgeColor','none');
structure_wire_h = patch('Vertices',structure_wire.vertices*slice_spacing, ...
    'Faces',structure_wire.faces, ...
    'FaceColor','none','EdgeColor',target_structure_color);


% Set up the probe

[ml_max,dv_max,ap_max] = size(tv);


probe_ref_top = [ap_max,bregma(3),0];
probe_ref_bottom = [0,bregma(3),0];
probe_vector = [probe_ref_top',probe_ref_bottom'];
probe_ref_h = line(probe_vector(1,:),probe_vector(2,:),probe_vector(3,:),'linewidth',4,'color','k','linestyle','--');

probe_angle_ml = 10;
probe_angle_dv = 45;

rotate(probe_ref_h,[0,0,1],probe_angle_ml,probe_ref_top)
rotate(probe_ref_h,[0,1,0],probe_angle_dv,probe_ref_top)


% set up probe-angle plane

% coordinates in plot: [ap,ml,dv]

max_slice_size = ceil(sqrt(size(tv,1)^2+size(tv,2)^2));

hsp = surf(1:max_slice_size,1:max_slice_size,zeros(max_slice_size)+probe_ref_top(end),'EdgeColor','none');
rotate(hsp,[0,0,1],probe_angle_ml,[ap_max,bregma(3),0])
rotate(hsp,[0,1,0],probe_angle_dv,[ap_max,bregma(3),0])

xd = round(hsp.XData);
yd = round(hsp.YData);
zd = round(hsp.ZData);

use_xd = xd > 0 & xd < size(tv,1);
use_yd = yd > 0 & yd < size(tv,3);
use_zd = zd > 0 & zd < size(tv,2);
use_slice_px = use_xd & use_yd & use_zd;
curr_slice_idx = sub2ind(size(tv),xd(use_slice_px),zd(use_slice_px),yd(use_slice_px));

curr_slice_isbrain = false(max_slice_size);
curr_slice_isbrain(use_slice_px) = isBrain(curr_slice_idx);

grab_pix_idx = sub2ind(size(tv),xd(curr_slice_isbrain),zd(curr_slice_isbrain),yd(curr_slice_isbrain));

curr_slice = nan(max_slice_size);
curr_slice(curr_slice_isbrain) = tv(grab_pix_idx);

set(hsp,'CData',curr_slice);


% move the probe to a point
ap_offset = 20;
ml_offset = 20;

set(probe_ref_h,'XData',get(probe_ref_h,'XData') + ap_offset);
set(hsp,'XData',get(hsp,'XData') + ap_offset)

set(probe_ref_h,'YData',get(probe_ref_h,'YData') + ml_offset);
set(hsp,'YData',get(hsp,'YData') + ml_offset)





