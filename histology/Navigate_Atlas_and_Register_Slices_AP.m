% ------------------------------------------------------------------------
%          Run Allen Atlas Browser
% ------------------------------------------------------------------------


%% ENTER FILE LOCATION AND PROBE-SAVE-NAME


% directory of histology
processed_images_folder = 'C:\Users\Andrew\Desktop\test_histology\processed\ps_slices'; 

% name the saved probe points, to avoid overwriting another set of probes going in the same folder
probe_save_name_suffix = 'striatum'; 

% directory of reference atlas files
annotation_volume_location = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF\annotation_volume_10um_by_index.npy';
structure_tree_location = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF\structure_tree_safe_2017.csv';
template_volume_location = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF\template_volume_10um.npy';




%% GET PROBE TRAJECTORY POINTS

% load the reference brain and region annotations
if ~exist('av','var') || ~exist('st','var') || ~exist('tv','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
    tv = readNPY(template_volume_location);
end

% create Atlas viewer figure
f = figure('Name','Atlas Viewer'); 

%%% AP: These two GUIs should be nested to work properly

% show histology in Slice Viewer
try; figure(slice_figure_browser); title('');
catch; slice_figure_browser = figure('Name','Slice Viewer'); end
reference_size = size(tv);
sliceBrowser(slice_figure_browser, processed_images_folder, f, reference_size);


% use application in Atlas Transform Viewer
% use this function if you have a processed_images_folder with appropriately processed .tif histology images
f = AtlasTransformBrowser(f, tv,av,st, slice_figure_browser, processed_images_folder, probe_save_name_suffix); 


% use the simpler version, which does not interface with processed slice images
% just run these two lines instead of the previous 5 lines of code
% 
% save_location = processed_images_folder;
% f = allenAtlasBrowser(tv,av,st, save_location, probe_save_name_suffix);



%% AP version 




% directory of histology
processed_images_folder = 'C:\Users\Andrew\Desktop\test_histology\processed\ps_slices'; 

% name the saved probe points, to avoid overwriting another set of probes going in the same folder
probe_save_name_suffix = 'striatum'; 

% directory of reference atlas files
annotation_volume_location = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF\annotation_volume_10um_by_index.npy';
structure_tree_location = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF\structure_tree_safe_2017.csv';
template_volume_location = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF\template_volume_10um.npy';




%% GET PROBE TRAJECTORY POINTS

if ~exist('av','var') || ~exist('st','var') || ~exist('tv','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
    tv = readNPY(template_volume_location);
end


AtlasTransformBrowser_AP(tv,av,st, processed_images_folder, probe_save_name_suffix); 

% new function
AP_align_histology_ccf(tv,av,st, processed_images_folder); 











