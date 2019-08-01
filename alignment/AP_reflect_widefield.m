function im_reflected = AP_reflect_widefield(im)
% Reflect image about the midline
% (currently only uses master alignment)
% this doesn't work perfectly so something's off - I've spent way too much
% time on it and it's too close to be worth not just running with at the
% moment

% Load bregma and master CCF tform
bregma = allenCCFbregma;
bregma(3) = bregma(3) + 0.5;
alignment_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment';
ccf_tform_fn = [alignment_path filesep 'ccf_tform.mat'];
load(ccf_tform_fn);

midline_align = [1,0,0]*ccf_tform.T;
midline_vector = [midline_align([1,2])./norm(midline_align([1,2])),0];

um2pixel = 20.6;
bregma_resize = bregma*(10/um2pixel);
bregma_align = [bregma_resize([3,1]),1]*ccf_tform.T;

% Translate bregma to origin, then reflect, then translate back
I = eye(3);

reflection_tform = I - 2*(midline_vector'*midline_vector);

origin_translation_tform = I;
origin_translation_tform(3,1:2) = -bregma_align(1:2);

bregma_translation_tform = I;
bregma_translation_tform(3,1:2) = bregma_align(1:2);

tform_matrix = origin_translation_tform*reflection_tform*bregma_translation_tform;

tform = affine2d;
tform.T = tform_matrix;
im_reflected = imwarp(im,tform,'Outputview',imref2d(size(im)));








