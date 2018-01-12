function [tform_matrix,im_aligned] = AP_align_widefield(animal,days)
% [tform_matrix,im_aligned] = AP_align_widefield(animal,days)
% 
% Get transform matrix to align widefield images across days
% (this is not comprehensively tested and isn't always good)

% what looks best at the moment: align by purple with large border and no
% blurring

avg_im_blue = cell(length(days),1);
for curr_day = 1:length(days)    
    [img_path,img_exists] = AP_cortexlab_filename(animal,days{curr_day},[],'imaging');    
    avg_im_blue{curr_day} = readNPY([img_path filesep 'meanImage_blue.npy']);
end

avg_im_purple = cell(length(days),1);
for curr_day = 1:length(days)
    [img_path,img_exists] = AP_cortexlab_filename(animal,days{curr_day},[],'imaging');
    avg_im_purple{curr_day} = readNPY([img_path filesep 'meanImage_purple.npy']);
end

% Choose reference image
ref_im_full = avg_im_purple{1};
ref_size = size(ref_im_full);

border_pixels = 50;

% Remove border for aligning
im_align = cellfun(@(x) x(border_pixels:end-border_pixels+1,border_pixels:end-border_pixels+1),avg_im_purple,'uni',false);
ref_im = ref_im_full(border_pixels:end-border_pixels+1,border_pixels:end-border_pixels+1);

fprintf('Registering average images:')
fprintf('\n');

tform_matrix = cell(length(avg_im_blue),1);

im_aligned = nan(ref_size(1),ref_size(2),length(days));

for curr_session = 1:length(days);
    
    % OLD
    
    %     [optimizer, metric] = imregconfig('monomodal');
    %     optimizer.MaximumStepLength = 0.02;
    %     %optimizer.MaximumStepLength = 0.0001;
    %     optimizer.RelaxationFactor = 0.1;
    %     optimizer.GradientMagnitudeTolerance = 1e-5;
    %     optimizer.MaximumIterations = 300;
    %     %         optimizer = registration.optimizer.OnePlusOneEvolutionary();
    %     %         optimizer.MaximumIterations = 500;
    %     %         optimizer.GrowthFactor = 1.00001;
    %     %         optimizer.InitialRadius = 1e-5;
    %
    %     % Perform the registration on the maximum image
    %     tformEstimate = imregtform(avg_im{curr_session},avg_im{1},'affine',optimizer,metric);
    %     %tformEstimate = imregcorr(summed_max(:,:,curr_session),summed_max(:,:,1),'similarity');
    %
    %     curr_im_reg = imwarp(avg_im{curr_session},tformEstimate,'Outputview',imref2d(size(avg_im{1})));
    %
    %     tform_matrix{curr_session} = tformEstimate.T;
    %     summed_max_reg(:,:,curr_session) = max_im_reg;
    %     summed_mean_reg(:,:,curr_session) = mean_im_reg;
     
    
    %%%%%%%%%%%%%
    
    % This is to do correlation, then affine (if above doesn't work)
    [optimizer, metric] = imregconfig('monomodal');
    optimizer = registration.optimizer.OnePlusOneEvolutionary();
    optimizer.MaximumIterations = 200;
    optimizer.GrowthFactor = 1+1e-6;
    optimizer.InitialRadius = 1e-4;
    
%     % Register 1) correlation
%     tformEstimate_corr = imregcorr(im_align{curr_session},im_align{ref_im_num},'similarity');
%     curr_im_reg_corr = imwarp(im_align{curr_session},tformEstimate_corr,'Outputview',imref2d(size(im_align{ref_im_num})));
%     
%     % Register 2) affine
%     tformEstimate_affine = imregtform(curr_im_reg_corr,im_align{ref_im_num},'affine',optimizer,metric);
%     
%     tformEstimate_combined = tformEstimate_corr;
%     tformEstimate_combined.T = tformEstimate_affine.T*tformEstimate_corr.T;
%     
%     curr_im_reg = imwarp(avg_im{curr_session},tformEstimate_combined,'Outputview',imref2d(size(avg_im{ref_im_num})));
%     
%     tform_matrix{curr_session} = tformEstimate_combined.T;

    %%% for just affine
    tformEstimate_affine = imregtform(im_align{curr_session},ref_im,'affine',optimizer,metric);
    curr_im_reg = imwarp(avg_im_blue{curr_session},tformEstimate_affine,'Outputview',imref2d(ref_size));
    tform_matrix{curr_session} = tformEstimate_affine.T;
    %%%%
    
    im_aligned(:,:,curr_session) = curr_im_reg;
    
    AP_print_progress_fraction(curr_session,length(days))
    
end
















% error('trying a second round of alignment: fix this first');
% 
% % Choose reference image
% ref_im_full = nanmean(im_aligned,3);
% ref_size = size(ref_im_full);
% 
% border_pixels = 50;
% 
% % Remove border for aligning
% im_aligned_cell = mat2cell(im_aligned,size(im_aligned,1),size(im_aligned,2),ones(size(im_aligned,3),1));
% im_align = cellfun(@(x) x(border_pixels:end-border_pixels+1,border_pixels:end-border_pixels+1),im_aligned_cell,'uni',false);
% ref_im = ref_im_full(border_pixels:end-border_pixels+1,border_pixels:end-border_pixels+1);
% 
% fprintf('Registering average images:')
% fprintf('\n');
% 
% tform_matrix = cell(length(avg_im_blue),1);
% 
% im_aligned = nan(ref_size(1),ref_size(2),length(days));
% 
% for curr_session = 1:length(days);
%     
%     % OLD
%     
%     %     [optimizer, metric] = imregconfig('monomodal');
%     %     optimizer.MaximumStepLength = 0.02;
%     %     %optimizer.MaximumStepLength = 0.0001;
%     %     optimizer.RelaxationFactor = 0.1;
%     %     optimizer.GradientMagnitudeTolerance = 1e-5;
%     %     optimizer.MaximumIterations = 300;
%     %     %         optimizer = registration.optimizer.OnePlusOneEvolutionary();
%     %     %         optimizer.MaximumIterations = 500;
%     %     %         optimizer.GrowthFactor = 1.00001;
%     %     %         optimizer.InitialRadius = 1e-5;
%     %
%     %     % Perform the registration on the maximum image
%     %     tformEstimate = imregtform(avg_im{curr_session},avg_im{1},'affine',optimizer,metric);
%     %     %tformEstimate = imregcorr(summed_max(:,:,curr_session),summed_max(:,:,1),'similarity');
%     %
%     %     curr_im_reg = imwarp(avg_im{curr_session},tformEstimate,'Outputview',imref2d(size(avg_im{1})));
%     %
%     %     tform_matrix{curr_session} = tformEstimate.T;
%     %     summed_max_reg(:,:,curr_session) = max_im_reg;
%     %     summed_mean_reg(:,:,curr_session) = mean_im_reg;
%      
%     
%     %%%%%%%%%%%%%
%     
%     % This is to do correlation, then affine (if above doesn't work)
%     [optimizer, metric] = imregconfig('monomodal');
%     optimizer = registration.optimizer.OnePlusOneEvolutionary();
%     optimizer.MaximumIterations = 200;
%     optimizer.GrowthFactor = 1+1e-6;
%     optimizer.InitialRadius = 1e-4;
%     
% %     % Register 1) correlation
% %     tformEstimate_corr = imregcorr(im_align{curr_session},im_align{ref_im_num},'similarity');
% %     curr_im_reg_corr = imwarp(im_align{curr_session},tformEstimate_corr,'Outputview',imref2d(size(im_align{ref_im_num})));
% %     
% %     % Register 2) affine
% %     tformEstimate_affine = imregtform(curr_im_reg_corr,im_align{ref_im_num},'affine',optimizer,metric);
% %     
% %     tformEstimate_combined = tformEstimate_corr;
% %     tformEstimate_combined.T = tformEstimate_affine.T*tformEstimate_corr.T;
% %     
% %     curr_im_reg = imwarp(avg_im{curr_session},tformEstimate_combined,'Outputview',imref2d(size(avg_im{ref_im_num})));
% %     
% %     tform_matrix{curr_session} = tformEstimate_combined.T;
% 
%     %%% for just affine
%     
%     new_tform = imregtform(im_align{curr_session},ref_im,'affine',optimizer,metric);
%     combined_tform_t = new_tform.T*tform_matrix{curr_session};
%     new_tform.T = combined_tform_t;
%     
%     curr_im_reg = imwarp(im_aligned_cell{curr_session},new_tform,'Outputview',imref2d(ref_size));
%     tform_matrix{curr_session} = new_tform.T;
%     %%%%
%     
%     im_aligned(:,:,curr_session) = curr_im_reg;
%     
%     AP_print_progress_fraction(curr_session,length(days))
%     
% end
% 
% 





