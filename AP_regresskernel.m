function [k,predicted_signals,explained_var] = AP_regresskernel(regressors,signals,t_shifts,lambdas,zs,cvfold)
% [k,predicted_signals,explained_var] = AP_regresskernel(regressors,signals,t_shifts,lambdas,zs,cvfold)
%
% Linear regression of kernel from regressors to outputs
% 
% Inputs:
% regressors - dim x time
% signals - dim x time, one kernel will be returned for each dim
% t_shifts - time shifts of regressors
% lambda - ridge regression regularization (can be one for each regressor)
% zs - zscore the regressors to return beta weights
% cvfold - fold cross-validation 
%
% FOR MULTIPLE MODALITIES: make regressors and t_shifts cell arrays
%
% Outputs: 
% k - kernel
% predicted_signals - regressors*k
% explained_var - .total (total model), 
%                 .reduced (unique - regressors x signals)
%
% NOTE IF GPU ERROR: 
% The GPU is by default set to time out quickly if it can't update the
% monitor because it's busy. If there is a GPU error, this time (TDR) can be
% extended by adding a registry key: regedit, HKEY_LOCAL_MACHINE > SYSTEM >
% CurrentControlSet > Control > GraphicsDrivers, add new REG_DWORD called
% "TdrDelay" and having a value of something like 30 (seconds)

% Convert regressors and t_shifts to cells if not already
if ~iscell(regressors)
    regressors = {regressors};
end
if ~iscell(t_shifts)
    t_shifts = {t_shifts};
end

% Z-score all regressors and signals to get beta weights if selected
if ~exist('zs','var') || isempty(zs)
    zs = false;
end
if zs
    regressors = cellfun(@(x) zscore(x,[],2),regressors,'uni',false);
    signals = zscore(signals,[],2);
end

% Set cross-validation to 1-fold if not entered
if ~exist('cvfold','var') || isempty(cvfold)
   cvfold = 1; 
end

% Create design matrix of all time-shifted regressors
regressor_design = cellfun(@(regressors,t_shifts) repmat(regressors', ...
    [1,1,length(t_shifts)]),regressors,t_shifts,'uni',false);

% Temporally shift each page
for curr_regressors = 1:length(regressor_design)
    for curr_kernel_frame = 1:length(t_shifts{curr_regressors});
        regressor_design{curr_regressors}(:,:,curr_kernel_frame) = ...
            circshift(regressor_design{curr_regressors}(:,:,curr_kernel_frame), ...
            [t_shifts{curr_regressors}(curr_kernel_frame),0,0]);
    end
end

regressor_design = cell2mat(cellfun(@(regressor_design) ...
    reshape(regressor_design,[],size(regressor_design,2)*size(regressor_design,3)), ...
    regressor_design,'uni',false));

% Ridge regression for reducing noise: add offsets to design matrix to penalize k
if exist('lambdas','var') && any(lambdas)
    if length(lambdas) == 1
        ridge_matrix = lambdas*eye(size(regressor_design,2));
    elseif length(lambdas) == length(regressors)
        lambda_vector = cell2mat(reshape(cellfun(@(reg,t,lam) repmat(lam,size(reg,1)*length(t),1), ...
            regressors,t_shifts,num2cell(lambdas),'uni',false),[],1));
        ridge_matrix = bsxfun(@times,eye(size(regressor_design,2)),lambda_vector);
    else
        error('Number of lambdas doesn''t match regressor groups');
    end
else
    ridge_matrix = [];
end

% Send everything to the GPU
regressors_gpu = gpuArray([regressor_design;ridge_matrix]);
signals_gpu = gpuArray([signals';zeros(length(ridge_matrix),size(signals,1))]);

% Regression (and cross validation if selected)
cv_partition = round(linspace(1,cvfold,size(regressor_design,1)))';
k_cv = zeros(size(regressors_gpu,2),size(signals,1),cvfold,'single');
predicted_signals = zeros(size(signals));
if length(regressors) > 1
    predicted_signals_reduced = zeros(length(regressors),size(signals,1),size(signals,2));
end

for curr_cv = 1:cvfold
    
    if cvfold == 1
        train_idx = true(size(regressors_gpu,1),1);
        test_idx = true(size(regressors_gpu,1),1);
    elseif cvfold > 1
        train_idx = [cv_partition ~= curr_cv;true(size(ridge_matrix,1),1)];
        test_idx = cv_partition == curr_cv;
    end
    
    % Used to do manual inv(fluor_gpu'*fluor_gpu)*fluor_gpu'*spikes_gpu
    % but looks like \ works fine on GPU
    k_cv(:,:,curr_cv) = ...
        gather(regressors_gpu(train_idx,:)\signals_gpu(train_idx,:));
    
    predicted_signals(:,test_idx(1:size(signals,2))) = ...
        gather(regressors_gpu(test_idx(1:size(signals,2)),:)*gpuArray(k_cv(:,:,curr_cv)))';
    
    % Predict on test set for each regressor modality (if > 1)
    if length(regressors) > 1
        regressor_split_size = cellfun(@(x) size(x,1),regressors).*cellfun(@length,t_shifts);
        regressor_split_idx = cellfun(@(x) setdiff(1:size(regressors_gpu,2),x), ...
            mat2cell(1:size(regressors_gpu,2),1,regressor_split_size),'uni',false);
        for curr_regressor = 1:length(regressors)
            curr_predicted = ...
                gather(regressors_gpu(test_idx,regressor_split_idx{curr_regressor})* ...
                gpuArray(k_cv(regressor_split_idx{curr_regressor},:,curr_cv)));
            
            predicted_signals_reduced(curr_regressor,:,:) = ...
                permute(curr_predicted(1:size(signals,2),:),[3,2,1]);
        end
    end  
end

% Total explained variance
sse_signals = sum(signals.^2,2);
sse_total_residual = sum(bsxfun(@minus,predicted_signals,signals).^2,2);
explained_var.total = (sse_signals - sse_total_residual)./sse_signals;

% Partial explained variance
if length(regressors) > 1
    sse_partial_residual = sum(bsxfun(@minus,predicted_signals_reduced,permute(signals,[3,1,2])).^2,3);
    explained_var.reduced = bsxfun(@rdivide,bsxfun(@minus,(sse_signals - sse_total_residual), ...
        bsxfun(@minus,sse_signals,sse_partial_residual)),(sse_signals - sse_total_residual));
end

% Get the final k from averaging
k = mean(k_cv,3);





%%%%%%
%%%%%%
%%%%%% 
% OLD w/o cross validation

% 
% % Spell out what to do - looks the same result and time as
% % fluor_gpu\spikes_gpu
% k = gather(inv(fluor_gpu'*fluor_gpu)*fluor_gpu'*spikes_gpu);
% 
% % Get predicted signals
% predicted_signals = gather(fluor_gpu*gpuArray(k));
% predicted_signals = predicted_signals(1:size(signals,2),:);
% 
% sse_signals = sum(signals.^2,2)';
% sse_total_residual = sum(bsxfun(@minus,predicted_signals,signals').^2,1);
% explained_var_total = (sse_signals - sse_total_residual)./sse_signals;
% 
% % Get explained variance for each regressor group (reduced: leave-one-out)
% regressor_split_size = cellfun(@(x) size(x,1),regressors).*cellfun(@length,t_shifts);
% regressor_split_idx = cellfun(@(x) setdiff(1:size(fluor_gpu,2),x), ...
%     mat2cell(1:size(fluor_gpu,2),1,regressor_split_size),'uni',false);
% predicted_signals_reduced = zeros(length(regressors),size(signals,1),size(signals,2));
% for curr_regressor = 1:length(regressors)
%     curr_predicted = ...
%         gather(fluor_gpu(:,regressor_split_idx{curr_regressor})* ...
%         gpuArray(k(regressor_split_idx{curr_regressor},:)));
%     
%     predicted_signals_reduced(curr_regressor,:,:) = ...
%         permute(curr_predicted(1:size(signals,2),:),[3,2,1]);
% end
% 
% sse_partial_residual = sum(bsxfun(@minus,predicted_signals_reduced,permute(signals,[3,1,2])).^2,3);
% explained_var_reduced = bsxfun(@rdivide,bsxfun(@minus,(sse_signals - sse_total_residual), ...
%     bsxfun(@minus,sse_signals,sse_partial_residual)),(sse_signals - sse_total_residual));
% 
% % Package explained variance
% explained_var.total = explained_var_total;
% explained_var.reduced = explained_var_reduced;


