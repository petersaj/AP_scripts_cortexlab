function [k,predicted_signals,explained_var,predicted_signals_reduced] = ...
    AP_regresskernel(regressors,signals,t_shifts,lambdas,zs,cvfold,return_constant)
% [k,predicted_signals,explained_var,predicted_signals_reduced] = AP_regresskernel(regressors,signals,t_shifts,lambdas,zs,cvfold)
%
% Linear regression of kernel from regressors to outputs
% (note constant term is included and used for prediction, but not output) 
%
% Inputs:
% regressors - dim x time
% signals - dim x time, one kernel will be returned for each dim
% t_shifts - time shifts of regressors
% lambda - ridge regression regularization (can be one for each regressor)
% zs - length 2 logical vector: zscore 1) regressors, 2) signals
% (default is [false,true])
% cvfold - fold cross-validation 
% return_constant - if true, include constant in kernel (false by default)
%
% FOR MULTIPLE MODALITIES: make regressors and t_shifts cell arrays
%
% Outputs: 
% k - kernel
% predicted_signals - regressors*k
% explained_var - .total (total model), 
%                 .reduced (unique - regressors x signals)
% predicted_signals_reduced - predicted signals with reduced k
%
% NOTE IF GPU ERROR: 
% The GPU is by default set to time out quickly if it can't update the
% monitor because it's busy. If there is a GPU error, this time (TDR) can be
% extended by adding a registry key: regedit, HKEY_LOCAL_MACHINE > SYSTEM >
% CurrentControlSet > Control > GraphicsDrivers, add new REG_DWORD called
% "TdrDelay" and having a value of something like 30 (seconds)

% Use no temporal delay if none specified
if ~exist('t_shifts','var') || isempty(t_shifts)
    t_shifts = 0;
end

% Convert regressors and t_shifts to cells if not already
if ~iscell(regressors)
    regressors = {regressors};
end
if ~iscell(t_shifts)
    t_shifts = {t_shifts};
end

% Standardize orientations
regressors = reshape(regressors,1,[]);
t_shifts = reshape(t_shifts,1,[]);

% Z-score all regressors and signals to get beta weights if selected
if ~exist('zs','var') || isempty(zs)
    zs = [false,true];
end
if zs(1)
    regressors = cellfun(@(x) zscore(x,[],2),regressors,'uni',false);
end
if zs(2)
    signals = zscore(signals,[],2);
end

% Set cross-validation to 1-fold if not entered
if ~exist('cvfold','var') || isempty(cvfold)
   cvfold = 1; 
end

% Set return_constant false if not entered
if ~exist('return_constant','var') || isempty(cvfold)
   return_constant = false; 
end

% Create design matrix of all time-shifted regressors
regressor_design = cellfun(@(regressors,t_shifts) repmat(regressors', ...
    [1,1,length(t_shifts)]),regressors,t_shifts,'uni',false);

% Temporally shift each page
for curr_regressors = 1:length(regressor_design)
    for curr_kernel_frame = 1:length(t_shifts{curr_regressors})
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
        ridge_matrix = lambdas*eye(size(regressor_design,2)+1);
    elseif length(lambdas) == length(regressors)
        lambda_vector = cell2mat(reshape(cellfun(@(reg,t,lam) repmat(lam,size(reg,1)*length(t),1), ...
            regressors,t_shifts,num2cell(lambdas),'uni',false),[],1));
        ridge_matrix = bsxfun(@times,eye(size(regressor_design,2)+1),lambda_vector);
    else
        error('Number of lambdas doesn''t match regressor groups');
    end
else
    ridge_matrix = [];
end

% Prepare column of 1's to have a constant term
constant = ones(size(regressor_design,1),1);

% Send everything to the GPU
regressors_gpu = gpuArray([[regressor_design,constant];ridge_matrix]);
signals_gpu = gpuArray([signals';zeros(length(ridge_matrix),size(signals,1))]);

% Regression (and cross validation if selected)
k_cv = zeros(size(regressors_gpu,2),size(signals,1),cvfold,'single');
% (use randomly distributed time points for cross validation)
cv_partition_ordered = round(linspace(1,cvfold,size(regressor_design,1)))';
cv_partition = cv_partition_ordered(randperm(length(cv_partition_ordered)));

predicted_signals = nan(size(signals));
predicted_signals_reduced = nan(size(signals,1),size(signals,2),length(regressors));
for curr_cv = 1:cvfold
    
    % Get training/test sets
    if cvfold == 1
        train_idx = true(size(regressors_gpu,1),1);
        test_idx = true(size(regressor_design,1),1);
    else
        train_idx = [cv_partition ~= curr_cv;true(size(ridge_matrix,1),1)];
        test_idx = cv_partition == curr_cv;
    end
    
    % If regressors for training fold are empty, warning
    if ~all(any(regressors_gpu(train_idx,:),1))
        warning('Regressors in fold unfilled (not enough trials?)');
    end
    
    % Used to do manual inv(fluor_gpu'*fluor_gpu)*fluor_gpu'*spikes_gpu
    % but looks like \ works fine on GPU
    k_cv(:,:,curr_cv) = ...
        gather(regressors_gpu(train_idx,:)\signals_gpu(train_idx,:));
    
    predicted_signals(:,test_idx) = ...
        gather(regressors_gpu(test_idx,:)*gpuArray(k_cv(:,:,curr_cv)))';
    
    % Predict on test set for each regressor modality (if > 1)
    if length(regressors) > 1
        
        regressor_split_size = [cellfun(@(x) size(x,1),regressors).*cellfun(@length,t_shifts),1];
        
%         (use everything BUT particular regressor)
        regressor_split_idx = cellfun(@(x) setdiff(1:size(regressors_gpu,2),x), ...
            mat2cell(1:size(regressors_gpu,2),1,regressor_split_size),'uni',false);

        % (use ONLY particular regressor)
%         regressor_split_idx = mat2cell(1:size(regressors_gpu,2),1,regressor_split_size);
        
        for curr_regressor = 1:length(regressors)
            
            curr_predicted = ...
                gather(regressors_gpu(test_idx,regressor_split_idx{curr_regressor})* ...
                gpuArray(k_cv(regressor_split_idx{curr_regressor},:,curr_cv)));
            
            predicted_signals_reduced(:,test_idx,curr_regressor) = ...
                curr_predicted';
            
        end
    end  
end

% Throw errors for bad numbers in kernel or predicted signals
if any(isinf(k_cv(:))) || any(isinf(predicted_signals(:))) || ...
        any(isnan(k_cv(:))) || any(isnan(predicted_signals(:)))
    error('Inf/NaN in kernel or predicted signal')
end

% Total explained variance
sse_signals = sum(signals.^2,2);
sse_total_residual = sum(bsxfun(@minus,predicted_signals,signals).^2,2);
explained_var.total = (sse_signals - sse_total_residual)./sse_signals;

% Partial explained variance
sse_complete_model = sum(predicted_signals.^2,2);
if length(regressors) > 1
    sse_partial = cell2mat(arrayfun(@(x) ...
        sum(predicted_signals_reduced(:,:,x).^2,2),1:length(regressors),'uni',false));
    sse_partial_residual = cell2mat(arrayfun(@(x) ...
        sum(bsxfun(@minus,predicted_signals_reduced(:,:,x),signals).^2,2),1:length(regressors),'uni',false));
    % (I think this gets unique + confounded expl var?)
%     explained_var.reduced = bsxfun(@rdivide,bsxfun(@minus,sse_signals,sse_partial_residual),sse_signals);
    % (this method + leave-one-out above should be unique expl var?)
    explained_var.reduced = bsxfun(@rdivide,bsxfun(@minus,sse_complete_model,sse_partial),sse_signals);
end

% Get the final k from averaging
% (to remove the constant term)
if ~return_constant
    k = mean(k_cv(1:end-1,:,:),3);
    % (to keep constant term)
elseif return_constant
    k = mean(k_cv,3);
end










