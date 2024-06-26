function [k,predicted_signals,explained_var,predicted_signals_reduced] = ...
    AP_regresskernel(regressors,signals,t_shifts,lambdas,zs,cvfold,return_constant,use_constant,discontinuities)
% [k,predicted_signals,explained_var,predicted_signals_reduced] = ...
%     AP_regresskernel(regressors,signals,t_shifts,lambdas,zs,cvfold,return_constant,use_constant,discontinuities)
%
% Linear regression of kernel from regressors to outputs
%
% Inputs:
% regressors - dim x time
% signals - dim x time, one kernel will be returned for each dim
% t_shifts - time shifts of regressors
% lambda - ridge regression regularization (can be one for each regressor)
% zs - length 2 logical vector: zscore 1) regressors, 2) signals
% (default is [false,true])
% cvfold - fold cross-validation 
% return_constant - if true, return constant in output (false by default)
% use_constant - if true, include constant in kernel (true by default)
% discontinuities - samples at the start of discontinuities
%
% FOR MULTIPLE MODALITIES: make regressors and t_shifts cell arrays
%
% Outputs: 
% k - kernel (regressor x t_shifts x signals, cell array if multiple
% t_shifts)
% predicted_signals - regressors*k
% explained_var - .total (total model), 
%                 .partial (2 columns: 1 = alone, 2 = omitted)
% predicted_signals_reduced - predicted signals with reduced k
%
% NaN handling: signals with all NaNs are ignored, common NaNs across
% signals are ignored, unique NaNs across signals gives an error (handling
% this would require separate matrix divisions)
%
% NOTE IF GPU ERROR: 
% The GPU is by default set to time out quickly if it can't update the
% monitor because it's busy. If there is a GPU error, this time (TDR) can be
% extended by adding a registry key: regedit, HKEY_LOCAL_MACHINE > SYSTEM >
% CurrentControlSet > Control > GraphicsDrivers, add new DWORD (32bit Value) called
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
lambdas = reshape(lambdas,1,[]);

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
if ~exist('return_constant','var') || isempty(return_constant)
   return_constant = false; 
end

% Set use_constant true if not entered
if ~exist('use_constant','var') || isempty(use_constant)
   use_constant = true; 
end

% Set/check discontinuities
if ~exist('discontinuities','var') || isempty(discontinuities)
   discontinuities = zeros(1,size(signals,2)); 
elseif length(discontinuities) ~= size(signals,2)
    error('Discontinuities vector doesn''t match signals length')        
end
% (the binning and end are always discontinuities)
discontinuities([1,end]) = 1;
% (force discontinuities to be a column matrix)
discontinuities = reshape(discontinuities,[],1);

% Create design matrix of all time-shifted regressors
regressor_design = cellfun(@(regressors,t_shifts) repmat(regressors', ...
    [1,1,length(t_shifts)]),regressors,t_shifts,'uni',false);

% Temporally shift each page (regressors and discontinuities)
for curr_regressors = 1:length(regressor_design)
    
    % Set up discontinuities matrix (none at t shift == 0)
    curr_discontinuities = repmat(discontinuities,1, ...
        size(regressor_design{curr_regressors},2), ...
        size(regressor_design{curr_regressors},3));
    curr_discontinuities(:,:,t_shifts{curr_regressors} == 0) = 0;
    
    for curr_kernel_frame = 1:length(t_shifts{curr_regressors})
        regressor_design{curr_regressors}(:,:,curr_kernel_frame) = ...
            circshift(regressor_design{curr_regressors}(:,:,curr_kernel_frame), ...
            [t_shifts{curr_regressors}(curr_kernel_frame),0,0]);
        
        curr_discontinuities(:,:,curr_kernel_frame) = ...
            circshift(curr_discontinuities(:,:,curr_kernel_frame), ...
            [t_shifts{curr_regressors}(curr_kernel_frame),0,0]);
    end
    
    % Cumulative sum the discontinuities in appropriate direction
    curr_discontinuities_cumulative = curr_discontinuities > 0;
    curr_discontinuities_cumulative(:,:,t_shifts{curr_regressors} < 0) = ...
        cumsum(curr_discontinuities(:,:,t_shifts{curr_regressors} < 0),3,'reverse') > 0;
    curr_discontinuities_cumulative(:,:,t_shifts{curr_regressors} > 0) = ...
        cumsum(curr_discontinuities(:,:,t_shifts{curr_regressors} > 0),3) > 0;
    
    % Zero regressors predicting discontinuous / invalid locations
    regressor_design{curr_regressors}(curr_discontinuities_cumulative) = 0;
    
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

% Prepare column of 1's to have a constant term (if selected)
if use_constant
    constant = ones(size(regressor_design,1),1);
    % if there's a ridge matrix, add another row and column of zeros
    if ~isempty(ridge_matrix)
        ridge_matrix(end+1,end+1) = 0;
    end
else
    constant = [];
end

% If regressors matrix is sparse, make sparse matrix class (way faster)
sparse_regressors = sum(regressor_design(:) == 0)/(numel(regressor_design)) > 0.5;
if sparse_regressors
    regressors_gpu = sparse(double([[regressor_design,constant];ridge_matrix]));
    signals_gpu = sparse(double([signals';zeros(length(ridge_matrix),size(signals,1))]));
else
    % Otherwise, do on the GPU
    regressors_gpu = gpuArray([[regressor_design,constant];ridge_matrix]);
    signals_gpu = gpuArray([signals';zeros(length(ridge_matrix),size(signals,1))]);
end

% Regression (and cross validation if selected)

% Check for NaNs: 
% Any signals which are all NaN won't be regressed
predictable_signals = ~all(isnan(signals),2);
% If there are NaNs not shared across the remaining signals, error out
if ~all(all(isnan(signals(predictable_signals,:)),1) | ...
        all(~isnan(signals(predictable_signals,:)),1))
    error('NaN values vary across signals')
end
    
% Predictable samples: non-zero/nan regressors, non-nan signals
predictable_samples = ...
    any(regressor_design,2) & ...
    ~any(isnan(regressor_design),2) & ...
    ~any(isnan(signals(predictable_signals,:)),1)';
cv_partition = nan(size(regressor_design,1),1);
cv_partition(predictable_samples) = min(floor(linspace(1,cvfold+1,sum(predictable_samples))),cvfold)';

% (to shuffle CV points instead of using consecutive chunks)
% (NOT VALID FOR TIME SERIES - but necessary if not)
% cv_partition(predictable_samples) = AP_shake(round(linspace(1,cvfold,sum(predictable_samples)))');

k_cv = nan(size(regressors_gpu,2),size(signals,1),cvfold,'double');
predicted_signals = nan(size(signals));
predicted_signals_reduced = nan(size(signals,1),size(signals,2),length(regressors));
predicted_signals_partial = nan(size(signals,1),size(signals,2),length(regressors),2);
for curr_cv = 1:cvfold
    
    % Get training/test sets
    if cvfold == 1
        train_idx = [predictable_samples;true(size(ridge_matrix,1),1)];
        test_idx = predictable_samples;
    else
        train_idx = [cv_partition ~= curr_cv & predictable_samples;true(size(ridge_matrix,1),1)];
        test_idx = cv_partition == curr_cv;
    end
    
    % If regressors for training fold are empty, warning
    if ~all(any(regressors_gpu(train_idx,:),1))
        warning('Regressors in fold unfilled (not enough trials?)');
    end
    
    % Used to do manual inv(regressors_gpu'*regressors_gpu)*regressors_gpu'*spikes_gpu
    % but looks like \ works fine on GPU
    k_cv(:,predictable_signals,curr_cv) = ...
        full(gather(regressors_gpu(train_idx,:)\signals_gpu(train_idx,predictable_signals)));
    
    if sparse_regressors
        predicted_signals(predictable_signals,test_idx) = ...
            full(regressors_gpu(test_idx,:)*sparse(double(k_cv(:,predictable_signals,curr_cv))))';
    else
        predicted_signals(predictable_signals,test_idx) = ...
            gather(regressors_gpu(test_idx,:)*gpuArray(k_cv(:,predictable_signals,curr_cv)))';
    end
    
    % Reduced predictions on test set for each regressor modality (if > 1)
    if length(regressors) > 1
        
        if use_constant
            regressor_split_size = [cellfun(@(x) size(x,1),regressors).*cellfun(@length,t_shifts),1];
        else
            regressor_split_size = [cellfun(@(x) size(x,1),regressors).*cellfun(@length,t_shifts)];
        end
        
        % Get indicies to use regressor alone or omitted
        regressor_alone_idx = mat2cell(1:size(regressors_gpu,2),1,regressor_split_size);
        regressor_omitted_idx = cellfun(@(x) setdiff(1:size(regressors_gpu,2),x), ...
            mat2cell(1:size(regressors_gpu,2),1,regressor_split_size),'uni',false);
        
        % Predict the signal using kernel with current regressor omitted
        % (for use in isolating responses, e.g. stim response = full
        % response - predicted move response)
        for curr_regressor = 1:length(regressors)
            
            if sparse_regressors
                curr_predicted = ...
                    full(regressors_gpu(test_idx,regressor_omitted_idx{curr_regressor})* ...
                    sparse(double(k_cv(regressor_omitted_idx{curr_regressor},predictable_signals,curr_cv))));
            else
                curr_predicted = ...
                    gather(regressors_gpu(test_idx,regressor_omitted_idx{curr_regressor})* ...
                    gpuArray(k_cv(regressor_omitted_idx{curr_regressor},predictable_signals,curr_cv)));
            end
            
            predicted_signals_reduced(predictable_signals,test_idx,curr_regressor) = ...
                curr_predicted';
            
        end
               
        % Refit the kernel and predict for each regressor
        % (for calculating partial explained variance, both alone and
        % omitted)
        for curr_regressor = 1:length(regressors)
            for curr_partial = 1:2
                
                switch curr_partial
                    case 1
                        regressor_split_idx = regressor_alone_idx;
                    case 2
                        regressor_split_idx = regressor_omitted_idx;
                end
               
            k_cv_reduced = ...
                full(gather(regressors_gpu(train_idx,regressor_split_idx{curr_regressor})\ ...
                signals_gpu(train_idx,predictable_signals)));
                       
            if sparse_regressors
                curr_predicted = ...
                    full(regressors_gpu(test_idx,regressor_split_idx{curr_regressor})* ...
                    sparse(double(k_cv_reduced)));
            else
                curr_predicted = ...
                    gather(regressors_gpu(test_idx,regressor_split_idx{curr_regressor})* ...
                    gpuArray(k_cv_reduced));
            end
            
            predicted_signals_partial(predictable_signals,test_idx,curr_regressor,curr_partial) = ...
                curr_predicted';    
            
            end
        end
    end  
end

% Throw errors for bad numbers in kernel or predicted signals
if ...
        any(reshape(isnan(k_cv(:,predictable_signals,:)),[],1)) || ...
        any(reshape(isnan(predicted_signals(predictable_signals,predictable_samples)),[],1)) || ...
        any(reshape(isinf(k_cv(:,predictable_signals,:)),[],1)) || ...
        any(reshape(isinf(predicted_signals(predictable_signals,:)),[],1)) ...
    error('Inf/NaN in kernel or predicted signal')
end

% Total explained variance (R^2)
sse_residual = sum((signals(predictable_signals,predictable_samples) - ...
    predicted_signals(predictable_signals,predictable_samples)).^2,2);
sse_total = sum((signals(predictable_signals,predictable_samples) - ...
    nanmean(signals(predictable_signals,predictable_samples),2)).^2,2);
explained_var.total = nan(size(signals,1),1);
explained_var.total(predictable_signals) = 1 - (sse_residual./sse_total);

% Partial (unique if exclusion) explained variance
if length(regressors) > 1    
    sse_residual_reduced = sum((signals(predictable_signals,predictable_samples) - ...
        predicted_signals_partial(predictable_signals,predictable_samples,:,:)).^2,2);
    sse_residual_full = sum((signals(predictable_signals,predictable_samples) - ...
        predicted_signals(predictable_signals,predictable_samples,:)).^2,2);   
    explained_var.partial = nan(size(signals,1),length(regressors),2);
    explained_var.partial(predictable_signals,:,:) = ...
        permute(1 - (sse_residual_full./sse_residual_reduced),[1,3,4,2]);
end

% Get the final k from averaging
% (to remove the constant term)
if ~return_constant
    k_vector = mean(k_cv(1:end-use_constant,:,:),3);
    % (to keep constant term)
elseif return_constant
    k_vector = mean(k_cv,3);
end

% Reshape kernel from vector to be regressor x t_shifts x signals
% (as cell array for each regressor/time shift pair)
if ~return_constant
    k = cellfun(@(x,t) reshape(x,[],length(t),size(signals,1)), ...
        mat2cell(k_vector,cellfun(@(x) ...
        size(x,1),regressors).*cellfun(@length,t_shifts),size(signals,1)), ...
        t_shifts','uni',false);
elseif return_constant
    k = cellfun(@(x,t) reshape(x,[],length(t),size(signals,1)), ...
        mat2cell(k_vector,[cellfun(@(x) size(x,1),regressors),1].*[cellfun(@length,t_shifts),1],size(signals,1)), ...
        [t_shifts,{1}]','uni',false);
end

% If kernel length is 1 (only 1 time shift) return as matrix
if length(k) == 1
    k = cell2mat(k);
end









