function AP_run_kilosort(data_filename,input_board,combined_bank)
% AP_run_kilosort(data_filename,input_board,combined_bank)
%
% data_filename = .dat flat binary file of all channels together
% input_board = 'whisper' or 'oe'
% combined_bank = true/false (false by default): if the input to OE is from
% two headstages to one channel, the probe section order is flipped
% 
% Runs kilosort (modified from master_file, NS)

%% Set options
if ~exist('combined_bank','var') || isempty(combined_bank)
    combined_bank = false;
end

[data_path,data_file,data_ext] = fileparts(data_filename);


ops.Nfilt               = 512 ; %  number of filters to use (512, should be a multiple of 32)
ops.Nrank               = 3;    % matrix rank of spike template model (3)
ops.nfullpasses         = 6;    % number of complete passes through data during optimization (6)
ops.whitening           = 'full'; % type of whitening (default 'full', for 'noSpikes' set options for spike detection below)
ops.maxFR               = 20000;  % maximum number of spikes to extract per batch (20000)
ops.fs                  = 30000; % sampling rate %%%% CHANGED: 25000
ops.fshigh              = 300;   % frequency for high pass filtering
ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection
ops.scaleproc           = 200;   % int16 scaling of whitened data
ops.verbose             = 1;     
ops.nNeighPC            = 12; %12; % number of channnels to mask the PCs, leave empty to skip (12)
ops.nNeigh              = 16; % number of neighboring templates to retain projections of (16)
ops.NT                  = 32*1024+ ops.ntbuff;% this is the batch size, very important for memory reasons. 
% should be multiple of 32 (or higher power of 2) + ntbuff

% these options can improve/deteriorate results. when multiple values are 
% provided for an option, the first two are beginning and ending anneal values, 
% the third is the value used in the final pass. 
ops.Th               = [6 12 12];    % threshold for detecting spikes on template-filtered data ([6 12 12]) ([4,8,8])
ops.lam              = [10 30 30];   % large means amplitudes are forced around the mean ([10 30 30]) ([5 5 5])
ops.nannealpasses    = 4;            % should be less than nfullpasses (4)
ops.momentum         = 1./[20 400];  % start with high momentum and anneal (1./[20 1000])
ops.shuffle_clusters = 1;            % allow merges and splits during optimization (1)
ops.mergeT           = .1;           % upper threshold for merging (.1)
ops.splitT           = .1;           % lower threshold for splitting (.1)

ops.nNeighPC    = 12; %12; % number of channnels to mask the PCs, leave empty to skip (12)
ops.nNeigh      = 32; % number of neighboring templates to retain projections of (16)

% new options
ops.initialize = 'no'; %'fromData' or 'no';

% options for initializing spikes from data (NOT used during optimization)
ops.spkTh           = -6;      % spike threshold in standard deviations (4)
ops.loc_range       = [3  1];  % ranges to detect peaks; plus/minus in time and channel ([3 1])
ops.long_range      = [30  6]; % ranges to detect isolated peaks ([30 6])
ops.maskMaxChannels = 5;       % how many channels to mask up/down ([5])
ops.crit            = .65;     % upper criterion for discarding spike repeates (0.65)
ops.nFiltMax        = 10000;   % maximum "unique" spikes to consider (10000)
dd                  = load('PCspikes2.mat'); % you might want to recompute this from your own data
ops.wPCA            = dd.Wi(:,1:7);   % PCs 

ops.fracse  = 0.1; % binning step along discriminant axis for posthoc merges (in units of sd)
ops.epu     = Inf;

ops.showfigures    = 0;
ops.nSkipCov       = 10; % compute whitening matrix from every N-th batch
ops.whiteningRange = 32; % how many channels to whiten together (Inf for whole probe whitening, should be fine if Nchan<=32)

ops.ForceMaxRAMforDat   = Inf; % if you want to force it to use a swap file on disk and less RAM, or no RAM at all (set to 0). 

%% Copy raw data to local SSD temporarily (only if not already there)
disp('Copying data to local SSD...');
local_path = 'C:\Users\Andrew\Documents\CarandiniHarrisLab\data\kilosort_temp';
local_file = [data_file data_ext];
local_data_filename = [local_path filesep local_file];
if ~exist(local_data_filename)
    copyfile(data_filename,local_data_filename);
else
    disp('Already copied');
end

% data paths
root        = local_path;
fname       = [local_path filesep local_file];
fnameTW     = 'temp_wh.dat'; % will be created. residual of whitened data (not fitting in RAM).
% if chanMap is a string, it will load the channel configuration from a file 
% if chanMap is an array, that will indicate the channel order indexing into the raw data channel order
% This can also be used to drop inactive channels.
switch input_board
    case 'whisper'
        ops.chanMap = [local_path filesep 'forPRBimecToWhisper.mat']; %
    case 'oe'
        if combined_bank
            ops.chanMap = [local_path filesep 'forPRBimecToOE_combinedBank_fixed']; 
        else
            ops.chanMap = [local_path filesep 'forPRBimecToOE.mat']; 
        end
end

chanMap = load(ops.chanMap);
ops.NchanTOT            = length(chanMap.connected);   % total number of channels
ops.Nchan               = sum(chanMap.connected);   % number of active channels

%% Loads data into RAM + residual data on SSD and picks out spikes by a threshold for initialization
load_data_and_PCproject;

%% Do scaled kmeans to initialize the algorith,
if strcmp(ops.initialize, 'fromData')
    optimizePeaks;
end

%% Run kilosort
clear initialized
run_reg_mu2; % iterate the template matching (non-overlapping extraction)

fullMPMU; % extracts final spike times (overlapping extraction)

%% Convert kilosort output to npy

fprintf('Time %3.0fs. Done. saving results... \n', toc);

save([local_path filesep data_file '_ks_results'], 'rez', 'ops');

phy_dir = [local_path filesep 'phy'];
mkdir(phy_dir);
rezToPhy(rez, phy_dir);

fprintf('Time %3.0fs. Done. \n', toc);

%% Copy output files to original path (should be basket)

disp('Copying sorted data to original path...')

% phy folder
local_phy_path = [local_path filesep 'phy'];
copyfile(local_phy_path,data_path);

% ks_results file
ks_results_filename = [data_file '_ks_results.mat'];
local_ks_results_filename = [local_path filesep ks_results_filename];
copyfile(local_ks_results_filename,[data_path filesep ks_results_filename]);


%% Delete temporary data on local SSD

disp('Deleting temporary files...')

delete(local_data_filename);
delete([local_path filesep fnameTW]);
rmdir(local_phy_path,'s');
delete(local_ks_results_filename);

disp('Done.')




