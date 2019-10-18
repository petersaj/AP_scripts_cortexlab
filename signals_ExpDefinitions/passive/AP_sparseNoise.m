function AP_sparseNoise(t,events,parameters,visStim,inputs,outputs,audio)
% AP_sparseNoise: sparse noise for 1 sec, n trials = n seconds
% (same as / modified sparseNoiseAsync_NS2)

%% Parameters
sparseNoiseTime = 1;
stimPersistence = 1/6;
stimulusRate = 0.12;
samplerFs = 20; % slow enough too be detected during flicker screen

nAltitude = 10;
nAzimuth = 36;

%% Stim times

% Stimulus timer from trial start to stops
trial_running = events.newTrial.delay(0).to(events.newTrial.delay(sparseNoiseTime));

sampler = keepWhen(skipRepeats(floor(t*samplerFs)),trial_running); % to run a command at a certain sampling rate

stimuliTracker = sampler.scan(...
    @sparseNoiseTrackBW_internal, ...
    zeros(nAltitude, nAzimuth), 'pars', stimPersistence, stimulusRate, samplerFs);

% Default to black (-1) on t start, display stimulus on expStart
stimuliOn = merge( ...
    at(-1.*ones(nAltitude,nAzimuth),skipRepeats(t > 0)), ...
    skipRepeats((stimuliTracker>0)*2-1));

%% Generate stimuli
myNoise = vis.checker6(t);
myNoise.pattern = stimuliOn;

visStim.myNoise = myNoise;

%% End trial after alloted time
events.endTrial = events.newTrial.delay(sparseNoiseTime);

%% Save events
events.stimuliOn = stimuliOn;
events.stimuliTracker = stimuliTracker;
events.sampler = sampler;

end


% Function for creating stimulus
function out = sparseNoiseTrackBW_internal(state, ~, stimPersistence, stimulusRate, samplerFs)
out = max(abs(state)-1/samplerFs, 0).*sign(state); % subtract the time since last frame (but don't go negative)
newlyOnW = rand(size(state))<stimulusRate/samplerFs/2; % these squares will turn on next
out(newlyOnW) = stimPersistence; % set the new squares to stimPersistence so they count down at the right time
end

