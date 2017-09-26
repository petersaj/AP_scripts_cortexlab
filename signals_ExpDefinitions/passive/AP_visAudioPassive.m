function AP_visAudioPassive(t, events, parameters, visStim, inputs, outputs, audio)
% AP 2017-09-22: present random tones and visual stimuli



%% Set up stimuli

stimTime = 0.5;
itiTimeRange= [1:0.1:2];

% Sound stim
%audioFreqRange = linspace(5000,20000,5);
audioFreqRange = linspace(440,5000,3);
audioSampleRate = 192e3;

audioRampDuration = 0.01;
audioChannels = 2;

% 
% toneSamples = onsetToneAmplitude*events.expStart.map(@(x) ...
%     aud.pureTone(onsetToneFreq,stimTime,audioSampleRate, ...
%     rampDuration,audioChannels));
% 
% 
% 
% 
% 
% audio_t = [0:1/audioSampleRate:stimTime];
% audio_freq = linspace(audio_freq_range(1),audio_freq_range(2),length(audio_t));
% 
% rampSamples = round(audioSampleRate*rampDuration);
% audio_vol = [linspace(0,1,rampSamples),ones(1,length(audio_t)-2*rampSamples),linspace(1,0,rampSamples)];
% 
% audioSweepUp = sin(audio_t.*audio_freq.*2*pi).*audio_vol;
% audioSweepDown = fliplr(audioSweepUp);
% 
% audioStim = [audioSweepUp;audioSweepDown];
% 
% % Visual stim
% 
% vis_azimuths = [-60,0,60];
% spatialFrequency = 0.01;
% contrast = 1;
% vis_azimuth_sigma = [45,360];
% vis_elevation_sigma = [360,45];
% stimFlickerFrequency = 5;
% orientation = 45;


%% Stim times

n_audio_reps = 2;
n_visual_reps = 50;

audio_freq_uniform = repmat(audioFreqRange,n_audio_reps,1);
audio_freq_shuffle = audio_freq_uniform(randperm(numel(audio_freq_uniform)));


audio_itiTimes = randsample(itiTimeRange,n_audio_reps,true);
audio_startTimes = cumsum(audio_itiTimes + stimTime);




% visual_itiTimes = randsample(itiTimeRange,n_visual_reps,true);
% visual_startTimes = cumsum(visual_itiTimes + stimTime);
% visual_stopTimes = visual_startTimes + stimTime;


%% Present stim


%audioOnset = t.map(@(x) find(x > audio_startTimes,1,'last'));
%audioOnset = t.scan(@advanceStim,audio_startTimes);
%audioOnset = t.ge(audio_startTimes);
audioOnset = t.map(@(t) sum(t > audio_startTimes)).skipRepeats;


toneSamples = audioOnset.map(@(x) ...
    aud.pureTone(audio_freq_shuffle(x),stimTime,audioSampleRate, ...
    audioRampDuration,audioChannels));
audio.playTone = toneSamples;

% (TRYING TO INDEX FOR FREQUENCY BUT NOT WORKING AT THE MOMENT)


% 
% beepFreq = beepOnset.scan(@(x,y)rand()*(maxFreq-minFreq)+minFreq,0);
% 
% toneSamples = beepAmplitude*...
%   mapn(beepFreq, beepDur, audioSR, 0.02, @aud.pureTone);


endTrial = t.ge(audio_startTimes(end)).skipRepeats;

%% Events

events.t = t;
events.audioOnset = audioOnset;
events.endTrial = endTrial;



% %% Set up trial parameters and events
% 
% % Choose random stimulus on each trial
% trialAzimuth = events.newTrial.mapn(@(x) randsample(stimAzimuths,1));
% 
% % Turn stim on for fixed interval
% stimOn = at(true,events.newTrial); 
% stimOff = stimOn.delay(stimTime);
% 
% % Start the ITI after the stimulus turns off
% trialITI = events.newTrial.mapn(@(x) randsample(itiTimes,1));
% endTrial = stimOff.delay(trialITI);
% 
% %% Visual stimuli
% 
% stimFlicker = mod(skipRepeats(floor((t - t.at(stimOn))/(1/stimFlickerFrequency))),2);
% %stim.contrast = trialContrast.at(stimOn)*stimFlicker;
% 
% stim = vis.grating(t, 'square', 'gaussian');
% stim.sigma = sigma;
% stim.spatialFrequency = spatialFrequency;
% stim.phase = pi*stimFlicker;
% stim.azimuth = trialAzimuth.at(stimOn);
% stim.contrast = contrast;
% stim.orientation = orientation;
% stim.show = stimOn.to(stimOff);
% 
% visStim.stim = stim;
% 
% %% Events
% 
% events.stimOn = stimOn;
% events.stimOff = stimOff;
% events.stimFlicker = stimFlicker;
% events.trialAzimuth = trialAzimuth;
% 
% events.endTrial = endTrial;



end

function stimNum = advanceStim(stimOrder,t)
stimNum = find(t > stimOrder,1,'last');
end



















