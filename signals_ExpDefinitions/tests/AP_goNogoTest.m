function AP_goNogoTest(t, events, parameters, visStim, inputs, outputs, audio)
% Testing go/nogo trial types

%% User parameters

startingAzimuth = 90;
responseDisplacement = 90;

quiescThreshold = 1;
noGoQuiescenceTime = 3;

itiTime = 0.5;

% (signals bug: some updates need small delay to not overlap/break)
signalsDelayTime = 0.1;


%% Initialize trial data

trialDataInit = events.expStart.mapn(@initializeTrialData).subscriptable;

%% Set up wheel 

wheel = inputs.wheelMM.skipRepeats();
wheelGain = 8; % deg/mm


%% Get new trial type at trial start

% Update performance
trialData = events.newTrial.scan(@updateTrialData,trialDataInit).subscriptable;

trialType = trialData.trialType;
goTrial = trialData.goTrial;
noGoTrial = trialData.noGoTrial;


%% Trial event times
% (this is set up to be independent of trial conditon, that way the trial
% condition can be chosen in a performance-dependent manner)

% Go trials:
% Stim onset
goStimOn = goTrial.delay(signalsDelayTime);
% Response
stimDisplacement = wheelGain*(wheel - wheel.at(goStimOn));
response = keepWhen(goStimOn.setTrigger(abs(stimDisplacement) ...
    >= responseDisplacement),goStimOn.to(events.newTrial));
% Stim offset
goStimOff = response.delay(signalsDelayTime);

% NoGo trials: 
% Stim onset
noGoStimOn = noGoTrial.delay(signalsDelayTime);
% Response
noGoQuiescenceTimeSig = at(noGoQuiescenceTime,noGoStimOn);
noGoQuiescence = sig.quiescenceWatch(noGoQuiescenceTimeSig, t, wheel, quiescThreshold); 
% Stim offset
noGoStimOff = noGoQuiescence.delay(signalsDelayTime);


% End trial at the end of the ITI
iti = merge(goStimOff.delay(itiTime),noGoStimOff.delay(itiTime));
endTrial = iti;


%% Visual stimulus

% Load images
allImgs{1} = imread('\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\noGoWorld\img1.jpeg');
allImgs{2} = imread('\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\noGoWorld\img2.jpeg');
allImgs{3} = imread('\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\noGoWorld\img3.jpeg');
trialStimImg =  cond( ...
    trialType == 1, allImgs{1}, ...
    trialType == 2, allImgs{2}, ...
    trialType == 3, allImgs{3});
    
goStimAzimuth = cond( ...
    goTrial.to(response), startingAzimuth + stimDisplacement, ...
    true, 0);

goStim = vis.image(t);
goStim.sourceImage = trialStimImg;
goStim.azimuth = goStimAzimuth;
goStim.show = goStimOn.to(goStimOff);

noGoStim = vis.image(t);
noGoStim.sourceImage = trialStimImg;
noGoStim.azimuth = 0;
noGoStim.show = noGoStimOn.to(noGoStimOff);

visStim.goStim = goStim;
visStim.noGoStim = noGoStim;

%% Display and save

events.trialType = trialType;
events.goTrial = goTrial;
events.noGoTrial = noGoTrial;

events.response = response;

events.goStimOn = goStimOn;
events.goStimOff = goStimOff;

events.noGoQuiescenceTimeSig = noGoQuiescenceTimeSig;
events.noGoQuiescence = noGoQuiescence;
events.noGoStimOn = noGoStimOn;
events.noGoStimOff = noGoStimOff;

events.endTrial = endTrial;


end

function trialDataInit = initializeTrialData(subject_info)
trialDataInit = struct;

trialDataInit.trialType = randi(3,1);
trialDataInit.goTrial = ismember(trialDataInit.trialType,[1,2]);
trialDataInit.noGoTrial = ismember(trialDataInit.trialType,[3]);
end

function trialData = updateTrialData(trialData,~)
trialData.trialType = randi(3,1);
trialData.goTrial = ismember(trialData.trialType,[1,2]);
trialData.noGoTrial = ismember(trialData.trialType,[3]);
end


