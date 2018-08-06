%% Analysis for Elina's 3-5Hz oscillation paper

%% Bin spikes to get MUA

% Discretize spikes into frames and count spikes per frame
use_spikes = spikeDepths > 0 & spikeDepths < 3500;

frame_edges = [frame_t,frame_t(end)+1/framerate];
frame_spikes = histcounts(spike_times_timeline(use_spikes),frame_edges);


%% Get fluorescence around probe

fluor_trace = AP_svd_roi(Udf,fVdf,avg_im);

%% Plot MUA and fluorescence together

figure; hold on
plot(frame_t,zscore(frame_spikes),'k','linewidth',1);
plot(frame_t,zscore(fluor_trace),'r','linewidth',1);
xlabel('Time');
ylabel('Z-scored activity');
legend({'MUA','Fluorescence'});


%% MUA/fluorescence cross-correlation and coherence

% Cross-correlation
corr_lag = 5; % in seconds
corr_lag_samples = round(corr_lag*framerate);
[mua_fluor_xcorr,lags] = xcorr(fluor_trace,frame_spikes,corr_lag_samples);
lags_t = lags./framerate;

% Coherence
[mua_fluor_coherence,f] = mscohere(frame_spikes',fluor_trace', ...
    hanning(round(framerate*10)),round(framerate*5),[],framerate);

% Plot
figure;

subplot(2,1,1);
plot(lags_t,mua_fluor_xcorr,'k','linewidth',2)
line([0,0],ylim,'linestyle','--','color','k');
xlabel('MUA Lag (s)');
ylabel('Cross-correlation')

subplot(2,1,2);
plot(f,mua_fluor_coherence,'k','linewidth',2)
xlabel('Freqency');
ylabel('Magnitude-squared coherence');


%% MUA/fluorescence spectral correlation

% Spectrogram settings
spect_overlap = 80;
window_length = 3; % in seconds
window_length_samples = window_length/(1/Fs);
N = window_length_samples; % window length
df = Fs/N; % frequency increment

% Time of traces
use_t = frame_t;
Fs = 1./median(diff(use_t));

% Power spectrum of MUA
use_trace = frame_spikes(1:end-1);

[s,f,t] = spectrogram(use_trace,window_length_samples, ...
    round(spect_overlap/100*window_length_samples),[],Fs);
s_squared = (s/Fs).*conj(s/Fs);  % Fs is used to normalize the FFT amplitudes
mua_power = s_squared*2*df;

% Power spectrum of fluorescence
use_trace = fluor_trace;

[s,f,t] = spectrogram(use_trace,window_length_samples, ...
    round(spect_overlap/100*window_length_samples),[],Fs);
s_squared = (s/Fs).*conj(s/Fs);  % Fs is used to normalize the FFT amplitudes
fluor_power = s_squared*2*df;

% Correlate power spectra of MUA and fluorescence
spectra_corr = mat2cell(corrcoef([mua_power',fluor_power']),...
    repmat(length(f),1,2),repmat(length(f),1,2));

figure; colormap(hot);
c = [min(reshape(cell2mat(spectra_corr),[],1)), ...
    max(reshape(cell2mat(spectra_corr),[],1))];

subplot(1,3,1); 
imagesc(f,f,spectra_corr{1,1});
xlabel('MUA frequency');
ylabel('MUA frequency');
axis square; caxis(c);

subplot(1,3,2); 
imagesc(f,f,spectra_corr{2,2});
xlabel('Fluor frequency');
ylabel('Fluor frequency');
axis square; caxis(c);

subplot(1,3,3); 
imagesc(f,f,spectra_corr{1,2});
xlabel('FLuor frequency');
ylabel('MUA frequency');
axis square; caxis(c);



%% ~~~~~~~~~~ OLD CODE TO DRAW FROM ~~~~~~~~~~


%% Spectral analysis
% Power spectrum
use_trace = r;
use_t = frame_t(1:end-1);

Fs = 1./median(diff(use_t));
L = length(use_trace);
NFFT = 2^nextpow2(L);
[P,F] = pwelch(double(use_trace)',[],[],NFFT,Fs);
Pc = smooth(P,50); 
figure;plot(F,log10(Pc),'k')
xlabel('Frequency');
ylabel('Log Power');

% Notch filter
freqs = [3 5];
for ff = 1:length(freqs)
    f = fdesign.notch('N,F0,Q',2,freqs(ff),10,Fs);
    h = design(f);
    % hfvt= fvtool(h,'Color','white');
    use_trace_filt = filter(h, use_trace')';
end

p = bandpower(use_trace,Fs,[3,5]);

% Spectrogram
spect_overlap = 50;
window_length = 5; % in seconds
window_length_samples = window_length/(1/Fs);
figure;spectrogram(use_trace,window_length_samples, ...
    round(spect_overlap/100*window_length_samples),[],Fs,'yaxis')
colormap(hot)

% Band power over time
spect_overlap = 80;
window_length = 3; % in seconds
window_length_samples = window_length/(1/Fs);

[s,f,t] = spectrogram(use_trace,window_length_samples, ...
    round(spect_overlap/100*window_length_samples),[],Fs);

N = window_length_samples; % window length
df = Fs/N; % frequency increment
s_squared = (s/Fs).*conj(s/Fs);  % Fs is used to normalize the FFT amplitudes
power_0_2 = 2*sum(s_squared(f >= 0.1 & f <= 2,:))*df; 
power_3_6 = 2*sum(s_squared(f >= 3 & f <= 6,:))*df; 
power_10_14 = 2*sum(s_squared(f >= 10 & f <= 14,:))*df; 

figure; hold on;
plot(power_0_2,'k');
plot(power_3_6,'b');
plot(power_10_14,'r');


%% Get kernel between spikes and fluorescence

% Skip the first n seconds to do this
skip_seconds = 10;
use_frames = frame_t > skip_seconds;

% Get fluorescence in ROI
roi_trace_full = AP_svd_roi(Udf,fVdf,avg_im);
roi_trace = roi_trace_full(use_frames);

% Get population spikes per frame
framerate = 1./nanmedian(diff(frame_t));
frame_edges = [frame_t,frame_t(end)+1/framerate];

use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > 0 & templateDepths < 1500)));

[frame_spikes_full,~,spike_frames] = histcounts(use_spikes,frame_edges);
frame_spikes = frame_spikes_full(use_frames);

corr_lags = 100;
[~,lags] = xcorr(ones(size(frame_spikes)),corr_lags);
lags_t = lags./framerate;
figure;
subplot(2,1,1); hold on;
plot(lags_t,xcov(frame_spikes,corr_lags,'coeff'),'k','linewidth',2)
plot(lags_t,xcov(roi_trace,corr_lags,'coeff'),'r','linewidth',2)
line([0,0],ylim,'linestyle','--','color','k');
legend({'Spikes','Fluorescence'});
xlabel('Lag (s)')
ylabel('Autocorrelation')
subplot(2,1,2);
plot(lags_t,xcov(roi_trace,frame_spikes,corr_lags),'k','linewidth',2)
line([0,0],ylim,'linestyle','--','color','k');
xlabel('Lag (s)');
ylabel('Cross-correlation')

% Non-normalized
x_nonorm = ifft((fft(roi_trace).*conj(fft(frame_spikes)))); % unnormalized

% I think this is what Krumin suggested? (not sure if he meant mean)
x = ifft((fft(roi_trace).*conj(fft(frame_spikes)))./(mean(fft(frame_spikes)).*mean(conj(fft(frame_spikes)))));

% This looks like from Nauhaus 2012?
soft_reg_factor = 1e6;
x_autonorm = ifft((fft(roi_trace).*conj(fft(frame_spikes)))./(soft_reg_factor+fft(frame_spikes).*conj(fft(frame_spikes))));

plot_frames = 35*5;

figure;

t_shift = [frame_t(end-plot_frames+1:end)-frame_t(end)-1/framerate,frame_t(1:plot_frames)-frame_t(1)];

p1 = subplot(2,1,1);
plot(t_shift,[x_nonorm(end-plot_frames+1:end),x_nonorm(1:plot_frames)],'k','linewidth',2);
xlabel('Time (s)');
ylabel('Impulse response')
title('Non-normalized: ifft(F*S'')');

p2 = subplot(2,1,2);
plot(t_shift,[x_autonorm(end-plot_frames+1:end),x_autonorm(1:plot_frames)],'k','linewidth',2);
xlabel('Time (s)');
ylabel('Impluse response');
title('Normalized: ifft(F*S''/S*S'')');

linkaxes([p1,p2],'x')

% Use the corrected impulse response for convolving kernel
gcamp_kernel = x_autonorm(1:plot_frames);
frame_spikes_conv_full = conv(frame_spikes_full,gcamp_kernel);
frame_spikes_conv = frame_spikes_conv_full(1:length(frame_spikes_full));

figure;hold on;
plot(frame_t,zscore(roi_trace_full),'k','linewidth',1);
plot(frame_t,zscore(frame_spikes_conv),'b','linewidth',1);
xlabel('Time (s)');
legend({'Fluorescence','Spikes conv'});



%% Coherence

frame_spikes = fake_frame_spikes;

skip_seconds = 4;
use_frames = (frame_t > skip_seconds);

% Downsample and get pixels to compare across
U_downsample_factor = 10;
Ud = imresize(U,1/U_downsample_factor,'bilinear');
px = reshape(svdFrameReconstruct(Ud,fV),[],size(fV,2));

% Coherence
[px_coherence,f] = mscohere(px(:,use_frames)',frame_spikes(use_frames)',hanning(round(framerate*10)),round(framerate*5),[],framerate);
coherence_map = reshape(sum(px_coherence,1),size(Ud,1),size(Ud,2));











