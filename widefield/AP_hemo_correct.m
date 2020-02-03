function [V_corrected, T] = AP_hemo_correct(U, V, Vaux, framerate, pixSpace)

% Subtract out means
zV = bsxfun(@minus,V,mean(V,2));
zVaux = bsxfun(@minus,Vaux,mean(Vaux,2));

zV = V;
zVaux = Vaux;

% Low pass filter to get rid of drift
highpassCutoff = 0.01;
[b, a] = butter(2, highpassCutoff/(framerate/2), 'high');
V_highpass = single(filtfilt(b,a,double(zV)')');
Vaux_highpass = single(filtfilt(b,a,double(zVaux)')');

% High pass filter to not fit high-frequency artifacts
lowpassCutoff = 10;
[b, a] = butter(2, lowpassCutoff/(framerate/2), 'low');
V_bandpass = single(filtfilt(b,a,double(V_highpass)')');
Vaux_bandpass = single(filtfilt(b,a,double(Vaux_highpass)')');

% Downsample spatial components to speed processing
Ud = imresize(U,1/pixSpace,'bilinear');
Usub = reshape(Ud,[],size(U,3));

pixTrace = Usub*V_bandpass;
pixAux = Usub*Vaux_bandpass;

n_px = size(pixTrace,1);
n_frames = size(pixTrace,2);

pixFit = nan(n_px,2);
pixCorrected = pixTrace;
for curr_px = 1:n_px
   curr_fit = [pixAux(curr_px,:)',ones(n_frames,1)]\pixTrace(curr_px,:)';
   curr_model = [pixAux(curr_px,:)',ones(n_frames,1)]*curr_fit;
   
   pixFit(curr_px,:) = curr_fit;
   pixCorrected(curr_px,:) = pixCorrected(curr_px,:) - curr_model';
   
   AP_print_progress_fraction(curr_px,n_px);
end

% Put fits into V space
V_add = pixFit(:,2)'/Usub';
V_mult = pinv(Usub)*diag(pixFit(:,1))*Usub;
V_corrected = V_highpass - bsxfun(@plus,Vaux_highpass'*V_mult,V_add)';




% turns out the scale and multiply doesn't work when done on filtered and
% applied to raw, maybe just filter the low freq stuff? that's not really
% kosher though. maybe this is all crap because somehow now kenneth's is
% looking like it's working on this data set even though it didn't
% before... 



