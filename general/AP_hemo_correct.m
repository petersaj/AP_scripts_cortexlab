function [V_corrected, T] = AP_hemo_correct(U, V, Vaux, fS, FreqRange, pixSpace)

% TO DO: low-pass filter V? 
fV = V;
fVaux = Vaux;

Ud = imresize(U,1/pixSpace,'bilinear');
Usub = reshape(Ud,[],size(U,3));

pixTrace = Usub*fV;
pixAux = Usub*fVaux;

n_px = size(pixTrace,1);
n_frames = size(pixTrace,2);
tic
pixFit = nan(n_px,2);
pixCorrected = pixTrace;
for curr_px = 1:n_px
   curr_fit = [pixAux(curr_px,:)',ones(n_frames,1)]\pixTrace(curr_px,:)';
   curr_model = [pixAux(curr_px,:)',ones(n_frames,1)]*curr_fit;
   
   pixFit(curr_px,:) = curr_fit;
   pixCorrected(curr_px,:) = pixCorrected(curr_px,:) - curr_model';
   
   AP_print_progress_fraction(curr_px,n_px);
end
toc

% Put fits into V space
V_add = pixFit(:,2)'/Usub';
V_mult = pinv(Usub)*diag(pixFit(:,1))*Usub;
V_corrected = fV - bsxfun(@minus,fVaux*V_mult,V_add);






