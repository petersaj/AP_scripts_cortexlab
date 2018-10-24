function out = sparseNoiseTrackBW(state, ~, stimPersistence, stimulusRate, samplerFs)

out = max(abs(state)-1/samplerFs, 0).*sign(state); % subtract the time since last frame (but don't go negative)
newlyOnW = rand(size(state))<stimulusRate/samplerFs/2; % these squares will turn on next
% newlyOnB = rand(size(state))<stimulusRate/samplerFs/2; % these squares will turn on next
out(newlyOnW) = stimPersistence; % set the new squares to stimPersistence so they count down at the right time
% out(newlyOnB) = -stimPersistence; % set the new squares to stimPersistence so they count down at the right time

end