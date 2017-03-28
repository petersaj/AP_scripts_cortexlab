function signals_tutorial_pong(t, events, pars, visStim, inputs, outputs, audio)
% 
% Setting up a simple version of pong is a great way to illustrate a lot of
% the common methods and strategies for setting up a Signals protocol. 
% 
% Create a wheel-controlled player paddle, a computer-controlled paddle
% which moves towards the ball, and a ball which moves around the screen
% and bounces off walls and paddles but resets when it gets missed by a
% paddle. 
%
% Some caution about a problem you'll run into:
% 
% This will bring up a really important issue of signals: truly
% co-dependent signals are impossible, because if signal1 and signal2 need
% information from each other, but you can only define one at a time, you
% cannot write 
% signal1 = f(signal2);
% signal2 = f(signal1); 
% because in the first line signal2 is not defined yet. 
% 
% You'll run into this when you try to make the computer paddle move
% towards the ball (computer paddle needs to know where the ball is to
% update it's position) while the ball needs to bounce off of the paddle
% (the ball needs to know where the paddle is to update it's position).
% This is solved by updating co-dependent parameters simultaneously using
% the scan function. In other words, try keeping all computer-controlled
% parameters in a structure which is intialized (seeded) with certain
% values and then updated via scan. Feed relevant player-data into this by
% having it be the source signal for the scan.









