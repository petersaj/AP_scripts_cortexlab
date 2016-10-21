function AP_save_sync_chan(data_filename, save_filename, n_chan, sync_chan)

% TO DO: this has to be locally loaded, takes way too long


% AP_save_sync_chan(data_filename, save_filename, numChans, sync_chan)
%
% Saves sync channel (always last channel?) from electrophysiology data
% Made from: saveExtraChans (NS)

disp(['Loading ' data_filename '...']);
thisDat = readDat(data_filename, n_chan, sync_chan);

disp(['Saving sync channel: ' save_filename '...']);
dat = thisDat;
save(save_filename, 'dat');

clear dat thisDat

disp(' Done')