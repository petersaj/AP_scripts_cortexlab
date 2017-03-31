function parsStruct = inferParameters(expdef)
%EXP.INFERPARAMETERS Infers the parameters required for experiment
%   Detailed explanation goes here

% create some signals just to pass to the definition function and track
% which parameter names are used

if ischar(expdef) && file.exists(expdef)
  expdeffun = fileFunction(expdef);
else
  expdeffun = expdef;
  expdef = which(func2str(expdef));
end

net = sig.Net;
e = struct;
e.t = net.origin('t');
e.events = net.subscriptableOrigin('events');
e.pars = net.subscriptableOrigin('pars');
e.pars.CacheSubscripts = true;
e.visual = net.subscriptableOrigin('visual');
e.audio = net.subscriptableOrigin('audio');
e.inputs = net.subscriptableOrigin('inputs');
e.outputs = net.subscriptableOrigin('outputs');

try
  expdeffun(e.t, e.events, e.pars, e.visual, e.inputs , e.outputs);
  paramNames = e.pars.Subscripts.keys';
  parsStruct = cell2struct(cell(size(paramNames)), paramNames);
  parsStruct.numRepeats = 0; % add 'numRepeats' parameter
  %%%% AP 2017-03-30 %%%%
  % Define the ExpPanel to use (automatically by name convention for now)
  [protocol_path,prototol_filename,protocol_extension] = fileparts(expdef);
  ExpPanel_name = [prototol_filename '_ExpPanel'];
  ExpPanel_fn = [protocol_path filesep ExpPanel_name protocol_extension];
  if exist(ExpPanel_fn,'file')
      cd(lower(protocol_path))
      parsStruct.expPanelFun = str2func(ExpPanel_name);
  end
  %%%%%%%%%%%%%%%%%%%%%%%
  parsStruct.defFunction = expdef;
  parsStruct.type = 'custom';
catch ex
  net.delete();
  rethrow(ex)
end

net.delete();


end