% Load in and plot behavior

animal = 'AP017';
% just get all days for now (eventually choose, supply date range, etc)
expInfo_path = ['\\zserver.cortexlab.net\Data\expInfo\' animal];
expInfo_dir = dir(expInfo_path);
days = {expInfo_dir(find([expInfo_dir(3:end).isdir])+2).name};

% Initialize the behavior structure
bhv = struct;

for curr_day = 1:length(days)
    
    day = days{curr_day};
    % In the event of multiple experiments, use the last one
    expDay_dir = dir([expInfo_path filesep days{curr_day}]);
    exp_nums = [expDay_dir(3:end).name];
    
    [block_filename, block_exists] = AP_cortexlab_filename(animal,day,exp_nums(end),'block');
   
    if ~block_exists
        continue
    end
    
    % Load the block file    
    load(block_filename)
    
    % (for now) Only use days with vanillaChoiceWorld
    [~,expDef] = fileparts(block.expDef);
    if ~strcmp(expDef,'vanillaChoiceworld');
        continue
    end
    
    % Time of session (in minutes)
    session_duration = block.duration/60;
    
    % Trial counts
    n_trials = length(block.paramsValues);
    total_water = sum(block.outputs.rewardValues);
           
    % Wheel movements/biases
    wheel_movement = diff(block.inputs.wheelValues);
    left_wheel_movement = abs(wheel_movement.*(wheel_movement < 0));
    right_wheel_movement = abs(wheel_movement.*(wheel_movement > 0));
    wheel_bias = (sum(right_wheel_movement)-sum(left_wheel_movement))/ ...
        (sum(right_wheel_movement)+sum(left_wheel_movement));
    
    % Performance (note that this excludes repeat on incorrect trials)
    performance = block.events.sessionPerformanceValues(:,end-10:end);  
        
    % Store in behavior structure   
    bhv.session_duration(curr_day) = session_duration;
    bhv.n_trials(curr_day) = n_trials;
    bhv.total_water(curr_day) = total_water;
    bhv.wheel_movement(curr_day) = sum(abs(wheel_movement));
    bhv.wheel_bias(curr_day) = wheel_bias;
    bhv.conditions = performance(1,:);
    bhv.n_trials_condition(curr_day,:) = performance(2,:);
    bhv.go_left_trials(curr_day,:) = performance(end,:);
    
    disp(['Behavior: day ' num2str(curr_day) '/' num2str(length(days))]);
    
end

% Plot summary
figure;
subplot(1,3,1);
yyaxis left
plot(1:length(days),bhv.n_trials./bhv.session_duration);
ylabel('Trials/min');
yyaxis right
plot(1:length(days),bhv.total_water);
ylabel('Total water');
xlabel('Session');

subplot(1,3,2);
yyaxis left
plot(1:length(days),bhv.wheel_movement./bhv.session_duration);
ylabel('Wheel movement / min');
yyaxis right
plot(1:length(days),bhv.wheel_bias);
ylim([-1,1]);
ylabel('Wheel bias');
xlabel('Session');

subplot(1,3,3);
imagesc(bhv.conditions,1:size(bhv.go_left_trials),bhv.go_left_trials./bhv.n_trials_condition);colormap(redblue)
color_step = diff(caxis)/size(colormap,1);
colormap([repmat(0.5,3);colormap]);
caxis([0-color_step*(255/size(colormap,1)),1]);
c = colorbar;
ylabel(c,'Go left (frac)');
xlabel('Condition');
ylabel('Session');










