% Load in and plot behavior

animal = 'AP013';
% just get all days for now (eventually choose, supply date range, etc)
expInfo_dir = dir(['\\zserver.cortexlab.net\Data\expInfo\' animal]);
days = {expInfo_dir(find([expInfo_dir(3:end).isdir])+2).name};

% Initialize the behavior structure
bhv = struct;

for curr_day = 1:length(days)
    
    day = days{curr_day};
    [block_filename, block_exists] = AP_cortexlab_filename(animal,day,1,'block');
    
    % For now: don't use days before a certain date when relevant variables
    % were added to vanillaChoiceWorld
    if ~block_exists || datenum(day) - datenum('2017-03-31') < 0
        disp(['Skipping day ' num2str(curr_day)]);
        continue
    end
    
    % Load the block file    
    load(block_filename)
    
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
    left_choice = performance(end,:)./performance(2,:);
        
    % Store in behavior structure   
    bhv.n_trials(curr_day) = n_trials;
    bhv.total_water(curr_day) = total_water;
    bhv.wheel_movement(curr_day) = sum(abs(wheel_movement));
    bhv.wheel_bias(curr_day) = wheel_bias;
    bhv.conditions = performance(1,:);
    bhv.left_choice(curr_day,:) = left_choice;
    
    disp(['Behavior: day ' num2str(curr_day) '/' num2str(length(days))]);
    
end

% Plot summary
figure;
subplot(1,3,1);
yyaxis left
plot(1:length(days),bhv.n_trials);
ylabel('Number of trials');
yyaxis right
plot(1:length(days),bhv.total_water);
ylabel('Total water');
xlabel('Days');

subplot(1,3,2);
yyaxis left
plot(1:length(days),bhv.wheel_movement);
ylabel('Wheel movement');
yyaxis right
plot(1:length(days),bhv.wheel_bias);
ylim([-1,1]);
ylabel('Wheel bias');
xlabel('Days');


subplot(1,3,3);
imagesc(bhv.conditions,1:size(bhv.left_choice),bhv.left_choice);colormap(redblue)
color_step = diff(caxis)/size(colormap,1);
colormap([repmat(0.5,3);colormap]);
caxis([0-color_step*(255/size(colormap,1)),1]);
c = colorbar;
ylabel(c,'Left choice (frac)');
xlabel('Condition');
ylabel('Day');










