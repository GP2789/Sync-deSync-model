function [  ] = recall_experiment( filename )
% Perform learning and recall experiment by creating and simulating
% networks specified by set_parameters. Give a filename for directory location (string).

%% DECLARATIONS
% initiate timers;
timer = tic;

% initiate global variables and set parameters
global weight_matrix; global weight_matrix_STDP; global par;
if(isempty(par)==1)
   set_parameters;
end

% make directories for simulation
if(exist(filename,'dir')~=7); mkdir(filename); end
if(exist([filename '/Data'],'dir')~=7); mkdir([filename '/Data']); end

h1 = waitbar(0, 'Recall Experiment', 'Units', 'normalized', 'Position', [0.5 0.4 0.2 0.1]);

%% MAIN SIMULATION

% random phase for static oscillations
par.Hip_r_phase_n = zeros(1, length(par.sim_order));
par.NC_r_phase_n = zeros(1, length(par.sim_order));

if(strcmp(par.rand_phase, 'on')==1)
    deg = 0:30:360; % 0-360 in 30deg intervals
    par.Hip_r_phase = round(deg(ceil(rand()*length(deg)))*(250/360))/1000; % choose random Theta phase
    par.NC_r_phase = round(deg(ceil(rand()*length(deg)))*(100/360))/1000; % choose random Alpha phase
else % set no random phases for static oscillations
    par.Hip_r_phase = 0; par.NC_r_phase = 0;
end

for n = 1:length(par.sim_order) % loop through each experiment phase
    % set sim_length and sim_type
    if(strcmp(par.sim_order{n},'idling')==1)
        par.sim_length = par.idling_length;
    else
        par.sim_length = par.pre_stim_length + par.stim_length;
    end
    par.sim_type = par.sim_order{n};
    waitbar(n/length(par.sim_order),h1);
    
    % create network
    create_network(); 
    par.Hip_r_phase_n(n) = par.Hip_r_phase; par.NC_r_phase_n(n) = par.NC_r_phase;
    
    if(n>1) % carry forward weight matrix from previous phase
        weight_matrix = repmat(carried_wm, [1 1 par.sim_length]); 
        weight_matrix_STDP = carried_wm_STDP;
    else % store weight matrix from first phase
        carried_wm_STDP = weight_matrix_STDP;
    end
    [data.sim_stats] = simulate_network(); % simulate network
    carried_wm = weight_matrix(:,:,end); 
    save([filename '/Data/' par.sim_order_n{n} '.mat'],'data')
end

close(h1);
display(['trial timer: ' int2str(toc(timer)) ' seconds']) % delete to suppress timer

end


