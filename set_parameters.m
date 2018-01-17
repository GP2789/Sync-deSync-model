function [ ] = set_parameters( )
% Create set of parameters used for simulating network.

%% SET PARAMETER SPACE FUNCTION
    global par;
     if(isempty(par)==1)
        par = struct;
     end

    %% INDEPENDENT PARAMETERS
    parameter_settings = ...
        ... % NEURONS PROPERTIES
        {'V_th', -55, 'E', -70, 'V_m_initial', -60,'C_m', 240, 'tau_m', 20,... % intracellular
        'NC_V_ref', 2,  'NC_noise', 41800, 'Hip_V_ref', 2, 'Hip_noise', 3900,... % refractory period and noise
        'NC_amp', 21, 'NC_freq', 10, 'NC_adp_amp', 0, 'NC_adp_freq', 10, ... % Neo-Cortex Alpha
        'Hip_amp', 28, 'Hip_freq', 4, 'Hip_adp_amp', 100, 'Hip_adp_freq', 4,... % Hippocampal Theta
        'rand_phase', 'on' ... % random phase for oscillation cosine waves ('on' or 'off')
        ... % NEURON GROUPS PROPERTIES
        'n_Items', 2, 'NC_per_item', 10, 'Hip_per_item', 5,... % group properties
        'NC_Hip_syn_str', 5, 'NC_Hip_R_conn', 100, 'NC_Hip_delay', 2, 'NC_Hip_tau_syn', 5,... % NC -> Hip connectivity
        'Hip_NC_syn_str', 80, 'Hip_NC_R_conn', 100, 'Hip_NC_delay', 2, 'Hip_NC_tau_syn', 1.5,... % Hip -> NC connectivity
        'NC_NC_syn_str', 120, 'NC_NC_R_conn', 25, 'NC_NC_delay', 2, 'NC_NC_tau_syn', 1.5,... % NC connectivity
        'Hip_Hip_syn_str', 0,'Hip_Hip_R_conn', 40, 'Hip_Hip_delay', 2, 'Hip_Hip_tau_syn', 5,... % Hip connectivity
        ... % HIPPOCAMPAL LEARNING PROPERTIES
        'NC_STDP', 'no', 'Hip_STDP', 'yes', ... % enable STDP ('yes','no') for neuron groups
        'weight_change_type', 'DEFW', ... % change weight change type ('LINW','EXPW','SIGW', or 'DEFW' for default)
        'weight_max', 120, 'weight_decay', 20, 'pos_LR', 1, ... % max weights, decay & learning rates
        ... % STIMULUS PRESENTATION PROPERTIES
        'idling_length', 1000, 'pre_stim_length', 2000, 'stim_length', 2000,... % length of simulation periods (ms)
        'stim_TS', 250, 'stimulus_strength', 80000 }; % length and strength of input
   
    % SET INDEPENDENT par IF NOT ALREADY SUPPLIED
    for i=1:length(parameter_settings)/2
        if(any(strcmp(parameter_settings{(i-1)*2+1},fieldnames(par)))==0)
            par.([parameter_settings{(i-1)*2+1}]) = parameter_settings{(i-1)*2+2};
        end
    end
    
    %% DEPENDENT PARAMETERS   
    % GROUP SIZES & IDs
    par.n_Hip = par.Hip_per_item * par.n_Items;
    par.n_NC = par.NC_per_item * par.n_Items;
    par.network_size = par.n_NC + par.n_Hip;
    par.nG = {'NC', 'Hip'};
    
    % STDP LEARNING RATES
    par.learning_pos_max = par.weight_max; % STDP LTP increases: A+
    par.learning_neg_max = par.weight_max; % STDP LTD decreases: A-
    neg_LR = par.pos_LR*1.1; % A- / A+ = 1.1 (Song, Miller & Abbot STDP equation)
    par.learning_rate_pos = par.learning_pos_max * par.pos_LR;
    par.learning_rate_neg = par.learning_neg_max * neg_LR;
    
    % ORDER OF THE SIMULATION
    par.sim_order = {'idling','NP presentation','P presentation','learning','NP presentation','P presentation','idling'}; 
    par.sim_order_n = {'idle-BL','NP-BL','P-BL','DL','NP-AL','P-AL','idle-AL'}; % used for filenames

end



