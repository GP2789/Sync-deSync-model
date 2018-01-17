function [ sim_stats ] = simulate_network()
% Simulate neural network defined in create_network and set_parameters. 

%% DECLARATIONS & INITIALISATIONS
% Declare global parameters
global par; global V_m; global I; global ADP; global osc; 
global MUA; global delay; global ref_period;
global weight_matrix; global weight_matrix_STDP; global wc_pos; global wc_neg;

h1 = waitbar(0, 'Simulating Network', 'Units', 'normalized', 'Position', [0.5 0.25 0.2 0.1]);

%% PRE-SIMULATION CALCULATIONS
% stores time since last spike for each neuron (used for ADP function)
last_spikes = ones(par.network_size,1); 

% defines Theta STDP multiplier, 1 for trough 0 for peak
theta = 1-(osc.Hip/max(osc.Hip) + 1)/2;

% find unique synaptic-time-constants in network and pre-calculate each EPSP
psp_u = unique(par.tau_syn); psp = cell(length(psp_u),1);
for i=1:length(psp_u)
    psp{i} = transpose(compute_psp(psp_u(i)));
end

% RECORD INPUT TO GROUPS OF NEURONS
I_REC.NC_NC = zeros(1, par.sim_length); I_REC.NC_HIP = zeros(1, par.sim_length);
I_REC.HIP_NC = zeros(1, par.sim_length); I_REC.HIP_HIP_WIT = zeros(1, par.sim_length); 
I_REC.HIP_HIP_BET = zeros(1, par.sim_length); I_REC.ADP_HIP = zeros(1, par.sim_length);
I_REC.BG_NC = mean(I.SYN(1:par.n_NC,:)) + mean(I.OSC(1:par.n_NC,:)); 
I_REC.BG_HIP = mean(I.SYN(par.n_NC+1:end,:)) + mean(I.OSC(par.n_NC+1:end,:));
 

%% MAIN SIMULATION
for t=2:par.sim_length
    waitbar(t/par.sim_length,h1)
    %% DETERMINE WEIGHT CHANGE TYPE
    % Define weight change multiplier based on current weight. Set to Default for JoN submission. 
    if(strfind(par.weight_change_type,'LINW')==1) % Linear
        weight_diff = ((par.weight_max - weight_matrix(:,:,t))./par.weight_max); 
    elseif(strfind(par.weight_change_type,'EXPW')==1) % Exponential
        weight_diff = exppdf(weight_matrix(:,:,t),par.weight_max/5)*par.weight_max/5;
    elseif(strfind(par.weight_change_type,'SIGW')==1) % Sigmoidal
        weight_diff = sigmf(weight_matrix(:,:,t),[-5/par.weight_max par.weight_max/2]);
    elseif(strfind(par.weight_change_type,'DEFW')==1) % Default - i.e. no modifier
        weight_diff = ones(size(weight_matrix(:,:,t)));
    end
    a_pos = weight_diff*par.learning_rate_pos; % define + synaptic rates of change
    a_neg = weight_diff*par.learning_rate_neg; % define - synaptic rates of change
    clear('weight_diff');

    %% WEIGHT DECAY 
    % Decay weight matrix exponentially to prune small weights
    wm = weight_matrix(:,:,t-1);
    decay = exppdf(weight_matrix(:,:,t-1),par.weight_decay).* weight_matrix_STDP .* (1-theta(t));
    wm = wm - decay; wm(wm<0)=0; % Piecewise linear function (no weights below 0)
    weight_matrix(:,:,t) = wm; clear('wm');
    
    % Decay potential +/- weight changes for STDP
    wc_pos(:,:,t) = round(wc_pos(:,:,t-1) - wc_pos(:,:,t-1)./(par.tau_syn)); % stores current potential + weight change
    wc_neg(:,:,t) = round(wc_neg(:,:,t-1) - wc_neg(:,:,t-1)./par.tau_syn);  % stores current potential - weight change
    
    %% HIPPOCAMPAL ADP ADDITION
    I.ADP(par.n_NC+1 : par.n_NC+par.n_Hip, t) = ADP.Hip(min((t) - ...
        last_spikes(par.n_NC+1 : par.n_NC+par.n_Hip),length(ADP.Hip)));
    I_REC.ADP_HIP(:,t) = mean(ADP.Hip(min((t) - last_spikes(par.n_NC+1 : par.n_NC+par.n_Hip),length(ADP.Hip))));
    
    %% VOLTAGE CHANGE
    % equation for integrate & fire voltage updates
    dv_dt = (par.E - V_m(:,t-1))/par.tau_m + (I.ADP(:,t) + I.SYN(:,t) + I.OSC(:,t))/par.C_m;
    % update voltages if neuron not in refractory
    V_m(:,t) = V_m(:,t-1) + dv_dt .* ref_period(:,t); 
    
    %% ADD SPIKES
    x = find(V_m(:,t) > -55); % find neurons over threshold
    MUA(x,t) = MUA(x,t) + 1; last_spikes(x) = t;  % add spike to MUA and update last spike time
    ref_period(MUA(:,t)==1, t+1:min(par.sim_length,t+par.V_ref)) = 0; % update refractory period (clamp)
    V_m(MUA(:,t)==1, t) = par.E; % reset voltage to E
    
    x = find(MUA(:,t)>0); % find spike events
    wm = weight_matrix(x,:,t); [c, y]=find(wm>0); clear('wm'); % find active synapses at spiking neurons
    for i = 1:length(c) % add spikes to input of connected neurons
        % find relevant EPSP based on synaptic time constant of spiking neurons
        epsp = psp{find(psp_u==par.tau_syn(x(c(i)),y(i)))}*weight_matrix(x(c(i)),y(i),t);
        del = t + delay(x(c(i)),y(i)); % find delay between relevant neurons
        % add EPSP to future input of downstream neurons
        I.SYN(y(i), del:min((del-1)+length(epsp),par.sim_length)) = ...
            I.SYN(y(i), del:min((del-1)+length(epsp),par.sim_length)) + ...
            epsp(1:min(length(epsp),par.sim_length-del+1)); 
        
        % RECORD INPUT STREAMS FOR GROUPS OF NEURONS (HIP & NC) - only necessary for Figure 3Cii
        if(x(c(i)) <= par.n_NC && y(i) <= par.n_NC) % NC -> NC
            I_REC.NC_NC(:, del:min((del-1)+length(epsp),par.sim_length)) = ...
                I_REC.NC_NC(:, del:min((del-1)+length(epsp),par.sim_length)) + ...
                epsp(1:min(length(epsp),par.sim_length-del+1));
        elseif(x(c(i)) <= par.n_NC && y(i) > par.n_NC) % NC -> HIP
            I_REC.NC_HIP(:, del:min((del-1)+length(epsp),par.sim_length)) = ...
                I_REC.NC_HIP(:, del:min((del-1)+length(epsp),par.sim_length)) + ...
                epsp(1:min(length(epsp),par.sim_length-del+1));
        elseif(x(c(i)) > par.n_NC && y(i) <= par.n_NC) % HIP -> NC
            I_REC.HIP_NC(:, del:min((del-1)+length(epsp),par.sim_length)) = ...
                I_REC.HIP_NC(:, del:min((del-1)+length(epsp),par.sim_length)) + ...
                epsp(1:min(length(epsp),par.sim_length-del+1));
        elseif(x(c(i)) > par.n_NC && y(i) > par.n_NC) % HIP -> HIP
            if(x(c(i)) <= par.n_NC + par.n_Hip/2 && y(i) <= par.n_NC + par.n_Hip/2 || ...
                    x(c(i)) > par.n_NC + par.n_Hip/2 && y(i) > par.n_NC + par.n_Hip/2) % WITHIN GROUPS
                I_REC.HIP_HIP_WIT(:, del:min((del-1)+length(epsp),par.sim_length)) = ...
                    I_REC.HIP_HIP_WIT(:, del:min((del-1)+length(epsp),par.sim_length)) + ...
                    epsp(1:min(length(epsp),par.sim_length-del+1));
            elseif(x(c(i)) <= par.n_NC + par.n_Hip/2 && y(i) > par.n_NC + par.n_Hip/2 || ...
                    x(c(i)) > par.n_NC + par.n_Hip/2 && y(i) <= par.n_NC + par.n_Hip/2) % BETWEEN GROUPS
                I_REC.HIP_HIP_BET(:, del:min((del-1)+length(epsp),par.sim_length)) = ...
                    I_REC.HIP_HIP_BET(:, del:min((del-1)+length(epsp),par.sim_length)) + ...
                    epsp(1:min(length(epsp),par.sim_length-del+1));
            end
        end
        clear('epsp');
    end
    
    %% STDP LEARNING
    if(sum(MUA(:,t))>0) 
        x = find(MUA(:,t)==1); % find spike events
        % add to potential weight changes
        wc_pos(x,:,t) = wc_pos(x,:,t) + a_pos(x,:).*weight_matrix_STDP(x,:); % potential + weight change
        wc_neg(:,x,t) = wc_neg(:,x,t) - a_neg(:,x).*weight_matrix_STDP(:,x); % potential - weight change
        % modified spike-time-dependent-plasticity rule
        weight_matrix(:,x,t) = weight_matrix(:,x,t) + wc_pos(:,x,t) * theta(t); % add to contributing weights
        weight_matrix(x,:,t) = weight_matrix(x,:,t) + wc_neg(x,:,t) * theta(t); % take from competing weights
        % piecewise linear bounding function (no weights below 0 or above max weight)
        wm = weight_matrix(:,:,t); wm(wm<0) = 0; wm(wm>par.weight_max) = par.weight_max;
        weight_matrix(:,:,t) = wm; clear('wm');
    end
end

%% SAVE SIMULATION DATA
% SPIKE DETECTOR
[s, d] = find(MUA>0); sim_stats.spike_detector = [s d];
ID = transpose(1:par.network_size); 
% VOLTMETER
sim_stats.voltmeter = zeros(par.sim_length*par.network_size, 3);
for t=1:par.sim_length
   sim_stats.voltmeter((t-1)*length(ID)+1:t*length(ID),:) = [ID t*ones(par.network_size,1) V_m(:,t)];
end
% INPUT RECORDING STREAMS
I_REC.NC_NC = I_REC.NC_NC / par.n_NC; I_REC.HIP_NC = I_REC.HIP_NC / par.n_NC; 
I_REC.HIP_HIP_WIT = I_REC.HIP_HIP_WIT / par.n_Hip; I_REC.HIP_HIP_BET = I_REC.HIP_HIP_BET / par.n_Hip; 
I_REC.NC_HIP = I_REC.NC_HIP / par.n_Hip; 
sim_stats.I_REC = I_REC; clear('I_REC');
% REFRACTORY PERIOD
sim_stats.RP = ref_period; clear('ref_period');
% WEIGHT MATRIX AND ACTIVE SYNAPSES
sim_stats.weight_matrix = weight_matrix; sim_stats.weight_matrix_STDP = weight_matrix_STDP;

close(h1)

end

