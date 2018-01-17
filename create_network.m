function [  ] = create_network( )
% Create and connect neural network based on parameters specified in set_parameters.

%% DECLARE GLOBAL VARIABLES
global par; global V_m; global I; global ADP; global osc; 
global MUA; global delay; global ref_period; 
global weight_matrix; global weight_matrix_STDP; global wc_pos; global wc_neg;

%% INITIALISE DATA STRUCTURES
% neuron data storage
V_m = zeros(par.network_size, par.sim_length); % keeps track of voltages through time
par.V_ref = zeros(par.network_size, 1); % stores refractory period data for neurons
par.tau_syn = ones(par.network_size, par.network_size); % stores synaptic time constants
ref_period = ones(par.network_size, par.sim_length); % tracks current refractory period status through time
I.SYN = zeros(par.network_size, par.sim_length); % input for spikes
I.OSC = zeros(par.network_size, par.sim_length); % input for static oscillations
I.ADP = zeros(par.network_size, par.sim_length); % input for ADP functions
MUA = zeros(par.network_size, par.sim_length); % tracks spike events through time

% synapse data storage
weight_matrix = zeros(par.network_size, par.network_size, par.sim_length); % tracks weights through time
wc_pos = zeros(par.network_size, par.network_size, par.sim_length); % tracks decaying + weight change through time
wc_neg = zeros(par.network_size, par.network_size, par.sim_length); % tracks decaying - weight change through time
weight_matrix_STDP = zeros(par.network_size, par.network_size); % informs on whether synapses have learning enabled (binary)
delay = zeros(par.network_size, par.network_size); % stores delay information on synapses

% set synaptic time constants
par.tau_syn(1:par.n_NC,1:par.n_NC) = par.NC_NC_tau_syn;
par.tau_syn(1:par.n_NC,1+par.n_NC:par.network_size) = par.NC_Hip_tau_syn;
par.tau_syn(1+par.n_NC:par.network_size,1:par.n_NC) = par.Hip_NC_tau_syn;
par.tau_syn(1+par.n_NC:par.network_size,1+par.n_NC:par.network_size) = par.Hip_Hip_tau_syn;

V_m(:,1) = par.V_m_initial; % initial voltage starting value

t = 0:1/1000:par.sim_length/1000-0.001;
c = 0; nG = par.nG;

%% CREATE AND CONNECT NEURONS
% loops through neuron groups to extract and set intrinsic and network properties defined in set_parameters.m
for g=1:length(nG)
    par.V_ref(c+1:c+par.(['n_' nG{g}])) = par.([nG{g} '_V_ref']); % group refractory periods
    %% create & connect oscillation sine wave
    if(strcmp(par.rand_phase,'on')==1) % random phase option
        par.([nG{g} '_r_phase']) =  round(rand()*(1000/par.([nG{g} '_freq'])))/1000;
    else
        par.([nG{g} '_r_phase']) = 0;
    end
    % create static oscillations
    osc.(nG{g})  = cos(2*pi*par.([nG{g} '_freq'])*(t+par.([nG{g} '_r_phase'])))*par.([nG{g} '_amp']);
    % add oscillation to I input for selected neurons
    I.OSC(c+1:c+par.(['n_' nG{g}]),:) = repmat(osc.(nG{g}), [par.(['n_' nG{g}]) 1]);
    
    %% create ADP function: equation from Jensen/Lisman/Idiart 96
    adp = 1:1:(1000/par.([nG{g} '_adp_freq']));
    adp =  adp.*exp(1-adp/(1000/par.([nG{g} '_adp_freq'])));
    ADP.(nG{g}) = (adp/length(adp)) * par.([nG{g} '_adp_amp']);
    
    %% create & connect background noise
    % create EPSP for spike events
    psp = transpose(compute_psp(par.([nG{g} '_' nG{g} '_tau_syn'])));
    for i=c+1:c+par.(['n_' nG{g}])
        % create new spike train for each neuron
       spikes = spike_train(par.([nG{g} '_noise']),par.sim_length,par.([nG{g} '_noise'])/100);
       % add EPSP * spike events to I input of selected neurons
       for s = 1:length(spikes)
          I.SYN(i,s:min((s-1)+length(psp),par.sim_length)) = I.SYN(i,s:min((s-1)+length(psp),par.sim_length)) + ...
              spikes(s)*psp(1:min(length(psp),par.sim_length-s+1)); 
       end
    end
    
    %% connect network groups together
    for g2 = 1:length(nG)
         c1 = c;
        for i = 1:par.n_Items
            c2 = (g2-1)*par.(['n_' nG{1}]);
            for i2 = 1:par.n_Items
                % find X->Y neuron IDs
                from = c1+1:c1+par.(['n_' nG{g}])/par.n_Items; 
                to = c2+1:c2+par.(['n_' nG{g2}])/par.n_Items;
                % set off/on synapses based on X->Y connectivity parameters
                syn = double(par.([nG{g} '_' nG{g2} '_R_conn'])/100 > rand(length(from),length(to)));
                str = par.([nG{g} '_' nG{g2} '_syn_str']);
                % set up STDP matrix to distinguish synapses with enabled learning
                if(strcmp(par.([nG{g} '_STDP']),'yes')==1 && strcmp(par.([nG{g2} '_STDP']),'yes')==1) 
                    weight_matrix_STDP(from,to) = syn;
                    % set pre-defined routes - half of active synapses to max weight, half to 0
                    if(i == i2)
                        x = find(syn==1); % find active synapses
                        x2 = zeros(floor(length(x)/2),1);
                        for xi = 1:length(x2)
                           xr = ceil(rand()*length(x));
                           x2(xi) =  x(xr); x(xr) = [];
                        end
                        str = zeros(size(syn)); str(x2) = par.weight_max;
                    end 
                end
                % input weighted synapses and delays for X->Y
                if(i == i2); weight_matrix(from,to,:) = repmat(syn .* str,[1 1 par.sim_length]); end
                delay(from,to) = syn * par.([nG{g} '_' nG{g2} '_delay']); 
                c2 = c2 + par.(['n_' nG{g2}])/par.n_Items; % counting variable
            end
            c1 = c1 + par.(['n_' nG{g}])/par.n_Items; % counting variable
        end
    end
    c = c + par.(['n_' nG{g}]); % counting variable
end
%% no self connections
for i=1:size(weight_matrix,1); for j=1:size(weight_matrix,2); if(i==j); weight_matrix(i,j,:) = 0; weight_matrix_STDP(i,j) = 0; delay(i,j) = 0;end; end; end 

%% create & connect stimulus input
if(isempty(strfind(par.sim_type,'idling'))==1)
    % create stimulus alpha function
    T = 1:1:par.stim_length; spikes = zeros(1,par.sim_length);
    stim = [zeros(1, par.pre_stim_length) (exp(1)*T./par.stim_TS).*exp(-T./par.stim_TS).*heaviside(T)];
    % create spike train based on parameters
    spikes(par.pre_stim_length+1:end) = ...
        spike_train(par.stimulus_strength ,par.stim_length, ceil(par.stimulus_strength/100));
    spikes = spikes .* stim; % multiply generated spike events by alpha function
    
    if(strfind(par.sim_type,'P')==1) % P stimulation
        n = 1 : par.NC_per_item;
    elseif(strfind(par.sim_type,'NP')==1) % NP stimulation
        n = par.NC_per_item+1 : par.NC_per_item*2;
    else
        n = 1:par.n_NC; % learning (P & NP stimulation)
    end
    % add generated spike events to future input
    for s = 1:length(spikes)
      I.SYN(n,s:min((s-1)+length(psp),par.sim_length)) = I.SYN(n,s:min((s-1)+length(psp),par.sim_length)) + ...
          repmat(spikes(s)*psp(1:min(length(psp),par.sim_length-s+1)),[length(n) 1]); 
    end
end

end

