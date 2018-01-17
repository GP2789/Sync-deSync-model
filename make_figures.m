function [  ] = make_figures( filename, p )
% Create all of the graphs that are part of the Sync/deSync Model. Give a
% filename as a string, and a modifier for variable sampling & number of trials (0<p<=1).

    filename = [filename '_' num2str(p)];
    mkdir(filename);

    %% DEFAULT PARAMETERS MANY TRIALS - WEIGHT CHANGE / ACTIVITY / TFA GRAPHS
    trials = ceil(p^2 * 1000);
    vars = evaluate_network({'Default'}, trials, filename); % simulate Figures 2-6
    analyse_network(vars.test_f, {'FR', 'WC', 'TFA'}); % analyse data for Figures 2-6
    
    %% VARYING STIMULUS STRENGTH - LEARNING / POWER @ ENCODING V STIMULUS STRENGTH GRAPH
    x = round(logspace(0, 6, ceil(p*10)*5)); % create logorithmic scale
    trials = ceil(p * 20);
    vars = evaluate_network({'stimulus_strength',{x}}, trials, filename); % simulate Figure 7
    analyse_network(vars.test_f, {'TFA'}); % analyse data for Figure 7
    
    %% VARYING POSITIVE LEARNING RATE - LEARNING V POWER @ RECALL GRAPHS
    if(p>=0.66); r=0.1; elseif(p>=0.33); r=1/6; else; r = 1/3; end
    x = 0:r:1; % parameter scale
    trials = ceil(p * 100);
    vars = evaluate_network({'pos_LR',{x}}, trials, filename); % simulate Figure 8
    analyse_network(vars.test_f, {'TFA'}); % analyse data for Figure 8 
    
end

