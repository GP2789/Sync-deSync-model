function [ vars ] = evaluate_network( variables, trials, filename )
%% Make multiple simulations to assess up to 3 different parameters:
% evaluate_network({par1 name, {[par1vars]}}, {par2 name, {[par2 vars]}}, {par3
% name, {[par3 vars]}}, trials). To run default parameters with multiple
% trials run evaluate_network({'Default'}, trials). See set_parameters.m to check parameter
% options.

%% SET PARAMETER SPACE
clearvars -global par; global par; set_parameters;

%% INITIALISE DIRECTORY AND VARIABLES FROM USER INPUT
vars.trials = trials; vars.test_f = [filename '/'];
for i = 1:3
    vars.(['v' int2str(i) '_name']) = []; % set var name to empty if no variable presented
    if(length(variables) > (i-1)*2 && strcmp(variables{i}, 'Default')~=1) 
        vars.(['v' int2str(i) '_parameter']) = variables{(i-1)*2+1}; 
        vars.(['v' int2str(i) '_vars']) = variables{(i-1)*2+2}; 
        % make shorter var name for simulation directory 
        n = textscan(vars.(['v' int2str(i) '_parameter']),'%s', 'delimiter','_'); 
        if(length(n{1})==1); x = 4; elseif(length(n{1})==2); x = 2; else; x = 1; end
        for a = 1:length(n{1})
            vars.(['v' int2str(i) '_name']) = ...
                [vars.(['v' int2str(i) '_name']) n{1}{a}(1:min(x,length(n{1}{a})))];
        end
        vars.(['v' int2str(i) '_name']) = upper(vars.(['v' int2str(i) '_name']));
        if(ischar(vars.(['v' int2str(i) '_vars']){1})~=1)
            vars.(['v' int2str(i) '_vars']) = num2cell(vars.(['v' int2str(i) '_vars']){1});
        end
    else
        % set var parameter name and variables to empty if none presented
        vars.(['v' int2str(i) '_parameter']) = 'None'; vars.(['v' int2str(i) '_vars']) = [];
    end
    if(isempty(vars.(['v' int2str(i) '_name']))~=1) 
        vars.test_f = [vars.test_f vars.(['v' int2str(i) '_name']) '_']; 
    end
end
v1_L = max(1,length(vars.v1_vars)); 
v2_L = max(1,length(vars.v2_vars));
v3_L = max(1,length(vars.v3_vars));
if(v1_L==1 && v2_L==1 && v3_L ==1); vars.test_f = [vars.test_f 'DEFAULT_']; end % if default settings
vars.test_f = [vars.test_f int2str(trials) 'T']; t = vars.test_f;
for i=1:10 % add a version ID if directory already exists
    if(exist(t,'dir')==7); t = [vars.test_f '_V' int2str(i)]; else; vars.test_f = t; break; end
end

% create default directory names when no variables given
if(isempty(vars.v3_name)==1); vars.v3_name = 'Simulations'; end
if(isempty(vars.v1_name)==1); vars.v1_name = 'Default'; end
h1 = waitbar(0, 'Parameter Progression', 'Units', 'normalized', 'Position', [0.5 0.55 0.2 0.1]);

%% SIMULATE
for k=1:v3_L % loop variable 3
    if(isempty(vars.v3_vars)~=1)
        par.([vars.v3_parameter]) = vars.v3_vars{k}; % change variable 3
        if(ischar(vars.v3_vars{k})~=1); v3 = num2str(vars.v3_vars{k}); 
        else; v3 = vars.v3_vars{k}; end
    else; v3 = []; 
    end; filename_1 = [vars.test_f '/' vars.v3_name v3];
    
    for i=1:v1_L % loop variable 1
        if(isempty(vars.v1_vars)~=1) 
            par.([vars.v1_parameter]) = vars.v1_vars{i}; % change variable 1
            if(ischar(vars.v1_vars{i})~=1); v1 = num2str(vars.v1_vars{i}); 
            else; v1 = vars.v1_vars{i}; end
        else; v1 = []; 
        end; filename_2 = [filename_1 '/' vars.v1_name v1];
        
        for j=1:v2_L % loop variable 2
            if(isempty(vars.v2_vars)~=1)
                par.([vars.v2_parameter]) = vars.v2_vars{j}; % change variable 2
                if(ischar(vars.v2_vars{j})~=1); v2 = num2str(vars.v2_vars{j}); 
                else; v2 = vars.v2_vars{j}; end
                filename_3 = [filename_2 '_' vars.v2_name v2];
            else; filename_3 = filename_2;
            end 
            
            for t=1:trials % loop trials
                % make directory
                filename = [filename_3 '/T' int2str(t)]; mkdir(filename);
                waitbar(((k-1)*vars.trials*v1_L*v2_L + (i-1)*vars.trials*v2_L + ...
                    (j-1)*vars.trials + t) / (vars.trials * v3_L * v1_L * v2_L),h1);
                % set parameters, run simulation, save current parameter set
                set_parameters(); recall_experiment(filename); vars.par{k}{i,j}{t} = par;
            end
        end
    end
end

close(h1);

%% SAVE VARIABLES
save([vars.test_f '/variables.mat'],'vars')

end


