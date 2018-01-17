function [ ] = analyse_network( filename, type )
% Analyse simulation data in the directory specified by filename. Analyse
% by type ('FR', 'WC', 'TFA' for Activity, Weight Change & Time Frequency
% Analysis, respectively).

%% INITIALISE DATA STRUCTURES
load([filename '/variables.mat']); par = vars.par;
v1_L = max(1,length(vars.v1_vars)); 
v2_L = max(1,length(vars.v2_vars));
v3_L = max(1,length(vars.v3_vars));
N = length(par{1}{1,1}{1}.sim_order_n);
n_sim = cell(1,N); n_group = {'NC', 'Hip'};
n_syn_g = {'P','NP','NP_P','P_NP'}; n_stage = {'BL','DL','AL'};
n_sub_g = {'NC_P','NC_NP','Hip_P','Hip_NP','NC_All','Hip_All'}; 
n_BS = 1000;

for n=1:N
    n_sim{n} = strrep(par{1}{1,1}{1}.sim_order_n{n},'-','_');
%     for p=1:4
%         if(p<3); FR.([n_sim{n}]).([n_subG{p}]) = cell(v3_L,1); end
%         if(p<3); FR.HIST.([n_sim{n}]).([n_subG{p}]) = cell(v3_L,1); end
%         if(p<3); FR.NC.HIST.([n_sim{n}]).([n_subG{p}]) = cell(v3_L,1); end
%         if(p<3); FR.PHASE.([n_sim{n}]).([n_subG{p}]) = cell(v3_L,1); end
%         WC.([n_sim{n}]).([n_subG{p}]) = cell(length(vars.v3_vars),1); 
%         t_w = cell(v2_L,1);
%         for k=1:v3_L % loop variable 1
%             if(p<3); FR.([n_sim{n}]).([n_subG{p}]){k} = zeros(v1_L,v2_L); end
%             if(p<3); FR.HIST.([n_sim{n}]).([n_subG{p}]){k} = cell(v1_L,v2_L); end
%             if(p<3); FR.NC.HIST.([n_sim{n}]).([n_subG{p}]){k} = cell(v1_L,v2_L); end
%             if(p<3); FR.PHASE.([n_sim{n}]).([n_subG{p}]){k} = cell(v1_L,v2_L); end
%             WC.([n_sim{n}]).([n_subG{p}]){k} = cell(v1_L,v2_L);
%             t_w{k} = cell(v1_L,v2_L);
%             for i=1:v1_L % loop variable 2
%                 for j=1:v2_L % loop variable 3
%                     if(isempty(vars.v2_vars)~=1)
%                         if(ischar(vars.v2_vars{j})~=1); v2 = int2str(vars.v2_vars{j}); 
%                         else; v2 = vars.v2_vars{j}; end
%                         filename_t = [vars.test_f '_' vars.v2_name v2];
%                     end
%                     if(p<3); FR.HIST.([n_sim{n}]).([n_subG{p}]){k}{i,j} = zeros(vars.trials,100); end
%                     %if(p<3); FR.NC.HIST.([namesF{n}]).([namesG{p}]){k}{i,j} = zeros(vars.trials,100); end
%                     if(p<3); FR.PHASE.([n_sim{n}]).([n_subG{p}]){k}{i,j} = cell(vars.trials,1); end
%                     WC.([n_sim{n}]).([n_subG{p}]){k}{i,j} = [];
%                     t_w{k}{i,j} = cell(vars.trials,1);
%                     for t=1:vars.trials % loop trials
%                         t_w{k}{i,j}{t} = ones(4,1);
%                     end
%                 end
%             end
%         end
%     end
end
n_sim{N+1} = 'All';

%% LOAD, CALCULATE & STORE SPIKE & WEIGHT DATA
% check to see if function has already been performed
if(exist([vars.test_f '/Analysis/firing-rate-data.mat'],'file')~=2 || exist([vars.test_f '/Analysis/weight-change-data.mat'],'file')~=2) 
    h1 = waitbar(0, 'Parameter Progression', 'Units', 'normalized', 'Position', [0.5 0.55 0.2 0.1]);
    for t=1:vars.trials % loop trials
       % LOOP THROUGH FILENAMES
       for k=1:v3_L % loop variable 1
           if(isempty(vars.v3_vars)~=1) 
                if(ischar(vars.v3_vars{k})~=1); v3 = num2str(vars.v3_vars{k}); 
                else; v3 = vars.v3_vars{k}; end
           else; v3 = []; 
            end; filename_1 = [vars.test_f '/' vars.v3_name v3];
            for i=1:v1_L % loop variable 2
                if(isempty(vars.v1_vars)~=1) 
                    if(ischar(vars.v1_vars{i})~=1); v1 = num2str(vars.v1_vars{i}); 
                    else; v1 = vars.v1_vars{i}; end
                else; v1 = []; 
                end; filename_2 = [filename_1 '/' vars.v1_name v1];
                for j=1:v2_L % loop variable 3
                    if(isempty(vars.v2_vars)~=1)
                        if(ischar(vars.v2_vars{j})~=1); v2 = num2str(vars.v2_vars{j}); 
                        else; v2 = vars.v2_vars{j}; end
                        filename_3 = [filename_2 '_' vars.v2_name v2];
                    else; filename_3 = filename_2;
                    end
                    filename = [filename_3 '/T' int2str(t)];
                    waitbar(((t-1)*v3_L*v1_L*v2_L + (k-1)*v1_L*v2_L + ...
                        (i-1)*v2_L + j) / (vars.trials * v3_L * v1_L * v2_L),h1); % increment progress bar

                    NCPI = par{k}{i,j}{t}.NC_per_item; HPI = par{k}{i,j}{t}.Hip_per_item; % extract network properties
                    for n=1:N % loop through stages of simulation
                         load([filename '/Data/' par{k}{i,j}{t}.sim_order_n{n} '.mat']) % load data
                         % extract spikes/weight/synapse data
                         spikes = data.sim_stats.spike_detector; weights_t = data.sim_stats.weight_matrix; 
                         syn = data.sim_stats.weight_matrix_STDP; 
                         sim_length = unique(data.sim_stats.voltmeter(:,2));
                         
                         %% EXTRACT INPUT RECORDING DATA
                         rec = fieldnames(data.sim_stats.I_REC);
                         for r = 1:length(rec)
                             if(t==1); FR.I.([n_sim{n}]).([rec{r}]){k}{i,j} = zeros(1, length(sim_length)); end % initialise data
                             FR.I.([n_sim{n}]).([rec{r}]){k}{i,j}  = FR.I.([n_sim{n}]).([rec{r}]){k}{i,j}  + data.sim_stats.I_REC.([rec{r}]); % increment data
                             if(t==vars.trials); FR.I.([n_sim{n}]).([rec{r}]){k}{i,j}  = (FR.I.([n_sim{n}]).([rec{r}]){k}{i,j}  / vars.trials) / par{k}{i,j}{t}.C_m; end% average data
                         end
                         clear('rec');
                         
                         %% EXTRACT ALL HIPPOCAMPAL RELATED DATA
                         for p = 1:2 
                             %% EXTRACT WEIGHT CHANGE OF HIP GROUPS
                             for w=1:2 % loop through weight groups
                                 if(w==p); v=0; else; v=1; end % extract all weight changes in group
                                 S = syn(NCPI*2 + HPI*v+1:NCPI*2 + HPI*(v+1) , NCPI*2 + HPI*(p-1)+1:NCPI*2 + HPI*p);
                                 wc = weights_t(NCPI*2 + HPI*v+1:NCPI*2 + HPI*(v+1) , NCPI*2 + HPI*(p-1)+1:NCPI*2 + HPI*p,:);
                                 [wy, wx] = find(S==1); % find non-zeros/active synapses in simulation
                                 W = zeros(length(sim_length),length(wx));
                                 for l=1:length(wx)
                                     W(:,l) = squeeze(wc(wy(l),wx(l),:));
                                 end % add to temporary matrix 
                                 if(t == 1)
                                     WC.([n_sim{n}]).([n_syn_g{p+(w-1)*2}]){k}{i,j} = [];
                                 end
                                 WC.([n_sim{n}]).([n_syn_g{p+(w-1)*2}]){k}{i,j} = ...
                                     [WC.([n_sim{n}]).([n_syn_g{p+(w-1)*2}]){k}{i,j} mean(W,2)]; % concatenate matrix
                                 if(t == vars.trials)
                                     WC.([n_sim{n}]).([n_syn_g{p+(w-1)*2}]){k}{i,j} = mean(WC.([n_sim{n}]).([n_syn_g{p+(w-1)*2}]){k}{i,j},2);
                                 end
                             end
                             
                             %% EXTRACT Y-DATA FOR HIPPOCAMPAL ACTIVITY PLOTS
                             % extract P/NP spikes
                             FR.ACT.bin_width = 20;
                             hip_spikes = spikes(spikes(:,1) > NCPI*2 + HPI *(p-1),:);
                             hip_spikes = hip_spikes(hip_spikes(:,1) <= HPI*p + NCPI*2,:);
                             if(t==1);FR.ACT.([n_sim{n}]).([n_syn_g{p}]){k}{i,j}.S = zeros(vars.trials,round(length(sim_length)/FR.ACT.bin_width));end
                             FR.ACT.([n_sim{n}]).([n_syn_g{p}]){k}{i,j}.S(t, :) = histcounts(hip_spikes(:,2),length(sim_length)/FR.ACT.bin_width,'binlimits',[0 length(sim_length)]) ...
                                 / (par{k}{i,j}{t}.n_Hip/2) * (1000 / FR.ACT.bin_width);
                             if(t==vars.trials)
                                 if(sum(cellfun(@sum,strfind(type,'FR')))>0)
                                 [FR.ACT.([n_sim{n}]).([n_syn_g{p}]){k}{i,j}.CI_L, FR.ACT.([n_sim{n}]).([n_syn_g{p}]){k}{i,j}.CI_U, ...
                                     FR.ACT.([n_sim{n}]).([n_syn_g{p}]){k}{i,j}.M] = ...
                                     bootstrap(FR.ACT.([n_sim{n}]).([n_syn_g{p}]){k}{i,j}.S, n_BS, []);
                                 end
                                 FR.ACT.([n_sim{n}]).([n_syn_g{p}]){k}{i,j} = ...
                                     rmfield(FR.ACT.([n_sim{n}]).([n_syn_g{p}]){k}{i,j},'S'); % remove field for data size
                             end
                             
                             %% EXTRACT AND CALCULATE DEGREES/RADIANS PHASE DATA FOR HIPPOCAMPAL NEURONS
                             hip_spikes(:,2) = hip_spikes(:,2) + par{k}{i,j}{t}.Hip_r_phase_n(n)*1000; % adjust for random Theta phase
                             phase_deg = [hip_spikes(:,1) (hip_spikes(:,2)/250 - floor(hip_spikes(:,2)/250))*360]; % degrees
                             phase_rad = [hip_spikes(:,1) phase_deg(:,2)*(pi/180)]; % radians
                             FR.PHASE.([n_sim{n}]).([n_syn_g{p}]){k}{i,j}{t} = [phase_deg phase_rad(:,2)]; % store data
                         end
                    end
                    %% EXTRACT AFTER-LEARNING (AL) WEIGHTS FOR HIPPOCAMPAL GROUPS
                    for p=1:4
                        WC.AL_W.([n_syn_g{p}]){k}{i,j}(t,1) = mean(WC.DL.([n_syn_g{p}]){k}{i,j}(end));
                        WC.REC_W.([n_syn_g{p}]){k}{i,j}(t,1) = mean([WC.P_AL.([n_syn_g{p}]){k}{i,j}(end) WC.NP_AL.([n_syn_g{p}]){k}{i,j}(end)]);
                    end
                end
            end
       end
    end
    close(h1)
    % SAVE ANALYSIS DATA
    mkdir([vars.test_f '/Analysis'])
    save([vars.test_f '/Analysis/firing-rate-data.mat'],'FR','-v7.3') 
    save([vars.test_f '/Analysis/weight-change-data.mat'],'WC','-v7.3')
end

%% AVERAGE FIRING RATE PLOTS OVER TRIALS
if(sum(cellfun(@sum,strfind(type,'FR')))>0)
    load([vars.test_f '/Analysis/firing-rate-data.mat']); % load activity data
    % loop through filenames
    for k=1:v3_L % loop variable 1
        if(isempty(vars.v3_vars)~=1)
            if(ischar(vars.v3_vars{k})~=1); v3 = num2str(vars.v3_vars{k}); 
            else; v3 = vars.v3_vars{k}; end
        else; v3 = []; 
        end; filename_1 = [vars.test_f '/' vars.v3_name v3];
        for i=1:v1_L % loop variable 2
            if(isempty(vars.v1_vars)~=1)
                if(ischar(vars.v1_vars{i})~=1); v1 = num2str(vars.v1_vars{i}); 
                else; v1 = vars.v1_vars{i}; end
            else; v1 = []; 
            end; filename_2 = [filename_1 '/' vars.v1_name v1];
            for j=1:v2_L % loop variable 3
                if(isempty(vars.v2_vars)~=1)
                    if(ischar(vars.v2_vars{j})~=1); v2 = num2str(vars.v2_vars{j}); 
                    else; v2 = vars.v2_vars{j}; end
                    filename_3 = [filename_2 '_' vars.v2_name v2];
                else; filename_3 = filename_2;
                end
                filename = filename_3;
                mkdir([filename '/Analysis/Firing Rate']);
                
                %% EXTRACT & SMOOTH DATA
                x1 = -par{k}{i,j}{t}.pre_stim_length + 1: 1 : par{k}{i,j}{t}.stim_length;
                sim_length = par{k}{i,j}{t}.pre_stim_length + par{k}{i,j}{t}.stim_length; bx = 10;
                I = FR.I.DL.BG_HIP{k}{i,j} + FR.I.DL.NC_HIP{k}{i,j} ...
                    + FR.I.DL.HIP_HIP_WIT{k}{i,j} + FR.I.DL.HIP_HIP_BET{k}{i,j} ...
                    + FR.I.DL.ADP_HIP{k}{i,j};
                I_F = zeros(1, sim_length);
                I_NC_F = zeros(1, sim_length); I_BG_F = zeros(1, sim_length); I_ADP_F = zeros(1, sim_length);
                I_Hip_W_F = zeros(1, sim_length); I_Hip_B_F = zeros(1, sim_length);
                I_Hip_BL = zeros(1, sim_length); I_Hip_AL = zeros(1, sim_length);
                for l = 1:sim_length % smooth data with box 
                    I_F(l) = mean(I(max(1,l-bx):min(sim_length,l+bx)));
                    I_BG_F(l) = mean(FR.I.DL.BG_HIP{k}{i,j}(max(1,l-bx):min(sim_length,l+bx)));
                    I_NC_F(l) = mean(FR.I.DL.NC_HIP{k}{i,j}(max(1,l-bx):min(sim_length,l+bx)));
                    I_ADP_F(l) = mean(FR.I.DL.ADP_HIP{k}{i,j}(max(1,l-bx):min(sim_length,l+bx)));
                    I_Hip_W_F(l) = mean(FR.I.DL.HIP_HIP_WIT{k}{i,j}(max(1,l-bx):min(sim_length,l+bx)));
                    I_Hip_B_F(l) = mean(FR.I.DL.HIP_HIP_BET{k}{i,j}(max(1,l-bx):min(sim_length,l+bx)));
                    I_Hip_BL(l) = mean(FR.I.NP_BL.HIP_HIP_BET{k}{i,j}(max(1,l-bx):min(sim_length,l+bx))) + ...
                        mean(FR.I.P_BL.HIP_HIP_BET{k}{i,j}(max(1,l-bx):min(sim_length,l+bx)));
                    I_Hip_AL(l) = mean(FR.I.NP_AL.HIP_HIP_BET{k}{i,j}(max(1,l-bx):min(sim_length,l+bx))) + ...
                        mean(FR.I.P_AL.HIP_HIP_BET{k}{i,j}(max(1,l-bx):min(sim_length,l+bx)));
                end
                
                x = -par{k}{i,j}{t}.pre_stim_length : FR.ACT.bin_width : par{k}{i,j}{t}.stim_length-FR.ACT.bin_width;
                x = x / 1000; x1 = x1/1000;
                recognition_BL = [(FR.ACT.NP_BL.NP{k}{i,j}.CI_L + FR.ACT.P_BL.P{k}{i,j}.CI_L) / 2; ...
                    (FR.ACT.NP_BL.NP{k}{i,j}.M + FR.ACT.P_BL.P{k}{i,j}.M) / 2; ...
                    (FR.ACT.NP_BL.NP{k}{i,j}.CI_U + FR.ACT.P_BL.P{k}{i,j}.CI_U) / 2];
                recognition_AL = [(FR.ACT.NP_AL.NP{k}{i,j}.CI_L + FR.ACT.P_AL.P{k}{i,j}.CI_L) / 2; ...
                    (FR.ACT.NP_AL.NP{k}{i,j}.M + FR.ACT.P_AL.P{k}{i,j}.M) / 2; ...
                    (FR.ACT.NP_AL.NP{k}{i,j}.CI_U + FR.ACT.P_AL.P{k}{i,j}.CI_U) / 2];
                recall_BL = [(FR.ACT.NP_BL.P{k}{i,j}.CI_L + FR.ACT.P_BL.NP{k}{i,j}.CI_L) / 2; ...
                    (FR.ACT.NP_BL.P{k}{i,j}.M + FR.ACT.P_BL.NP{k}{i,j}.M) / 2; ...
                    (FR.ACT.NP_BL.P{k}{i,j}.CI_U + FR.ACT.P_BL.NP{k}{i,j}.CI_U) / 2];
                recall_AL = [(FR.ACT.NP_AL.P{k}{i,j}.CI_L + FR.ACT.P_AL.NP{k}{i,j}.CI_L) / 2; ...
                    (FR.ACT.NP_AL.P{k}{i,j}.M + FR.ACT.P_AL.NP{k}{i,j}.M) / 2; ...
                    (FR.ACT.NP_AL.P{k}{i,j}.CI_U + FR.ACT.P_AL.NP{k}{i,j}.CI_U) / 2];
                DL = [(FR.ACT.DL.P{k}{i,j}.CI_L + FR.ACT.DL.NP{k}{i,j}.CI_L) / 2; ...
                    (FR.ACT.DL.P{k}{i,j}.M + FR.ACT.DL.NP{k}{i,j}.M) / 2; ...
                    (FR.ACT.DL.P{k}{i,j}.CI_U + FR.ACT.DL.NP{k}{i,j}.CI_U) / 2];
                y_max = ceil(max(max([recognition_BL; recognition_AL; recall_BL; recall_AL; DL]))/10)*10; a = 0.2;
                
                %% PLOT recognition AND RECALL ACTIVITY
                f = figure; set(f, 'Position', [0 0 1080 800]); % NP presentations
                subplot(5,6,[1 2 7 8]); hold on; 
                text(0.02,0.98,'Ai','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',16,'Fontweight','bold')
                fill([x fliplr(x)],[recognition_BL(1,:) fliplr(recognition_BL(3,:))], 'b','edgecolor','b'); alpha(a); % plot mean distribution
                plot(x,recognition_BL(2,:),'color', 'b', 'linestyle', '-','linewidth',1); % plot average
                ylim([0 y_max]); xlim([max(-0.5,min(x)) min(1.5,max(x))]); title('Recognition BL');
                ylabel('Activity (Hz)'); ytickangle(90); xticklabels([]); ax = gca; ax.FontSize = 14;

                subplot(5,6,[13 14 19 20]); hold on; 
                text(0.02,0.98,'Aii','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',16,'Fontweight','bold')
                fill([x fliplr(x)],[recall_BL(1,:) fliplr(recall_BL(3,:))], [0.8 0.8 0.8],'edgecolor','k'); alpha(a); % plot mean distribution
                plot(x,recall_BL(2,:),'color', 'k', 'linestyle', '-','linewidth',1); % plot average
                xlim([max(-0.5,min(x)) min(1.5,max(x))]); ylim([0 y_max]);title('Recall BL'); 
                ylabel('Activity (Hz)'); ytickangle(90); xticklabels([]); ax = gca; ax.FontSize = 14;

                subplot(5,6,[3 4 9 10]); hold on; 
                text(0.02,0.98,'Bi','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',16,'Fontweight','bold')
                fill([x fliplr(x)],[recognition_AL(1,:) fliplr(recognition_AL(3,:))], 'r','edgecolor','r'); alpha(a); % plot mean distribution
                plot(x,recognition_AL(2,:),'color', 'r', 'linestyle', '-','linewidth',1); % plot average
                title('Recognition AL'); xlim([max(-0.5,min(x)) min(1.5,max(x))]); ylim([0 y_max]); 
                xticklabels([]); yticklabels([]); ax = gca; ax.FontSize = 14;

                subplot(5,6,[15 16 21 22]); hold on; 
                text(0.02,0.98,'Bii','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',16,'Fontweight','bold')
                fill([x fliplr(x)],[recall_AL(1,:) fliplr(recall_AL(3,:))], [0.2 1 0.2],'edgecolor',[0 0.6 0]); alpha(a); % plot mean distribution
                plot(x,recall_AL(2,:),'color', [0 0.6 0], 'linestyle', '-','linewidth',1); % plot average
                xlim([max(-0.5,min(x)) min(1.5,max(x))]); ylim([0 y_max]);title('Recall AL'); 
                xticklabels([]); yticklabels([]); xticklabels([]); ax = gca; ax.FontSize = 14;
                
                %% PLOT DURING LEARNING
                % ACTIVITY
                subplot(5,6,[5 6 11 12]); hold on;
                text(0.02,0.98,'Ci','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',16,'Fontweight','bold')
                fill([x fliplr(x)],[DL(1,:) fliplr(DL(3,:))], [0 0.4 0.6],'edgecolor',[0 0.4 0.6]); alpha(a); % plot mean distribution
                plot(x,DL(2,:),'color', [0 0.4 0.6], 'linestyle', '-','linewidth',1); % plot average
                xlim([max(-0.5,min(x)) min(1.5,max(x))]); ylim([0 y_max]); title('DL'); 
                xticklabels([]); yticklabels([]); ax = gca; ax.FontSize = 14;
                
                % INPUT BREAKDOWN
                subplot(5,6,[17 18 23 24]); hold on;
                text(0.02,0.98,'Cii','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',16,'Fontweight','bold')
                area(x1, I, 'facecolor', [0.1 0.1 0.3],'linestyle','none'); alpha(0.1);
                plot(x1, I_NC_F + I_BG_F, 'color', [0.8 0.3 0.3], 'linestyle', '-', 'marker', 'none','linewidth',2);
                plot(x1, I_ADP_F, 'color', [1 .8 0], 'linestyle', '-', 'marker', 'none','linewidth',2);
                plot(x1, I_Hip_W_F , 'color', [0.3 0.3 0.8], 'linestyle', '-', 'marker', 'none','linewidth',2);
                plot(x1, I_Hip_B_F , 'color', [0.3 0.3 0.8], 'linestyle', ':', 'marker', 'none','linewidth',2);
                title('Input Breakdown DL');  ylabel('I / C_m'); ytickangle(90); xticklabels([]);
                ylim([0 max(I_F)*1.1]); xlim([max(-0.5,min(x)) min(1.5,max(x))]); ax = gca; ax.FontSize = 14;
                legend('I', 'I_e_x_t', 'I_A_D_P', 'I_H', 'I_H_<_>_H','location','northeast');
                
                %% PLOT RASTERS
                r_trial = ceil(rand()*vars.trials); % choose random trial
                r_P =  par{k}{i,j}{t}.n_NC + ceil(rand() * par{k}{i,j}{t}.n_Hip/2); % choose random P neuron
                r_NP = par{k}{i,j}{t}.n_NC + ceil(rand() * par{k}{i,j}{t}.n_Hip/2) + par{k}{i,j}{t}.n_Hip/2; % choose random NP neuron
                % load BL data for trial
                load([vars.test_f '/Simulations/Default/T' int2str(r_trial) '/Data/P-BL.mat']); 
                spikes = data.sim_stats.spike_detector; spikes(:,2) = spikes(:,2) - par{k}{i,j}{t}.pre_stim_length;
                P_BL = spikes(spikes(:,1)==r_P,:); NP_BL = spikes(spikes(:,1)==r_NP,:); clear('data');
                P_BL(:,1) = P_BL(:,1) - r_P + 2; NP_BL(:,1) = NP_BL(:,1) - r_NP + 1;
                % load AL data for trial
                load([vars.test_f '/Simulations/Default/T' int2str(r_trial) '/Data/P-AL.mat']); 
                spikes = data.sim_stats.spike_detector; spikes(:,2) = spikes(:,2) - par{k}{i,j}{t}.pre_stim_length;
                P_AL = spikes(spikes(:,1)==r_P,:); NP_AL = spikes(spikes(:,1)==r_NP,:); clear('data');
                P_AL(:,1) = P_AL(:,1) - r_P + 2; NP_AL(:,1) = NP_AL(:,1) - r_NP + 1;
                % load DL data for trial
                load([vars.test_f '/Simulations/Default/T' int2str(r_trial) '/Data/DL.mat']); 
                spikes = data.sim_stats.spike_detector; spikes(:,2) = spikes(:,2) - par{k}{i,j}{t}.pre_stim_length;
                P_DL = spikes(spikes(:,1)==r_P,:); NP_DL = spikes(spikes(:,1)==r_NP,:); clear('data');
                P_DL(:,1) = P_DL(:,1) - r_P + 2; NP_DL(:,1) = NP_DL(:,1) - r_NP + 1;
                
                subplot(5,6,25:26); hold on; ylim([0.5 5.5]); % BL
                text(0.02,0.98,'Aiii','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',16,'Fontweight','bold')
                for s = 1:length(P_BL); plot([P_BL(s,2) P_BL(s,2)], [3 4],'b'); end
                for s = 1:length(NP_BL); plot([NP_BL(s,2) NP_BL(s,2)], [1 2],'k'); end
                yticks([1.5 3.5]); yticklabels({'NP','P'}); xlim([-500 1500]); 
                xticks(-500:500:1500); xticklabels({'-0.5','0','0.5','1','1.5'});
                xlabel('Time (s)'); title('P BL'); ax = gca; ax.FontSize = 14;
                
                subplot(5,6,27:28); hold on; ylim([0.5 5.5]); % AL
                text(0.02,0.98,'Biii','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',16,'Fontweight','bold')
                for s = 1:length(P_AL); plot([P_AL(s,2) P_AL(s,2)], [3 4],'r'); end
                for s = 1:length(NP_AL); plot([NP_AL(s,2) NP_AL(s,2)], [1 2],'color', [0 0.6 0]); end
                yticks([1.5 3.5]); yticklabels([]); xlim([-500 1500]); 
                xticks(-500:500:1500); xticklabels({'-0.5','0','0.5','1','1.5'});
                xlabel('Time (s)'); title('P AL'); ax = gca; ax.FontSize = 14;
                
                subplot(5,6,29:30); hold on; ylim([0.5 5.5]); % DL
                text(0.02,0.98,'Ciii','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',16,'Fontweight','bold')
                for s = 1:length(P_DL); plot([P_DL(s,2) P_DL(s,2)], [3 4],'color',[0 0.4 0.6]); end
                for s = 1:length(NP_DL); plot([NP_DL(s,2) NP_DL(s,2)], [1 2],'color',[0 0.4 0.6]); end
                yticks([1.5 3.5]); yticklabels([]); xlim([-500 1500]); 
                xticks(-500:500:1500); xticklabels({'-0.5','0','0.5','1','1.5'});
                xlabel('Time (s)'); title('DL'); ax = gca; ax.FontSize = 14;
                saveas(f,[filename '/Analysis/Firing Rate/Activity_all.tiff']); close(f);
                
                %% MATCH ISON et al. FIGURE 5
                f = figure; set(f, 'Position', [0 0 680 800]); % NP presentations
                sm_par = 1e-7; a = 0.2; x = x*1000;
                resp_BL_L = fit(transpose(x), transpose(recognition_BL(1,:)),'smoothingspline','smoothingparam',sm_par);
                resp_BL_M = fit(transpose(x), transpose(recognition_BL(2,:)),'smoothingspline','smoothingparam',sm_par);
                resp_BL_U = fit(transpose(x), transpose(recognition_BL(3,:)),'smoothingspline','smoothingparam',sm_par);
                
                resp_AL_L = fit(transpose(x), transpose(recognition_AL(1,:)),'smoothingspline','smoothingparam',sm_par);
                resp_AL_M = fit(transpose(x), transpose(recognition_AL(2,:)),'smoothingspline','smoothingparam',sm_par);
                resp_AL_U = fit(transpose(x), transpose(recognition_AL(3,:)),'smoothingspline','smoothingparam',sm_par);
                
                rec_BL_L = fit(transpose(x), transpose(recall_BL(1,:)),'smoothingspline','smoothingparam',sm_par);
                rec_BL_M = fit(transpose(x), transpose(recall_BL(2,:)),'smoothingspline','smoothingparam',sm_par);
                rec_BL_U = fit(transpose(x), transpose(recall_BL(3,:)),'smoothingspline','smoothingparam',sm_par);
                
                rec_AL_L = fit(transpose(x), transpose(recall_AL(1,:)),'smoothingspline','smoothingparam',sm_par);
                rec_AL_M = fit(transpose(x), transpose(recall_AL(2,:)),'smoothingspline','smoothingparam',sm_par);
                rec_AL_U = fit(transpose(x), transpose(recall_AL(3,:)),'smoothingspline','smoothingparam',sm_par);
                
                %f = figure; set(f, 'Position', [0 0 800 500]); % NP presentations
                subplot(5,4,[1 2 5 6]);  hold on; 
                text(0.02,0.98,'Di','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14,'Fontweight','bold')
                h(1) = fill([x fliplr(x)],[resp_BL_L(x); flipud(resp_BL_U(x))], 'b','edgecolor','b'); alpha(a); % plot mean distribution
                plot(x, resp_BL_M(x),'color', 'b', 'linestyle', '-','linewidth',2); % plot average
                h(2) = fill([x fliplr(x)],[resp_AL_L(x); flipud(resp_AL_U(x))], 'r','edgecolor','r'); alpha(a); % plot mean distribution
                plot(x, resp_AL_M(x),'color', 'r', 'linestyle', '-','linewidth',2); % plot average
                title('Recognition');  legend([h(1), h(2)], 'Before Learning', 'After Learning','location','northeast');
                xlim([max(-750,min(x)) min(1000,max(x))]); ylim([0 y_max]);  ylabel('Activity (Hz)'); ytickangle(90);
                
                subplot(5,4,[3 4 7 8]); hold on;
                text(0.02,0.98,'Dii','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14,'Fontweight','bold')
                h(1) = fill([x fliplr(x)],[rec_BL_L(x); flipud(rec_BL_U(x))], [0.8 0.8 0.8],'edgecolor','k'); alpha(a); % plot mean distribution
                plot(x, rec_BL_M(x),'color', 'k', 'linestyle', '-','linewidth',2); % plot average
                h(2) = fill([x fliplr(x)],[rec_AL_L(x); flipud(rec_AL_U(x))], [0.2 1 0.2],'edgecolor',[0 0.6 0]); alpha(a); % plot mean distribution
                plot(x, rec_AL_M(x),'color', [0 0.6 0], 'linestyle', '-','linewidth',2); % plot average
                title('Recall'); legend([h(1), h(2)], 'Before Learning', 'After Learning','location','northeast');
                xlim([max(-750,min(x)) min(1000,max(x))]); ylim([0 y_max]);  yticklabels([]);
                
                % ISON REAL DATA
                ison = imread('ison.jpg');
                ax = subplot(5,4,9:20); hold on
                text(0.02,0.98,'Diii','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14,'Fontweight','bold')
                RI = imref2d(size(ison)); RI.XWorldLimits = [0 RI.ImageSize(2)]; RI.YWorldLimits = [0 RI.ImageSize(1)]; 
                xlim([0 RI.ImageSize(2)]); ylim([0 RI.ImageSize(1)]);
                imshow(ison,RI, 'Border','tight'); xlabel({'Rapid Encoding of New Memories by Individual Neurons in the Human Brain'; 'Ison et al. 2015'});
                set(ax, 'box', 'on', 'Visible', 'on', 'xtick', [], 'ytick', []);
                saveas(f,[filename '/Analysis/Firing Rate/Activity_Ison.tiff']); close(f);
                
                %% EXTRACT DATA FOR POLAR HISTOGRAM GRAPHS
                for n=1:N
                    for p=1:2
                        %mkdir([filename '/Analysis/Phase Precession/' namesG{p}]);
                        phase_rad = []; fs = zeros(vars.trials,par{k}{i,j}{t}.Hip_per_item);
                        for t=1:vars.trials
                            phase_rad = [phase_rad; FR.PHASE.([n_sim{n}]).([n_syn_g{p}]){k}{i,j}{t}(:,3)];
                            x = unique(FR.PHASE.([n_sim{n}]).([n_syn_g{p}]){k}{i,j}{t}(:,1));
                            x_ID = x - par{k}{i,j}{t}.n_NC - par{k}{i,j}{t}.Hip_per_item*(p-1);
                            rad = [FR.PHASE.([n_sim{n}]).([n_syn_g{p}]){k}{i,j}{t}(:,1)...
                                FR.PHASE.([n_sim{n}]).([n_syn_g{p}]){k}{i,j}{t}(:,3)];
                             for l=1:length(x) % extract first spikes
                                 m = rad(rad(:,1)==x(l),2);
                                 m2 = min(m(m>180*(pi/180))); % take first spike after 180deg
                                 if(isempty(m2)==1); m2 = min(m); end % or take first spike after 0deg
                                 fs(t,x_ID(l)) = m2;
                             end
                        end
                        first_spike.([n_sim{n}]).([n_syn_g{p}]) = fs; 
                        phase_spikes.([n_sim{n}]).([n_syn_g{p}]) = phase_rad;
                    end
                end
                %% PLOT POLAR HISTOGRAM FOR THETA PHASE
                titles = {'Recall','Recognition'}; binwidth = pi/2; binnorm = 'probability'; a = 0.8; sp = {'A','B','C'};
                for p = 1:2
                    f(6) = figure(6); set(f(6), 'Position', [0 0 900 600]); hold on;
                    annotation('textbox', [0 0.9 1 0.1], 'String', [titles{p} ' All Spikes'], 'EdgeColor', 'none', 'FontSize', 16,'HorizontalAlignment', 'center')
                    annotation('textbox', [0 0.425 1 0.1], 'String', [titles{p} ' First Spikes'], 'EdgeColor', 'none', 'FontSize', 16,'HorizontalAlignment', 'center')
                    if(p==1); cN = {[0.3 0.3 0.3],[0 0.4 0.6],[0.2 1 0.2]}; else; cN = {'b',[0 0.4 0.6],'r'}; end
                    for i2 = 1:3
                        if(i2==1 || i2==3) % RECALL OR recognition DATA
                            FS = [first_spike.([n_syn_g{3-p} '_' n_stage{i2}]).P first_spike.([n_syn_g{p} '_' n_stage{i2}]).NP];
                            AS = [phase_spikes.([n_syn_g{3-p} '_' n_stage{i2}]).P; phase_spikes.([n_syn_g{p} '_' n_stage{i2}]).NP];
                        else % DURING LEARNING DATA
                            FS = [first_spike.DL.P first_spike.DL.NP]; AS = [phase_spikes.DL.P; phase_spikes.DL.NP];
                        end
                        FS = reshape(FS,[vars.trials*par{k}{i,j}{t}.Hip_per_item*2, 1]);

                        for i3=1:2 % LTP (green) & LTD (red) phases of Theta
                            subplot(2,3,i2+(i3-1)*3);
                            polarhistogram(2:4.5,'binwidth',pi*2,'Normalization',binnorm,'FaceColor', 'g','FaceAlpha',.1,'edgealpha',0); hold on;
                            polarhistogram(4.8:7.9,'binwidth',pi*2,'Normalization',binnorm,'FaceColor', 'r','FaceAlpha',.1,'edgealpha',0); hold off;
                        end
                        
                        subplot(2,3,i2); hold on; % plot all spikes polar hist
                        annotation('textbox', [0.115+0.275*(i2-1) 0.85 1 0.1], 'String', [sp{i2} 'i'], 'EdgeColor', 'none', 'FontSize', 18,'fontweight','bold')
                        annotation('textbox', [0.21+0.28*(i2-1) 0.475 1 0.1], 'String', n_stage{i2}, 'EdgeColor', 'none', 'FontSize', 14)
                        polarhistogram(AS,'binwidth',binwidth,'Normalization',binnorm,'FaceColor', cN{i2},'FaceAlpha',a,'edgecolor','none');
                        polarhistogram(AS,'binwidth',binwidth,'Normalization',binnorm,'displaystyle','stairs','edgecolor','k','linewidth',1,'linestyle','-');
                        rticks(0.25:0.25:1); rticklabels({'','0.5','','1'}); thetaticks([0 90 180 270]); thetaticklabels({'\pi/2','\pi','-\pi/2','0'});
                        pax = gca; pax.ThetaZeroLocation = 'top';
                        
                        subplot(2,3,i2+3); hold on; % plot first spikes polar hist
                        annotation('textbox', [0.115+0.275*(i2-1) 0.38 1 0.1], 'String', [sp{i2} 'ii'], 'EdgeColor', 'none', 'FontSize', 18,'fontweight','bold')
                        annotation('textbox', [0.21+0.28*(i2-1) 0 1 0.1], 'String', n_stage{i2}, 'EdgeColor', 'none', 'FontSize', 14)
                        polarhistogram(FS,'binwidth',binwidth,'Normalization',binnorm,'FaceColor', cN{i2},'FaceAlpha',a,'edgecolor','none');
                        polarhistogram(FS,'binwidth',binwidth,'Normalization',binnorm,'displaystyle','stairs','edgecolor','k','linewidth',1,'linestyle','-');
                        rticks(0.25:0.25:1); rticklabels({'','0.5','','1'}); thetaticks([0 90 180 270]); thetaticklabels({'\pi/2','\pi','-\pi/2','0'});
                        pax = gca; pax.ThetaZeroLocation = 'top';
                    end
                    saveas(f(6),[filename '/Analysis/Firing Rate/Phase Precession ' titles{p} '.tiff']); close(f(6));
                end
            end
        end
    end
end

%% AVERAGE WEIGHT CHANGE PLOTS
if(sum(cellfun(@sum,strfind(type,'WC')))>0)
    load([vars.test_f '/Analysis/weight-change-data.mat']);
    titles = {'Preferred Synapses','Non-Preferred Synapses','NP -> P Synapses', 'P -> NP Synapses'}; 
    f_names = {'Within_Groups','Between_Groups'}; % weight groups
    plots_r = {'-','--'}; plots_c = {'b','m'};
    % loop through filenames
    for k=1:v3_L % loop variable 1
        if(isempty(vars.v3_vars)~=1)
            if(ischar(vars.v3_vars{k})~=1); v3 = num2str(vars.v3_vars{k}); 
            else; v3 = vars.v3_vars{k}; end
        else; v3 = []; 
        end; filename_1 = [vars.test_f '/' vars.v3_name v3];
        for i=1:v1_L % loop variable 2
            if(isempty(vars.v1_vars)~=1)
                if(ischar(vars.v1_vars{i})~=1); v1 = num2str(vars.v1_vars{i}); 
                else; v1 = vars.v1_vars{i}; end
            else; v1 = []; 
            end; filename_2 = [filename_1 '/' vars.v1_name v1];
            for j=1:v2_L % loop variable 3
                if(isempty(vars.v2_vars)~=1)
                    if(ischar(vars.v2_vars{j})~=1); v2 = num2str(vars.v2_vars{j}); 
                    else; v2 = vars.v2_vars{j}; end
                    filename_3 = [filename_2 '_' vars.v2_name v2];
                else; filename_3 = filename_2;
                end
                filename = filename_3;
                mkdir([filename '/Analysis/Weight Change']);
                
                %% DEFINE BOX CO-ORDINATES FOR STIMULUS PRESENTATIONS
                Y = [-100,-100,par{k}{i,j}{t}.weight_max*1.2,par{k}{i,j}{t}.weight_max*1.2]; X = 0;
                PX = cell(1,N); c = cell(1,N);
                for n=2:N
                   if(strcmp(par{k}{i,j}{t}.sim_order{n},'idling')==1)
                       X = par{k}{i,j}{t}.idling_length/1000; 
                   else
                       pre = par{k}{i,j}{t}.pre_stim_length; 
                       if(strcmp(par{k}{i,j}{t}.sim_order{n},'P presentation')==1)
                           stim = par{k}{i,j}{t}.stim_length; c{n} = 'b'; 
                       elseif(strcmp(par{k}{i,j}{t}.sim_order{n},'NP presentation')==1)
                           stim = par{k}{i,j}{t}.stim_length; c{n} = 'm';
                       elseif(strcmp(par{k}{i,j}{t}.sim_order{n},'learning')==1)
                           stim = par{k}{i,j}{t}.stim_length; c{n} = 'g';
                       end
                       pre = pre/1000; stim=stim/1000; 
                       PX{n} = [pre, pre + stim, pre + stim, pre] + X; X = max(PX{n}); % shift box along to next stimulus
                   end
                end
                f(1) = figure(1); set(f(1), 'Position', [0 0 900 600]); % weight average figure
                %% EXTRACT AND PLOT HIPPOCAMPAL WEIGHT CHANGE
                sp = {'A', 'B'};
                for p=1:2 % loop hippocampal groups
                    subplot(2,1,p); hold on;
                    text(0.02,0.98,sp{p},'Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',22,'Fontweight','bold')
                    for w=1:2 % loop weight groups
                        average_y = [];
                        for n=2:N
                           % concatenate weight change data over simulations
                           average_y_t = transpose(WC.([n_sim{n}]).([n_syn_g{w+(p-1)*2}]){k}{i,j});
                           average_y = [average_y average_y_t];
                        end
                        x = (1:length(average_y))/1000;
                        p_h(w) = plot(x,average_y,plots_c{sqrt(((p-1)*3-w)^2)},'linewidth',2,'linestyle',plots_r{sqrt(((p-1)*3-w)^2)}); % plot average
                    end
                    for n = 2:N-1
                        if(isempty(PX{n})~=1); patch(PX{n},Y,c{n},'edgecolor','none'); end 
                        annotation('rectangle', [0.205 + 0.14765*(n-2) 0.875 - (p-1)*0.475 .071 .05], ...
                            'FaceColor',[1 1 1]);
                        annotation('textbox', [0.205 + 0.14765*(n-2) 0.875 - (p-1)*0.475 .071 .05], ...
                            'fontweight','bold','String', strrep(n_sim{n},'_',' '));
                    end; alpha(0.3);
                    xlim([0 max(x)]); ylim([0 max(Y)]); title(['Weight Change ' strrep(f_names{p},'_',' ')]);
                    if(p==1); xticklabels([]); else; xlabel('Time (seconds)'); end
                    ylabel('Synaptic Efficacy'); yticks(0:20:par{k}{i,j}{t}.weight_max); box on;
                    % save weight change figure
                    legend(p_h,titles(1+(p-1)*2:2*p),'location','southeast');
                end
                saveas(f(1),[filename '/Analysis/Weight Change/Weight_Change.tiff']); close(f(1));
            end
        end
    end
end

%% TIME FREQUENCY ANALYSIS
if(sum(cellfun(@sum,strfind(type,'TFA')))>0)
    %% DECLARATIONS
    if(length(vars.v1_vars) > 1 && strcmp('pos_LR', vars.v1_parameter)~=1)
        n_power = {'LowF', 'Theta', 'Alpha'}; f_power = {[3 30], [2 6], [8 12]}; 
    else
        n_power = {'Theta', 'Alpha'}; f_power = {[2 6], [8 12]}; 
    end
    
    
    n_stim = {'pre','post'}; 
    TFA_par.SR = 1000; 
    
    %% EXTRACT, ANALYSE AND SAVE DATA FOR TIME FREQUENCY ANALYSIS
    if(exist([vars.test_f '/Analysis/TFA Data'],'dir')~=7)
        mkdir([vars.test_f '/Analysis/TFA Data']);
        h1 = waitbar(0, 'Parameter Progression', 'Units', 'normalized', 'Position', [0.5 0.55 0.2 0.1]);
       for k=1:v3_L % loop variable 3
            if(isempty(vars.v3_vars)~=1)
                if(ischar(vars.v3_vars{k})~=1); v3 = num2str(vars.v3_vars{k}); 
                else; v3 = vars.v3_vars{k}; end
                pow_n1 = [vars.v3_name v3];
            else; v3 = []; pow_n1 = [];
            end; filename_1 = [vars.test_f '/' vars.v3_name v3];
            for i=1:v1_L % loop variable 1
                if(isempty(vars.v1_vars)~=1)
                    if(ischar(vars.v1_vars{i})~=1); v1 = num2str(vars.v1_vars{i}); 
                    else; v1 = vars.v1_vars{i}; end
                    pow_n2 = [vars.v1_name v1];
                else; v1 = []; pow_n2 = [];
                end; filename_2 = [filename_1 '/' vars.v1_name v1];
                for j=1:v2_L % loop variable 2
                    if(isempty(vars.v2_vars)~=1) 
                        if(ischar(vars.v2_vars{j})~=1); v2 = num2str(vars.v2_vars{j}); 
                        else; v2 = vars.v2_vars{j}; end
                        filename_3 = [filename_2 '_' vars.v2_name v2];
                        pow_n3 = [vars.v2_name v2];
                    else; filename_3 = filename_2; pow_n3 = [];
                    end
                    
                    pow_n = [pow_n1 pow_n2 pow_n3];
                    if(isempty(pow_n)==1); pow_n = 'TFA data'; end
                    
                    for t=1:vars.trials % loop trials
                        HPI = par{k}{i,j}{t}.Hip_per_item; NCPI = par{k}{i,j}{t}.NC_per_item;
                        n_NC = NCPI * par{k}{i,j}{t}.n_Items; n_Hip = par{k}{i,j}{t}.Hip_per_item * par{k}{i,j}{t}.n_Items; 
                        spike_sep = {[0,NCPI],[NCPI,NCPI*2], [n_NC,n_NC+HPI],[n_NC+HPI,n_NC+HPI*2],[0, n_NC],[n_NC, n_NC + n_Hip]};
                        
                        filename = [filename_3 '/T' int2str(t)]; % extend filename
                        waitbar(((k-1)*vars.trials*v1_L*v2_L...
                            + (i-1)*v2_L*vars.trials + ...
                            (j-1)*vars.trials + t) / ...
                            (vars.trials * v3_L * v1_L * v2_L),h1); % increment progress bar
                        sim_length = cell(1,N); spikes_all = [];
                        
                        for n=1:N % loop through simulation (put N+1 for all through simulation)
                            if(n<=N)
                                load([filename '/Data/' par{k}{i,j}{t}.sim_order_n{n} '.mat'])
                                TFA_par.sim_length = data.sim_stats.voltmeter(end,2);
                                spikes = data.sim_stats.spike_detector;
                                spikes_t = [spikes(:,1) spikes(:,2) + sum(cellfun(@sum,sim_length))];
                                spikes_all = [spikes_all; spikes_t]; sim_length{n} = TFA_par.sim_length;
                            else % use accumalitive data for whole sim
                                sim_length{N+1} = sum(cellfun(@sum,sim_length)); 
                                TFA_par.sim_length = sim_length{N+1}; spikes = spikes_all;
                            end
                            for s = 1:length(spike_sep) % loop through NC & Hip groups
                                % extract spike events of specific group
                                spikes_t = spikes(spikes(:,1) > spike_sep{s}(1),:); 
                                spikes_t = spikes_t(spikes_t(:,1) <= spike_sep{s}(2),:);
                                spikes_t(:,1) = spikes_t(:,1) - spike_sep{s}(1);
                                ns = length(spike_sep{s}(1)+1:spike_sep{s}(2));
                                if(strfind(n_sub_g{s},'NC')==1)
                                    r_phase = round(par{k}{i,j}{t}.NC_r_phase_n(n) * 1000); phase_end = 100; % alpha phase diff
                                elseif(strfind(n_sub_g{s},'Hip')==1)
                                    r_phase = round(par{k}{i,j}{t}.Hip_r_phase_n(n) * 1000); phase_end = 250; % theta phase diff
                                end
                                t_pre = [900 100]; t_post = [100 900];
                                for p=1:length(n_power)
                                    %% create local-field-potential based on spike data of group
                                    [lfp, pow_rc] = create_LFP(spikes_t, ns, f_power{p} ,TFA_par.sim_length, par{k}{i,j}{t}.pre_stim_length);
                            
                                    if(t==1); TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).RC = zeros(1, TFA_par.sim_length + 1); end
                                    TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).RC = TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).RC + pow_rc;
                                    if(t==vars.trials); TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).RC = ...
                                    TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).RC / vars.trials; end
                            
                                    %% initialise, increment & average TFA data structure over trials
                                    TFA_par.UBF = f_power{p}(1); TFA_par.OBF = f_power{p}(2);
                                    if(TFA_par.UBF<=30); TFA_par.Gamma=0.5; else; TFA_par.Gamma=pi*2; end
                                    [T, ~] = GaborFilter(lfp,TFA_par); pow = abs(T); % TFA analysis using LFP
                                    freq{n}{p} = TFA_par.UBF:(TFA_par.OBF - TFA_par.UBF)/(size(pow,1)-1):TFA_par.OBF;

                                    if(t==1) % initialise power data structures
                                        TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).pow = ...
                                            zeros(size(pow,1), TFA_par.sim_length + 1);
                                        if(length(vars.v1_vars)>1)
                                            if(i==1); for st = 1:2; v1_surfs.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).([n_stim{st}]){k}{j} = ...
                                                    zeros(v1_L, length(freq{n}{p}), vars.trials); end; end
                                        end
                                        TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).means = zeros(vars.trials, TFA_par.sim_length + 1);
                                    end
                                    TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).pow = ...
                                        TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).pow + pow; % increment
                                    TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).means(t,:) = mean(pow,1);
                                    
                                    if(length(vars.v1_vars)>1) % store data for multiple parameters analysis
                                        if(strfind(n_sim{n},'idle')==1); pre_stim = par{k}{i,j}{t}.idling_length; else; pre_stim = par{k}{i,j}{t}.pre_stim_length; end
                                        v1_surfs.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).pre{k}{j}(i,:,t) = ...
                                            transpose(mean(pow(:,t_pre(1)+1:pre_stim-t_pre(2)),2));
                                        v1_surfs.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).post{k}{j}(i,:,t) = ...
                                            transpose(mean(pow(:,pre_stim+t_post(1)+1:end-t_post(2)),2));
                                    end

                                    if(t==vars.trials) % average power data over trials
                                        TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).pow = ...
                                            TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).pow / vars.trials;
                                        if(v1_L <= 1) % BOOTSTRAP TFA LINES
                                            [TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).CI_L, ...
                                                TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).CI_U, ...
                                                TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).M] = ...
                                                bootstrap(TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).means, n_BS,[]);
                                        else % for speed reasons, only bootstrap for TFA plots in Figures 5-6
                                            TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).M = mean(TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).means);
                                            TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).CI_L = min(TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).means);
                                            TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).CI_U = max(TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).means);
                                        end
                                        TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]) = ...
                                            rmfield(TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]), 'means');
                                    end
                                    
                                    %% take random phase from LFPs after TFA has been performed & save in data structures
                                    lfp = lfp(r_phase + 1 : end-(phase_end - r_phase));
                                    if(t==1); TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).LFP = zeros(1, TFA_par.sim_length - phase_end + 1); end
                                    TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).LFP = TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).LFP + lfp;
                                    if(t==vars.trials); TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).LFP = ...
                                    TFA_data.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).LFP / vars.trials; end
                                    
                                    %% store power at recall (P/NP power to NP/P stimulus after-learning)
                                    if(strcmp(n_sim{n},'NP_AL')==1 || strcmp(n_sim{n},'P_AL')==1 || strcmp(n_sim{n},'P_BL')==1 || strcmp(n_sim{n},'NP_BL')==1)
                                        recall.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).post{k}{j}{i}(t,1) = ...
                                            mean(mean(pow(min(find(floor(freq{n}{p})==f_power{p}(1))):max(find(ceil(freq{n}{p})==f_power{p}(2))), ...
                                            par{k}{i,j}{t}.pre_stim_length+126:end-125)));
                                        recall.([n_sim{n}]).([n_sub_g{s}]).([n_power{p}]).pre{k}{j}{i}(t,1) = ...
                                            mean(mean(pow(min(find(floor(freq{n}{p})==f_power{p}(1))):max(find(ceil(freq{n}{p})==f_power{p}(2))), ...
                                            126:par{k}{i,j}{t}.pre_stim_length-125)));
                                    end
                                end
                            end
                        end
                    end
                    save([vars.test_f '/Analysis/TFA Data/Pow_' pow_n '.mat'],'TFA_data','-v7.3'); clearvars('TFA_data');
                end
            end
       end
       close(h1); 
       if(length(vars.v1_vars)>1); TFA_analysis.v1_surfs = v1_surfs; end
       TFA_analysis.freq = freq; TFA_analysis.recall = recall;
       save([vars.test_f '/Analysis/TFA Data/TFA_analysis.mat'],'TFA_analysis','-v7.3')
    else
        load([vars.test_f '/Analysis/TFA Data/TFA_analysis.mat']); 
    end

    %% TIME FREQUENCY ANALYSIS EACH SINGLE SET OF PARAMETERS
    for k=1:v3_L % loop variable 3
       if(isempty(vars.v3_vars)~=1)
            if(ischar(vars.v3_vars{k})~=1); v3 = num2str(vars.v3_vars{k}); 
            else; v3 = vars.v3_vars{k}; end
            pow_n1 = [vars.v3_name v3];
        else; v3 = []; pow_n1 = [];
        end; filename_1 = [vars.test_f '/' vars.v3_name v3];
        for i=1:v1_L % loop variable 1
            if(isempty(vars.v1_vars)~=1) 
                if(ischar(vars.v1_vars{i})~=1); v1 = num2str(vars.v1_vars{i}); 
                else; v1 = vars.v1_vars{i}; end
                pow_n2 = [vars.v1_name v1];
            else; v1 = []; pow_n2 = [];
            end; filename_2 = [filename_1 '/' vars.v1_name v1];
            for j=1:v2_L % loop variable 2
                if(isempty(vars.v2_vars)~=1) 
                    if(ischar(vars.v2_vars{j})~=1); v2 = num2str(vars.v2_vars{j}); 
                    else; v2 = vars.v2_vars{j}; end
                    filename_3 = [filename_2 '_' vars.v2_name v2];
                    pow_n3 = [vars.v2_name v2];
                else; filename_3 = filename_2; pow_n3 = [];
                end
                pow_n = [pow_n1 pow_n2 pow_n3];
                if(isempty(pow_n)==1); pow_n = 'TFA data'; end
                if(v1_L<=1)
                    mkdir([filename_3 '/Analysis/Time Frequency Analysis'])
                    load([vars.test_f '/Analysis/TFA Data/Pow_' pow_n '.mat']); % load TFA data from previous function   
                    %% EXTRACT PLOT PARAMETERS
                    pow = struct; rec_BL = struct; rec_AL = struct;
                    for b = 1:length(n_group)
                        for p=1:length(n_power) % extract Theta/Alpha power for each stage
                            pow.([n_group{b}]).([n_power{p}]).rec_BL = (TFA_data.NP_BL.([n_group{b} '_P']).([n_power{p}]).pow + ... % recall BL
                                TFA_data.P_BL.([n_group{b} '_NP']).([n_power{p}]).pow) / 2;
                            pow.([n_group{b}]).([n_power{p}]).resp_BL = (TFA_data.NP_BL.([n_group{b} '_NP']).([n_power{p}]).pow + ... % recognition BL
                                TFA_data.P_BL.([n_group{b} '_P']).([n_power{p}]).pow) / 2;
                            pow.([n_group{b}]).([n_power{p}]).DL = (TFA_data.DL.([n_group{b} '_NP']).([n_power{p}]).pow + ... % DL
                                TFA_data.DL.([n_group{b} '_P']).([n_power{p}]).pow) / 2; 
                            pow.([n_group{b}]).([n_power{p}]).rec_AL = (TFA_data.NP_AL.([n_group{b} '_P']).([n_power{p}]).pow + ... % recall AL
                                TFA_data.P_AL.([n_group{b} '_NP']).([n_power{p}]).pow) / 2;
                            pow.([n_group{b}]).([n_power{p}]).resp_AL = (TFA_data.NP_AL.([n_group{b} '_NP']).([n_power{p}]).pow + ... % recognition AL
                                TFA_data.P_AL.([n_group{b} '_P']).([n_power{p}]).pow) / 2;
                            
                            rec_BL.([n_group{b}]).([n_power{p}]) = [(TFA_data.NP_BL.([n_group{b} '_P']).([n_power{p}]).CI_U + ... % recall BL
                                TFA_data.P_BL.([n_group{b} '_NP']).([n_power{p}]).CI_U) / 2; ...
                                (TFA_data.NP_BL.([n_group{b} '_P']).([n_power{p}]).M + ... 
                                TFA_data.P_BL.([n_group{b} '_NP']).([n_power{p}]).M) / 2; ...
                                (TFA_data.NP_BL.([n_group{b} '_P']).([n_power{p}]).CI_L + ... 
                                TFA_data.P_BL.([n_group{b} '_NP']).([n_power{p}]).CI_L) / 2];
                            rec_AL.([n_group{b}]).([n_power{p}]) = [(TFA_data.NP_AL.([n_group{b} '_P']).([n_power{p}]).CI_U + ... % recall BL
                                TFA_data.P_AL.([n_group{b} '_NP']).([n_power{p}]).CI_U) / 2; ...
                                (TFA_data.NP_AL.([n_group{b} '_P']).([n_power{p}]).M + ... 
                                TFA_data.P_AL.([n_group{b} '_NP']).([n_power{p}]).M) / 2; ...
                                (TFA_data.NP_AL.([n_group{b} '_P']).([n_power{p}]).CI_L + ... 
                                TFA_data.P_AL.([n_group{b} '_NP']).([n_power{p}]).CI_L) / 2];
                        end
                    end
                    %% PLOT DATA
                    leg_loc = {'southwest','west'};
                    for b=1:2
                        %% PLOT TFA HEATMAPS
                        f = figure; set(f, 'Position', [0 0 1300 350]); % POWER TFA PLOTS
                        if(strcmp(n_group{b},'NC')==1); p = 2; else; p = 1; end
                        freq = TFA_analysis.freq{2}{p}; colormap('hot');
                        T = -1:1/1000:1; x = par{k}{i,j}{1}.pre_stim_length - 1000 : par{k}{i,j}{1}.pre_stim_length + 1000;  
                        % find min/max of group through simulation stages
                        z_max = max(max([max(TFA_data.NP_BL.([n_group{b} '_P']).([n_power{p}]).pow(:,x)); ...
                            max(TFA_data.NP_BL.([n_group{b} '_NP']).([n_power{p}]).pow(:,x)); max(TFA_data.P_BL.([n_group{b} '_P']).([n_power{p}]).pow(:,x)); ...
                            max(TFA_data.P_BL.([n_group{b} '_NP']).([n_power{p}]).pow(:,x)); max(TFA_data.DL.([n_group{b} '_P']).([n_power{p}]).pow(:,x)); ...
                            max(TFA_data.DL.([n_group{b} '_NP']).([n_power{p}]).pow(:,x)); max(TFA_data.P_AL.([n_group{b} '_P']).([n_power{p}]).pow(:,x)); ...
                            max(TFA_data.P_AL.([n_group{b} '_P']).([n_power{p}]).pow(:,x))]));
                        z_min = min(min([min(TFA_data.NP_BL.([n_group{b} '_P']).([n_power{p}]).pow(:,x)); ...
                            min(TFA_data.NP_BL.([n_group{b} '_NP']).([n_power{p}]).pow(:,x)); min(TFA_data.P_BL.([n_group{b} '_P']).([n_power{p}]).pow(:,x)); ...
                            min(TFA_data.P_BL.([n_group{b} '_NP']).([n_power{p}]).pow(:,x)); min(TFA_data.DL.([n_group{b} '_P']).([n_power{p}]).pow(:,x)); ...
                            min(TFA_data.DL.([n_group{b} '_NP']).([n_power{p}]).pow(:,x)); min(TFA_data.P_AL.([n_group{b} '_P']).([n_power{p}]).pow(:,x)); ...
                            min(TFA_data.P_AL.([n_group{b} '_P']).([n_power{p}]).pow(:,x))]));
                        
                        subplot(2,5,[1 6]); imagesc(T, freq, pow.([n_group{b}]).([n_power{p}]).rec_BL(:,x)); 
                        title('Recall BL'); ylabel('Frequency (Hz)'); caxis([z_min z_max]); ax = gca; ax.FontSize = 14;
                        annotation('rectangle', [0.13 0.825 .03 .1], 'FaceColor',[1 1 1]);
                        annotation('textbox', [0.13 0.825 1 0.1], 'String', 'Ai', 'EdgeColor', 'none', 'FontSize', 14,'fontweight','bold')
                        set(gca,'linewidth',5,'XColor','k','YColor','k');
                        
                        subplot(2,5,[2 7]); imagesc(T, freq, pow.([n_group{b}]).([n_power{p}]).resp_BL(:,x)); 
                        title('Recognition BL'); caxis([z_min z_max]); ax = gca; ax.FontSize = 14;
                        annotation('rectangle', [0.2935 0.825 .03 .1], 'FaceColor',[1 1 1]);
                        annotation('textbox', [0.2935 0.825 1 0.1], 'String', 'Aii', 'EdgeColor', 'none', 'FontSize', 14,'fontweight','bold')
                        
                        subplot(2,5,[3 8]); imagesc(T, freq, pow.([n_group{b}]).([n_power{p}]).DL(:,x)); 
                        title('DL'); xlabel('Time (s)'); caxis([z_min z_max]); ax = gca; ax.FontSize = 14;
                        annotation('rectangle', [0.4555 0.825 .03 .1], 'FaceColor',[1 1 1]);
                        annotation('textbox', [0.4555 0.825 1 0.1], 'String', 'Aiii', 'EdgeColor', 'none', 'FontSize', 14,'fontweight','bold')
                        
                        subplot(2,5,[4 9]); imagesc(T, freq, pow.([n_group{b}]).([n_power{p}]).rec_AL(:,x)); 
                        title('Recall AL'); caxis([z_min z_max]); ax = gca; ax.FontSize = 14;
                        annotation('rectangle', [0.618 0.825 .03 .1], 'FaceColor',[1 1 1]);
                        annotation('textbox', [0.618 0.825 1 0.1], 'String', 'Aiv', 'EdgeColor', 'none', 'FontSize', 14,'fontweight','bold')
                        set(gca,'linewidth',5,'XColor',[0 0.6 0],'YColor',[0 0.6 0]);
                        
                        subplot(2,5,[5 10]); imagesc(T, freq, pow.([n_group{b}]).([n_power{p}]).resp_AL(:,x)); 
                        title('Recognition AL'); caxis([z_min z_max]); ax = gca; ax.FontSize = 14;
                        annotation('rectangle', [0.782 0.825 .03 .1], 'FaceColor',[1 1 1]);
                        annotation('textbox', [0.782 0.825 1 0.1], 'String', 'Av', 'EdgeColor', 'none', 'FontSize', 14,'fontweight','bold')
                        saveas(f,[filename_3 '/Analysis/Time Frequency Analysis/' n_group{b} '_TFA.tiff']); close(f);
                        
                        %% PLOT POWER DIFFERENCES
                        f = figure; set(f, 'Position', [0 0 1300 350]); % POWER TFA PLOTS
                        x_pre = x(1:floor(length(x)/2)); if(strcmp(n_group{b},'NC')==1); p = 2; else; p = 1; end
                        if(b==1); yL = [-50 20]; elseif(b==2); yL = [-20 70]; end
                        
                        subplot(2,5,(1:3)); hold on
                        fill([T fliplr(T)], [rec_BL.([n_group{b}]).([n_power{p}])(1,x) fliplr(rec_BL.([n_group{b}]).([n_power{p}])(3,x))], [0.2 0.2 0.2],...
                            'edgecolor','k'); 
                        h(1) = plot(T, rec_BL.([n_group{b}]).([n_power{p}])(2,x), 'linewidth',2, 'linestyle','--', 'color','k');
                        fill([T fliplr(T)], [rec_AL.([n_group{b}]).([n_power{p}])(1,x) fliplr(rec_AL.([n_group{b}]).([n_power{p}])(3,x))], [0.2 0.8 0.2],...
                            'edgecolor',[0.1 0.5 0.1]); 
                        h(2) = plot(T, rec_AL.([n_group{b}]).([n_power{p}])(2,x), 'linewidth',2, 'linestyle','-', 'color',[0.1 0.5 0.1]);                    
                        title([n_group{b} ' ' n_power{p} ' recall power']); ax = gca; ax.FontSize = 14;
                        legend([h(1) h(2)], 'Before learning', 'After learning','location',leg_loc{b}); ylabel('Power');
                        text(0.02,0.98,'Bi','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14,'Fontweight','bold'); alpha(0.4);
                        xticklabels([]);
                        
                        subplot(2,5,(1:3)+5); hold on;
                        BL_change = ((rec_BL.([n_group{b}]).([n_power{p}])(:,x) - mean(rec_BL.([n_group{b}]).([n_power{p}])(2,x_pre))) / ...
                            mean(rec_BL.([n_group{b}]).([n_power{p}])(2,x_pre))) * 100;
                        AL_change = ((rec_AL.([n_group{b}]).([n_power{p}])(:,x) - mean(rec_AL.([n_group{b}]).([n_power{p}])(2,x_pre))) / ...
                            mean(rec_AL.([n_group{b}]).([n_power{p}])(2,x_pre))) * 100;
                        fill([T fliplr(T)], [BL_change(1,:) fliplr(BL_change(3,:))], [0.2 0.2 0.2], 'edgecolor','k'); 
                        plot(T, BL_change(2,:), 'linewidth',2, 'linestyle','--', 'color','k');
                        fill([T fliplr(T)], [AL_change(1,:) fliplr(AL_change(3,:))], [0.2 0.8 0.2],'edgecolor',[0.1 0.5 0.1]); 
                        plot(T, AL_change(2,:), 'linewidth',2, 'linestyle','-', 'color',[0.1 0.5 0.1]);
                        title([n_group{b} ' ' n_power{p} ' recall pre-post % change']); ax = gca; ax.FontSize = 14;
                        ylabel('% change'); xlabel('Time (s)'); ylim(yL);
                        text(0.02,0.98,'Bii','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14,'Fontweight','bold'); alpha(0.4);
                        
                        subplot(2,5,4:5); hold on
                        fill([T fliplr(T)], [rec_BL.([n_group{3-b}]).([n_power{p}])(1,x) fliplr(rec_BL.([n_group{3-b}]).([n_power{p}])(3,x))],[0.2 0.2 0.2],...
                            'edgecolor','k'); 
                        plot(T, rec_BL.([n_group{3-b}]).([n_power{p}])(2,x), 'linewidth',2, 'linestyle','--', 'color','k');
                        fill([T fliplr(T)], [rec_AL.([n_group{3-b}]).([n_power{p}])(1,x) fliplr(rec_AL.([n_group{3-b}]).([n_power{p}])(3,x))],[0.2 0.8 0.2],...
                            'edgecolor',[0.1 0.5 0.1]); 
                        plot(T, rec_AL.([n_group{3-b}]).([n_power{p}])(2,x), 'linewidth',2, 'linestyle','-', 'color',[0.1 0.5 0.1]);                    
                        title([n_group{3-b} ' ' n_power{p} ' recall power']); xticklabels([]); ax = gca; ax.FontSize = 14;
                        text(0.02,0.98,'Ci','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14,'Fontweight','bold'); alpha(0.4);
                        
                        subplot(2,5,(4:5)+5); hold on;
                        BL_change = ((rec_BL.([n_group{3-b}]).([n_power{p}])(:,x) - mean(rec_BL.([n_group{3-b}]).([n_power{p}])(2,x_pre))) / ...
                            mean(rec_BL.([n_group{3-b}]).([n_power{p}])(2,x_pre))) * 100;
                        AL_change = ((rec_AL.([n_group{3-b}]).([n_power{p}])(:,x) - mean(rec_AL.([n_group{3-b}]).([n_power{p}])(2,x_pre))) / ...
                            mean(rec_AL.([n_group{3-b}]).([n_power{p}])(2,x_pre))) * 100;
                        fill([T fliplr(T)], [BL_change(1,:) fliplr(BL_change(3,:))], [0.2 0.2 0.2],'edgecolor','k'); 
                        plot(T, BL_change(2,:), 'linewidth',2, 'linestyle','--', 'color','k');
                        fill([T fliplr(T)], [AL_change(1,:) fliplr(AL_change(3,:))], [0.2 0.8 0.2],'edgecolor',[0.1 0.5 0.1]); 
                        plot(T, AL_change(2,:), 'linewidth',2, 'linestyle','-', 'color',[0.1 0.5 0.1]);
                        title([n_group{3-b} ' ' n_power{p} ' recall pre-post % change']);  xlabel('Time (s)'); ylim(yL); ax = gca; ax.FontSize = 14;
                        text(0.02,0.98,'Cii','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14,'Fontweight','bold'); alpha(0.4);
                        
                        saveas(f,[filename_3 '/Analysis/Time Frequency Analysis/' n_group{b} '_power.tiff']); close(f);
                    end
                end
            end
        end
        
        %% TFA MULTIPLE PARAMETER ANALYSIS
        if(length(vars.v1_vars)>1)
            load([vars.test_f '/Analysis/weight-change-data.mat']); % load weight change data
            load([vars.test_f '/Analysis/firing-rate-data.mat']);
            %x = cell2mat(vars.v1_vars);
            x = cell2mat(vars.v1_vars);
            if(exist([filename_1 '/Analyse Variables'],'dir')~=7);  mkdir([filename_1 '/Analyse Variables']); end

            for j=1:v2_L % loop variable 2 
                if(strcmp(vars.v1_parameter, 'pos_LR')~=1)
                    %% EXTRACT 3D DATA IF VARYING 2 PARAMETERS
                    a_x = [8 12]; t_x = [3 5];
                    h1 = waitbar(0, 'Variables', 'Units', 'normalized', 'Position', [0.5 0.7 0.2 0.1]);
                   for s=1:length(n_sub_g) % loop through NC & Hip spikes
                      for n=2:N-1
                         waitbar(((n-1)+(s-1)*5)/(5*6),h1);
                         % PRE -> POST POWER CHANGE
                         z_pre = TFA_analysis.v1_surfs.([n_sim{n}]).([n_sub_g{s}]).LowF.pre{k}{j};
                         z_post = TFA_analysis.v1_surfs.([n_sim{n}]).([n_sub_g{s}]).LowF.post{k}{j};
                         z = (z_post ./ z_pre) * 100 - 100;
                         y = TFA_analysis.freq{n}{1};
                         % ALPHA
                         z_a = transpose(squeeze(mean(z(:,find(y==a_x(1)):find(y==a_x(2)),:),2)));
                         if(vars.trials>1)
                            [CI_L, CI_U, M] = bootstrap(z_a, n_BS,[]);
                         else
                             CI_L = z_a; CI_U = z_a; M = z_a;
                         end
                         pow_y.([n_sub_g{s}]).LowF.([n_sim{n}]).Alpha = [CI_L; M; CI_U];

                         % THETA
                         z_t = transpose(squeeze(mean(z(:,find(y==t_x(1)):find(y==t_x(2)),:),2)));
                          if(vars.trials>1)
                            [CI_L, CI_U, M] = bootstrap(z_t, n_BS,[]);
                         else
                             CI_L = z_t; CI_U = z_t; M = z_t;
                         end
                         pow_y.([n_sub_g{s}]).LowF.([n_sim{n}]).Theta = [CI_L; M; CI_U];
                      end
                   end
                   close(h1);

                   %% FIGURE 7 PLOT
                   n_S = {'Encoding', 'Recall'}; p_col = {[1 0.3 0.3], [0.3 0.3 1]}; x_marker = {'o', 'x', '^'};
                   sm_par_pow = 1e-7; sm_par_W = 1e-13; y_min = []; y_max = []; n_power = {'Alpha', 'Theta'};
                   label = {'A','B','C','D','E','F'}; label_i = {'i','ii'};
                   left_color = [0.5 0 0.5]; right_color = [0 0 0];
                   f(1) = figure(1); set(f(1), 'Position', [0 0 1170 585]); colormap('Jet')
                   set(f(1),'defaultAxesColorOrder',[left_color; right_color]);
                   for s = 1:length(n_S)
                       % extract P<->NP weights
                       if(strcmp('Encoding', n_S{s}) == 1)
                           w = ([WC.AL_W.NP_P{k}{:,j}]+[WC.AL_W.P_NP{k}{:,j}])/2;
                       elseif(strcmp('Recall', n_S{s}) == 1)
                           w = sqrt((([WC.REC_W.NP_P{k}{:,j}]+[WC.REC_W.P_NP{k}{:,j}])/2 - ([WC.AL_W.NP_P{k}{:,j}]+[WC.AL_W.P_NP{k}{:,j}])/2).^2);
                       end
                       % adapt to % weight change (0-max) per second
                       w = (w/(par{k}{i,j}{t}.stim_length/1000))/(par{k}{i,j}{t}.weight_max/2/100); W = [];
                       [W(1,:), W(3,:), W(2,:)] = bootstrap(w, n_BS,[]);
                   
                       % NC POWER PLOT
                       subplot(4, 6, 1 + (s-1)*6*length(n_S)); 
                       y = TFA_analysis.freq{4}{1}; 
                       if(strcmp('Encoding', n_S{s}) == 1)
                           z = mean(TFA_analysis.v1_surfs.DL.NC_All.LowF.post{k}{j},3);
                       elseif(strcmp('Recall', n_S{s}) == 1)
                           z = mean(TFA_analysis.v1_surfs.NP_AL.NC_P.LowF.post{k}{j},3) + ...
                               mean(TFA_analysis.v1_surfs.P_AL.NC_NP.LowF.post{k}{j},3);
                       end
                       pcolor(x,y,log10(transpose(z))); shading interp; set(gca,'XScale','log'); %colorbar;
                       annotation('rectangle', [0.13 0.845-(s-1)*0.4375 .0308 .0615], 'FaceColor',[1 1 1]);
                       annotation('textbox', [0.13 0.845-(s-1)*0.4375 .0308 .0615],'string',[label{1+(s-1)*3} 'i'],'FontSize',14,'Fontweight','bold');
                       ylabel('Frequency (Hz)');  ytickangle(90); xticklabels([]);
                       title(['NC ' n_S{s}]); 

                       % HIP POWER PLOT
                       subplot(4, 6, 7 + (s-1)*6*length(n_S)); 
                       y = TFA_analysis.freq{4}{1}; 
                       if(strcmp('Encoding', n_S{s}) == 1)
                           z = mean(TFA_analysis.v1_surfs.DL.Hip_All.LowF.post{k}{j},3);
                       elseif(strcmp('Recall', n_S{s}) == 1)
                           z = mean(TFA_analysis.v1_surfs.NP_AL.Hip_P.LowF.post{k}{j},3) + ...
                               mean(TFA_analysis.v1_surfs.P_AL.Hip_NP.LowF.post{k}{j},3);
                       end
                       pcolor(x,y,log10(transpose(z))); shading interp; set(gca,'XScale','log'); %colorbar;
                       annotation('rectangle', [0.13 0.625-(s-1)*0.4375 .0308 .0615], 'FaceColor',[1 1 1]);
                       annotation('textbox', [0.13 0.625-(s-1)*0.4375 .0308 .0615],'string',[label{1+(s-1)*3} 'ii'],'FontSize',14,'Fontweight','bold');
                       ylabel('Frequency (Hz)'); ytickangle(90);
                       if(s==2); xlabel(regexprep(vars.v1_parameter,'_',' ','emptymatch')); end
                       xticks([10^0 10^2 10^4 10^6]); xticklabels({'10^0', '10^2', '10^4', '10^6'}); 
                       title(['Hip ' n_S{s}]);

                       % WEIGHTS 
                       subplot(4, 6, [2:5 8:11] + (s-1)*6*length(n_S))
                       WU = fit( transpose(x), transpose(W(3,:)),'smoothingspline','smoothingparam',sm_par_W); % UPPER CI
                       WL = fit( transpose(x), transpose(W(1,:)),'smoothingspline','smoothingparam',sm_par_W); % LOWER CI
                       WM = fit( transpose(x), transpose(W(2,:)),'smoothingspline','smoothingparam',sm_par_W); % MEANS
                       yyaxis right; hold on; 
                       xlim([min(x) max(x)]); set(gca,'XScale','log');
                       fill([x fliplr(x)], [transpose(WL(x)) transpose(flipud(WU(x)))],[0.3 0.3 1],'linestyle','none','marker','none'); alpha(0.1);
                       h(1 + 12*(s-1)) = semilogx(x, WU(x)); set(h(1 + 12*(s-1)),'color','k','linestyle','-','marker','none'); 
                       h(2 + 12*(s-1)) = semilogx(x, WL(x)); set(h(2 + 12*(s-1)),'color','k','linestyle','-','marker','none');
                       h(3 + 12*(s-1)) = semilogx(x, WM(x)); set(h(3 + 12*(s-1)),'color','k','linestyle','-','marker','none','linewidth',2); 
                       ylim([0 max(W(3,:))]);
                       y_min(s,1,2) = min(min(W)); y_max(s,1,2) = max(max(W));

                       % THETA
                       if(strcmp('Encoding', n_S{s}) == 1)
                           y = pow_y.Hip_All.LowF.DL.Theta; % During learning
                       elseif(strcmp('Recall', n_S{s}) == 1)
                           y = (pow_y.Hip_NP.LowF.P_AL.Theta + pow_y.Hip_P.LowF.NP_AL.Theta)/2; % Recall
                       end
                       TU = fit(transpose(x),transpose(y(3,:)),'smoothingspline','smoothingparam',sm_par_pow);
                       TL = fit(transpose(x),transpose(y(1,:)),'smoothingspline','smoothingparam',sm_par_pow);
                       TM = fit(transpose(x),transpose(y(2,:)),'smoothingspline','smoothingparam',sm_par_pow);
                       yyaxis left; hold on; 
                       fill([x fliplr(x)], [transpose(TL(x)) transpose(flipud(TU(x)))],[0.3 0.3 1],'linestyle','none','marker','none'); alpha(0.1);
                       h(4 + 12*(s-1)) = semilogx(x,TU(x)); set(h(4 + 12*(s-1)),'color',[0.3 0.3 1], 'linestyle','-','marker','none');
                       h(5 + 12*(s-1)) = semilogx(x, TL(x)); set(h(5 + 12*(s-1)),'color',[0.3 0.3 1], 'linestyle','-','marker','none');
                       h(6 + 12*(s-1)) = semilogx(x, TM(x)); set(h(6 + 12*(s-1)),'color',[0.3 0.3 1], 'linestyle','-','marker','none','linewidth',2); 
                       y_min(s,1,1) = min(min(y)); y_max(s,1,1) = max(max(y));

                       % ALPHA
                       if(strcmp('Encoding', n_S{s}) == 1)
                           y = pow_y.NC_All.LowF.DL.Alpha; % During learning
                       elseif(strcmp('Recall', n_S{s}) == 1)
                           y = (pow_y.NC_NP.LowF.P_AL.Alpha + pow_y.NC_P.LowF.NP_AL.Alpha)/2; % Recall
                       end
                       AU = fit(transpose(x),transpose(y(3,:)),'smoothingspline','smoothingparam',sm_par_pow);
                       AL = fit(transpose(x),transpose(y(1,:)),'smoothingspline','smoothingparam',sm_par_pow);
                       AM = fit(transpose(x),transpose(y(2,:)),'smoothingspline','smoothingparam',sm_par_pow);
                       fill([x fliplr(x)], [transpose(AL(x)) transpose(flipud(AU(x)))],[1 0.3 0.3],'linestyle','none','marker','none'); alpha(0.1);
                       h(7 + 12*(s-1)) = semilogx(x, AU(x)); set(h(7 + 12*(s-1)),'color',[1 0.3 0.3],'linestyle','-','marker','none');
                       h(8 + 12*(s-1)) = semilogx(x, AL(x)); set(h(8 + 12*(s-1)),'color',[1 0.3 0.3],'linestyle','-','marker','none');
                       h(9 + 12*(s-1)) = semilogx(x, AM(x)); set(h(9 + 12*(s-1)),'color',[1 0.3 0.3],'linestyle','-','marker','none','linewidth',2); 
                       y_min(s,2,1) = min(min(y)); y_max(s,2,1) = max(max(y));

                       % LABELS AND ZERO LINE
                       plot(zeros(1, max(x)),'k:');
                       text(0.02,0.98,label{2+(s-1)*3},'Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14,'Fontweight','bold');
                       if(s==length(n_S)); xlabel(regexprep(vars.v1_parameter,'_',' ','emptymatch')); else; xticklabels([]); end
                       yyaxis left; ylabel('% change'); ytickangle(90);
                       yyaxis right;  ytickangle(270); %ylabel('% weight change/second');
                       title([n_S{s} ' v ' regexprep(vars.v1_parameter,'_',' ','emptymatch')]);
                       xlim([x(2) x(end)]); 
                       % plot lines to show which points are taken for I plots in graph below
                       X_i = [1 min(find(max(AM(x))==AM(x)==1)) min(find(max(TM(x))==TM(x)==1))];
                       X_i(X_i==1) = min(find(x>=10)); % find first > 10^1
                       X_i(X_i==length(x)) = length(x) - 1; % find one before last
                       if(X_i(2)==X_i(1)); X_i(2) = min(find(x>=10^3)); end
                       for a = 1:length(X_i)
                          h(9+a + 12*(s-1)) = plot(ones(101,1)*x(X_i(a)),0:100,'color',[0.1 .4 0.1],...
                              'linestyle','-','marker', x_marker{a},'linewidth',1); 
                          line_lab{a} = [int2str(floor(abs(x(X_i(a))) ./ 10.^floor(log10(abs(x(X_i(a))))))) 'x10^' int2str(floor(log10(abs(x(X_i(a)))+1)))];
                       end
                       legend([h(3 + 12*(s-1)) h(6 + 12*(s-1)) h(9 + 12*(s-1)) h(10 + 12*(s-1)) h(11 + 12*(s-1)) h(12 + 12*(s-1))],...
                           {'Weights', 'Hip Theta', 'NC Alpha',line_lab{1},line_lab{2},line_lab{3}},'position',[0.375  0.765-(s-1)*.45 0.065 0.1]);

                       plot_L = 1000; % ms 
                       g = {'NC', 'Hip'}; lfp_gap = [0.15 0.05];
                       l_styles = {'-','-','-'}; mean_NC = zeros(1,3); mean_HIP = zeros(1,3);
                       for p = 1:2
                           subplot(4, 6, 6 + (p-1)*6 + (s-1)*6*length(n_S)); hold on; % PLOT LFP ACTIVITY
                           text(0.02,0.98,[label{3+(s-1)*3} label_i{p}],'Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14,'Fontweight','bold');
                           for a = 1:length(X_i)
                               load([vars.test_f '/Analysis/TFA Data/Pow_' vars.v1_name int2str(x(X_i(a))) '.mat']);
                               if(strcmp('Encoding', n_S{s}) == 1)
                                   lfp = mean([TFA_data.DL.([g{p} '_P']).([n_power{p}]).LFP; TFA_data.DL.([g{p} '_NP']).([n_power{p}]).LFP],1);
                               elseif(strcmp('Recall', n_S{s}) == 1)
                                   lfp = mean([TFA_data.NP_AL.([g{p} '_P']).([n_power{p}]).LFP; TFA_data.P_AL.([g{p} '_NP']).([n_power{p}]).LFP],1);
                               end
                               phase_diff = (round(length(lfp)/1000)*1000+1) - length(lfp); % calculate different in phases for stimulus onset
                               lfp = lfp(par{k}{i,j}{t}.pre_stim_length+1 - phase_diff : par{k}{i,j}{t}.pre_stim_length - phase_diff+ plot_L);
                               if(a==1); x_a = 0; else; x_a = lfp_max(a-1); end
                               lfp = lfp + x_a + lfp_gap(p); lfp_min(a) = min(lfp); lfp_max(a) = max(lfp);
                               plot(lfp,'color', [0.1 0.4 0.1], 'linewidth',1,'linestyle','-','marker',x_marker{a},'markerindices',[1 1000]);
                           end
                           fill([1:phase_diff fliplr(1:phase_diff)], [ones(1, phase_diff)*-0.5 ones(1, phase_diff)*1.5],...
                               p_col{p}, 'linestyle','none','marker','none'); alpha(0.2); % box for possible SO period from varying phases
                           title([g{p} ' LFPs']);
                           ylim([min(lfp_min)-lfp_gap(p)/2 max(lfp_max)+lfp_gap(p)*1.75]); yticks([]);
                           if(p~=2); xticklabels([]); end
                           if(s==length(n_S) && p==2); xlabel('Time after SO (ms)'); end
                       end
                   end
                   
                   for s = 1:length(n_S)
                       subplot(2 * length(n_S), 6, [2:5 8:11] + (s-1)*6*length(n_S))
                       yyaxis left; ylim([min(min(y_min(:,:,1))) max(max(y_max(:,:,1)))]);
                       yyaxis right; ylim([0 ceil(max(max(y_max(:,:,2))))]);
                       for a=1:3; h(9+a + 12*(s-1)).MarkerIndices = round(1:floor(max(max(y_max(:,:,2))))/3-1:ceil(max(max(y_max(:,:,2))))); end
                   end
                   saveas(f(1),[filename_1 '/Analyse Variables/POWvW.tiff']); close(f(1));
                
                elseif(strcmp(vars.v1_parameter, 'pos_LR')==1)
                    %% FIGURE 8 PLOT
                    titles = {'Recognition','Recall'}; label = {'A','B'}; label_i = {'i','ii'}; sp = {1:3, 4:5};
                    x = cell2mat(vars.v1_vars);
                    weights_t = []; weights_m = []; weights_a = []; 
                    for i=1:v1_L
                        weights_t = [weights_t (WC.AL_W.P_NP{k}{i,j} + WC.AL_W.NP_P{k}{i,j})/2];
                        weights_m = [weights_m; mean(WC.AL_W.P_NP{k}{i,j}) mean(WC.AL_W.NP_P{k}{i,j})];
                        weights_a = [weights_a; WC.AL_W.P_NP{k}{i,j} WC.AL_W.NP_P{k}{i,j}];
                    end
                    W_M = mean(weights_t,1); 
                    if(vars.trials>1)
                        W_CI_L = W_M - var(weights_t); W_CI_U = W_M + var(weights_t); 
                    else
                        W_CI_L = W_M; W_CI_U = W_M; 
                    end
                    
                    for s = 1:2
                        %% FIND PRE - POST POWERS
                        f = figure; set(f, 'Position', [0 0 900 450]); 
                        for b = 1:length(n_group)
                            for p = 1:length(n_power)
                                pre = []; post = []; if(b==1); px = 3-p; else; px = p; end
                                for i=1:v1_L
                                    pre_BL = (TFA_analysis.recall.P_BL.([n_group{b} '_' n_syn_g{s}]).([n_power{px}]).pre{k}{j}{i} + ...
                                        TFA_analysis.recall.NP_BL.([n_group{b} '_' n_syn_g{3-s}]).([n_power{px}]).pre{k}{j}{i}) / 2;
                                    pre_AL = (TFA_analysis.recall.P_AL.([n_group{b} '_' n_syn_g{s}]).([n_power{px}]).pre{k}{j}{i} + ...
                                        TFA_analysis.recall.NP_AL.([n_group{b} '_' n_syn_g{3-s}]).([n_power{px}]).pre{k}{j}{i}) / 2;
                                    post_BL = (TFA_analysis.recall.P_BL.([n_group{b} '_' n_syn_g{s}]).([n_power{px}]).post{k}{j}{i} + ...
                                        TFA_analysis.recall.NP_BL.([n_group{b} '_' n_syn_g{3-s}]).([n_power{px}]).post{k}{j}{i}) / 2;
                                    post_AL = (TFA_analysis.recall.P_AL.([n_group{b} '_' n_syn_g{s}]).([n_power{px}]).post{k}{j}{i} + ...
                                        TFA_analysis.recall.NP_AL.([n_group{b} '_' n_syn_g{3-s}]).([n_power{px}]).post{k}{j}{i}) / 2;
                                    pre = [pre ((pre_AL - pre_BL)./pre_BL)*100]; post = [post ((post_AL - post_BL)./post_BL)*100]; 
                                end
                                
                                [pre_CI_L, pre_CI_U, pre_CI_M] = bootstrap(pre, n_BS,[]); % pre stim line
                                [post_CI_L, post_CI_U, post_CI_M] = bootstrap(post, n_BS, []); % post stim line
                                
                                subplot(2,7,sp{p} + (b-1)*7);  hold on;
                                fill([x fliplr(x)],[pre_CI_U fliplr(pre_CI_L)],'k', 'linestyle','none');
                                h(1) = plot(x, pre_CI_M,'k', 'linewidth',2);
                                fill([x fliplr(x)],[post_CI_U fliplr(post_CI_L)],'r', 'linestyle','none');alpha(0.2);
                                h(2) = plot(x, post_CI_M,'r','linewidth',2); 
                                y_min = min([-20 min([pre_CI_L post_CI_L])*1.1]);
                                y_max = max([20 max([pre_CI_U post_CI_U])*1.1]);
                                ylim([floor(y_min/10)*10 ceil(y_max/10)*10]);
                                if(p==1); ylabel('BL-AL % power change'); end 
                                if(b==2); xlabel('Epsilon'); end
                                if(b==1 && p==1); legend([h(1) h(2)],{'pre-stim', 'post-stim'}); end
                                title([n_group{b} ' ' n_power{px} ' ' titles{s}]);
                                text(0.02,0.98,[label{b} label_i{p}],'Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14,'Fontweight','bold');
                            end
                        end
                        subplot(2,7,[6 7]); hold on; box on;
                        fill([x fliplr(x)], [W_CI_L fliplr(W_CI_U)], [0.2 0.2 0.5], 'linestyle','none'); alpha(0.3);
                        plot(x, W_M, 'color', [0.2 0.2 0.5], 'linewidth', 1);
                        xlabel('Epsilon'); ylabel('P <-> NP weights');
                        ax = gca; ax.YAxisLocation = 'right'; ylim([0 ceil(max(W_CI_U)/10)*10]);
                        title('Weight Change');
                        text(0.02,0.98,'C','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14,'Fontweight','bold');
                        
                        saveas(f, [vars.test_f '/Simulations/Analyse Variables/' titles{s} '.tiff']); close(f);
                    end
                end
            end
        end
    end
end
end

