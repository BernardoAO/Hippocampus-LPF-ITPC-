%%%% TD speed Analysis
%%% Bernardo AO, adapted from ITPC_analysis.m
% Data needed: sessionInfo.mat, CSC(n).ncs
% Script dependencies: ITPC_analysis.m,  plotTrialPerCell, trialInfo2simple, readCRTsd
% 
addpath(genpath('W:\LABanalysis\SilviaProjectCode\AnalysisPerTrial\RunAnalysis'))
addpath("W:\Lorena\Analysis_scripts\Bernardo_code\")

disp('Starting ITPC Analysis')
clc

%% Define which animal / days to run 

load(fullfile('W:\Lorena\Analysis_scripts\DataOrganization', 'sessionInfo.mat'));
All_sessInfo = sessInfo; clear sessInfo
animal_numbers = unique([All_sessInfo.animal]);
all_rois = {'return';'delay';'stem';'choice';'reward'};
genotypes = ["Control", "CA1-APP"];
colors = [0.14,0.85,0.71; 0.85,0.14,0.28];
sp = "W:\Lorena\Analysis_scripts\Bernardo_code\plots";

% Parameters 
win = 640; % window to smooth the ratio around win / fs = +- 20 ms

% Outputs
problematicSessions = {};
results = struct([]);

wb = waitbar(0,'Running animal loop 0%');
a = 1;
%animal = 1433; sessInfo = All_sessInfo(1); block_c = sessInfo(i).sessDirs 
for animal = animal_numbers
    sessions = {All_sessInfo.animal};
    s_indices = find(cellfun(@(sessions) ~isempty(sessions) && ...
        sessions==animal, sessions));
    
    lfps = {}; 
    pideals = {};
    
    s = 0;
    for iii = s_indices

        sessInfo = All_sessInfo(iii);        
        ptempl = parsingtemplate('fig8mouse:rdscr');
                
        i = 1;
        try
            % Liniarize path
            [trialInfo, parsingInfo, pathData, pathDataIdeal, ...
                pathDataLin] = plotTrialPerCell.loadInfo(sessInfo, i);

            s = s + 1;
            for block_c = {'d10'}%sessInfo(i).sessDirs
                block = block_c{1};
                blockDir = fullfile(sessInfo(i).mainDir, block);
                LFP_dir = fullfile(blockDir,'LFP'); if ...
                    ~exist(LFP_dir,'dir'), mkdir(LFP_dir),end
    
                tinfo = trialInfo.(block);
                pideal = pathDataIdeal.(block);
    
                
                %% Collect LFPs
                channelInLayer = sprintf('CSC%d.ncs', sessInfo.cellLayerChann); % Picks the channel that is in the layer
                lfp = readCRTsd(fullfile(blockDir, channelInLayer));
    
                lfps{s} = lfp;
                pideals{s} = pideal;

            end
        catch
            problematicSessions = [problematicSessions; sessInfo(i).mainDir];
            continue;
        end
        clc
    end
    
    %% Run TD_ratio
    fprintf('\nRunning TD ratio\n');
    opt_td = struct(...
                'fs', EEG.Fs(lfp), ...
                'n_cycles', 10, ...
                'color', colors(sessInfo.genotype,:), ...
                'fig_name', int2str(animal), ...
                'win', win, ...
                'downsample', 2000, ...
                'save_path', "W:\Lorena\Analysis_scripts\Bernardo_code\plots", ...
                'ext', ".png");
 
    %try
        [vel_grid, ratio_grid] = theta_delta_ratio(pideals, lfps, opt_td);
        results(a).vel_grid = vel_grid;
        results(a).ratio_grid = ratio_grid;
    %catch
    %    problematicSessions = [problematicSessions; animal];
    %end

    %% Save
    results(a).animal = animal;
    results(a).genotype = sessInfo.genotype;
    
    r = a/length(animal_numbers);
    waitbar(r,wb,"Running " + num2str(r*100,2) + "%");
    a = a + 1;
end
clc
delete(wb);

%% Results

plot_mean_grid(results,"vel_grid", "velocity", opt_td)
plot_mean_grid(results,"ratio_grid", "ratio", opt_td)

%% Groups

for g = 1:2
    r = results([results.genotype] == g);
    plot_mean_grid(r,"vel_grid", "velocity " + genotypes(g), opt_td)
    plot_mean_grid(r,"ratio_grid", "ratio " + genotypes(g), opt_td)
end

function plot_mean_grid(results, gname,name, opt)
    
    x_edges = -23:23;
    y_edges = -35:35;
    x_centers = (x_edges(1:end-1) + x_edges(2:end)) / 2;
    y_centers = (y_edges(1:end-1) + y_edges(2:end)) / 2;
    
    g_shape = size(results(1).(gname));
    n_animals = numel(results);
    all_grid = zeros(g_shape(1), g_shape(2), n_animals);

    for a = 1:n_animals
        if name == "theta power" || name == "delta power"
            vel_grid = results(a).(gname);
            m = mean(vel_grid,"all","omitmissing");
            s = std(vel_grid, [],"all","omitmissing");
            all_grid(:,:,a) = (vel_grid - m) / s;
        else
            all_grid(:,:,a) = results(a).(gname);
        end
    end
    mean_grid = mean(all_grid,3);

    % Plot
    figure;
    imagesc(x_centers, y_centers, mean_grid);
    axis xy; % flip y-axis so low values are at the bottom
    t_name = " Average " + name + " heatmap";
    title(t_name);
    axis off
    colormap(flipud(hot))
    a = colorbar;
    a.Label.String = name;
    saveas(gcf, fullfile(opt.save_path, t_name + "all.png"));
end

