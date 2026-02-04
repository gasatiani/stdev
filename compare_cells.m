%% ==============================================================
% PIPELINE:
%   1. Configuration
%   2. File selection
%   3. Reference loading
%   4. Spectral processing
%   5. Optical redox computation
%   6. Plot generation
%% ==============================================================

close all; warning('off','all');

%% --------------------------------------------------------------
% GLOBAL PLOT SETTINGS
%% --------------------------------------------------------------
set(groot,'defaultTextInterpreter','none');
set(groot,'defaultAxesTickLabelInterpreter','none');
set(groot,'defaultLegendInterpreter','none');

addpath("Functions/")
format("shortE")

%% --------------------------------------------------------------
% CONFIGURATION
%% --------------------------------------------------------------
dataDir   = '/home/asatiani/Desktop/Thesis/Experiment/optics/';
save_root = '/home/asatiani/Desktop/Thesis/Figures';

min_wl = 430;
max_wl = 645;

normalize = 1;   % 1 = normalize by spectrum integral

% Fluorophore names used for indexing + plotting
fluor_names = { ...
    'NADH_bound', ...
    'NADH_free', ...
    'FAD', ...
    'PbFMN', ...
    'Lipo', ...
    'PpIX620', ...
    'PpIX636'};
num_fluor = numel(fluor_names);

if normalize
    norm_label = ' for yeast cells at optimal temperature (normalized)';
else
    norm_label = ' for yeast cells at optimal temperature';
end

% Start parallel workers if not already running
if isempty(gcp('nocreate'))
    parpool;
end


%% --------------------------------------------------------------
% FILE SELECTION
%% --------------------------------------------------------------
[files,path_selected] = uigetfile('*fluo.mat',...
    'Select one or more .fluo.mat files',...
    'MultiSelect','on',dataDir);

if isequal(files,0)
    error('No files selected');
end
if ~iscell(files)
    files = {files};
end
num_files = numel(files);

%% --------------------------------------------------------------
% LOAD REFERENCE SPECTRA
%% --------------------------------------------------------------
[FAD385,FAD405,lambda] = load_fluo("Reference/flavine_fluo.mat",min_wl,max_wl);
[NADH385,NADH405,~]    = load_fluo("Reference/NADH_fluo_4.mat",min_wl,max_wl);

water_file = '/home/asatiani/Desktop/Thesis/Experiment/optics/cell_cultures_37Degree/PBS_000_fluo.mat';
[water385,water405,~] = load_fluo(water_file,min_wl,max_wl);

fluorophores_385 = [FAD385; NADH385; water385];
fluorophores_405 = [FAD405; NADH405; water405];

%% --------------------------------------------------------------
% PREALLOCATION
%% --------------------------------------------------------------
means_385 = zeros(num_files,num_fluor);
stds_385  = zeros(num_files,num_fluor);
means_405 = zeros(num_files,num_fluor);
stds_405  = zeros(num_files,num_fluor);

optical_redox = zeros(num_files,2);
std_redox     = zeros(num_files,2);

file_labels = cell(num_files,1);

%% --------------------------------------------------------------
% MAIN PROCESSING LOOP
%% --------------------------------------------------------------
for i = 1:num_files

    file = files{i};

    % -------- Extract label from filename --------
    [~,fname,~] = fileparts(file);
    parts = strsplit(fname,'_');

    if contains(fname,"DEAD")
        file_labels{i} = strjoin(parts([1 2]),':');
    else
        file_labels{i} = strjoin(parts([1 2]),':');
    end

    % -------- Decompose fluorescence spectra --------
    [~,~,fluorophore,~,~,~,~,fluo385_each,fluo405_each] = ...
        extract_fluo_no_correction_cells( ...
        path_selected,file,min_wl,max_wl,...
        fluorophores_385,fluorophores_405);

    Nspec = size(fluo385_each,1);

    total_385 = sum(fluo385_each,2);
    total_405 = sum(fluo405_each,2);

    integrals_385 = zeros(Nspec,num_fluor);
    integrals_405 = zeros(Nspec,num_fluor);

    % -------- Compute per-spectrum fluorophore integrals --------
    parfor s = 1:Nspec
    
        norm385 = total_385(s);
        norm405 = total_405(s);
        if ~normalize
            norm385 = 1;
            norm405 = 1;
        end
    
        local_385 = zeros(1,num_fluor);
        local_405 = zeros(1,num_fluor);
    
        for f = 1:num_fluor
            fname_f = fluor_names{f};
    
            local_385(f) = ...
                sum(fluorophore.([fname_f '_385_each']){s}) ./ norm385;
    
            local_405(f) = ...
                sum(fluorophore.([fname_f '_405_each']){s}) ./ norm405;
        end
    
        integrals_385(s,:) = local_385;
        integrals_405(s,:) = local_405;
    end


    % -------- Compute statistics across spectra --------
    means_385(i,:) = mean(integrals_385,1);
    stds_385(i,:)  = std(integrals_385,0,1);
    means_405(i,:) = mean(integrals_405,1);
    stds_405(i,:)  = std(integrals_405,0,1);

    % -------- Optical redox ratio --------
    optical_redox(i,:) = [ ...
        means_385(i,3)/(means_385(i,2)+means_385(i,3)), ...
        means_405(i,3)/(means_405(i,2)+means_405(i,3))];

    cov385 = cov(integrals_385(:,3),integrals_385(:,2));
    cov405 = cov(integrals_405(:,3),integrals_405(:,2));

    std_redox(i,:) = [ ...
        compute_redox_std(optical_redox(i,1),means_385(i,2),means_385(i,3),stds_385(i,2),stds_385(i,3),cov385(1,2)), ...
        compute_redox_std(optical_redox(i,2),means_405(i,2),means_405(i,3),stds_405(i,2),stds_405(i,3),cov405(1,2))];
end

%% --------------------------------------------------------------
% SAVE DIRECTORY
%% --------------------------------------------------------------
save_dir = uigetdir(save_root,'Select folder to save the generated figures');
if save_dir==0
    error('No folder selected');
end

%% --------------------------------------------------------------
% PLOTTING
%% --------------------------------------------------------------
colors = lines(num_fluor);

for f = 1:num_fluor
    plot_compare(means_385(:,f),means_405(:,f),...
        stds_385(:,f),stds_405(:,f),...
        file_labels,fluor_names{f},norm_label,save_dir);    
end 

plot_full(means_385,stds_385,file_labels,fluor_names,colors,'385 nm',norm_label,save_dir);
plot_full(means_405,stds_405,file_labels,fluor_names,colors,'405 nm',norm_label,save_dir);

plot_redox(optical_redox(:,1),std_redox(:,1),file_labels,'385 nm',norm_label,save_dir);
plot_redox(optical_redox(:,2),std_redox(:,2),file_labels,'405 nm',norm_label,save_dir);

%% ==============================================================
% LOCAL FUNCTIONS
%% ==============================================================

function s = compute_redox_std(redox,NADH,FAD,stdNADH,stdFAD,covNF)
% Propagation of uncertainty for:
%   redox = FAD / (FAD + NADH)
s = abs(redox/(NADH+FAD))* ...
    sqrt((NADH^2)/(FAD^2)*stdFAD^2 + stdNADH^2 - 2*NADH/FAD*covNF);
end

function add_annotations(ax)
% Adds experimental condition annotations
x_norm = ax.Position(1) + ...
    (ax.XTick(end)-ax.XLim(1))/diff(ax.XLim)*ax.Position(3);
y_norm = ax.Position(2);

annotation('doublearrow',[ax.Position(1)+0.02,x_norm-0.06],[0.91,0.91],'Color','r');
annotation('textarrow',[x_norm+0.06,x_norm-0.01],...
    [y_norm-0.06,y_norm-0.06],'String',"Killed by heat",'Color','r');

text(0.35,0.95,"Kept at optimal temperature",...
    'Units',"normalized",'HorizontalAlignment','center','Color','r')
end

function plot_compare(m385,m405,s385,s405,labels,name,norm_label,save_dir)
figure('Position',[100 100 1000 700]); hold on;
data=[m385,m405]; stds=[s385,s405];
b=bar(data);

for k=1:2
    errorbar(b(k).XEndPoints,data(:,k),stds(:,k),'k','linestyle','none');
end

set(gca,'XTick',1:numel(labels),'XTickLabel',labels,'XTickLabelRotation',45);
ylabel('Fluorescence intensity (a.u.)');
title([name ' comparison 385 vs 405' norm_label]);
legend({[name ' 385'],[name ' 405']},'Location','bestoutside');

% add_annotations(gca); 
grid on;

fname=strrep([name '_385_vs_405' norm_label '.png'],' ','_');
saveas(gcf,fullfile(save_dir,fname));
end

function plot_full(means,stds,labels,names,colors,laser,norm_label,save_dir)
figure('Position',[100 100 1000 700]); hold on;
b=bar(means);

for k=1:numel(names)
    b(k).FaceColor=colors(k,:);
    errorbar(b(k).XEndPoints,means(:,k),stds(:,k),'k','linestyle','none');
end

set(gca,'XTick',1:numel(labels),'XTickLabel',labels,'XTickLabelRotation',45);
ylabel('Fluorescence intensity (a.u.)');
title(['All fluorophores - ' laser norm_label]);
legend(names,'Location','bestoutside');

% add_annotations(gca); 
grid on;

fname=strrep(['AllFluorophores_' strrep(laser,' ','') norm_label '.png'],' ','_');
saveas(gcf,fullfile(save_dir,fname));
end

function plot_redox(data,stds,labels,laser,norm_label,save_dir)
figure('Position',[100 100 1000 700]); hold on;
b=bar(data);

errorbar(b.XEndPoints,data,stds,'k','linestyle','none');

set(gca,'XTick',1:numel(labels),'XTickLabel',labels,'XTickLabelRotation',45);
ylabel('Redox ratio of FAD/(NADH+FAD)');
title(['Optical redox ratio - ' laser norm_label]);

% add_annotations(gca); 
grid on;

fname=strrep(['Optical_redox_ratio_' strrep(laser,' ','') norm_label '.png'],' ','_');
saveas(gcf,fullfile(save_dir,fname));
end