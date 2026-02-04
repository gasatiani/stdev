%% 
%% compare_cells.m
close all; warning('off','all');
set(groot, 'defaultTextInterpreter', 'none');
set(groot, 'defaultAxesTickLabelInterpreter', 'none');
set(groot, 'defaultLegendInterpreter', 'none');
addpath("Functions/")
dataDir = '/home/asatiani/Desktop/Thesis/Experiment/optics/';
format("shortE")

min_wl = 430; max_wl = 645;

% =======================
% NORMALIZATION SWITCH
% =======================
normalize = 0;     % 1 = normalize by total spectrum integral, 0 = no normalization

[files, path_selected] = uigetfile('*fluo.mat','Select one or more .fluo.mat files','MultiSelect','on', dataDir);
if isequal(files,0)
    error('No files selected');
end
if ~iscell(files)
    files = {files};
end
num_files = length(files);

% Load references
[FAD385total, FAD405total, lambda] = load_fluo("Reference/flavine_fluo.mat", min_wl, max_wl);
[NADH385total, NADH405total, ~] = load_fluo("Reference/NADH_fluo_4.mat", min_wl, max_wl);
water = '/home/asatiani/Desktop/Thesis/Experiment/optics/cell_cultures_37Degree/PBS_000_fluo.mat';
[water385total, water405total, ~] = load_fluo(water, min_wl, max_wl);

fluorophores_385 = [FAD385total; NADH385total; water385total];
fluorophores_405 = [FAD405total; NADH405total; water405total];

fields = {'NADH_bound','NADH_free','FAD','PbFMN','Lipo','PpIX620','PpIX636'};
num_fluor = length(fields);

means_385 = zeros(num_files, num_fluor); stds_385 = zeros(num_files, num_fluor);
means_405 = zeros(num_files, num_fluor); stds_405 = zeros(num_files, num_fluor);
optical_redox_ratios=zeros(num_files, 2); stds_optical_redox_ratios=zeros(num_files, 2);
file_labels = cell(num_files,1);

% Add label depending on normalization
if normalize
    norm_label = ' for cell cultures at optimal temperature (normalized)';
else
    norm_label = ' for cell cultures at optimal temperature';
end

for i = 1:num_files
    file = files{i};
    [~, fname, ~] = fileparts(file);
    parts = strsplit(fname, '_');
    if contains(fname, "DEAD")
        % fprintf(" _TRUE_ ")
        main_name = strjoin(parts([4 5 6]), ':');
    else
        % fprintf(" _FALSE_ ")
        main_name = strjoin(parts([4 5]), ':');

    end
    
    file_labels{i} = main_name;

    % Fit each spectrum in the file
    [~, ~, fluorophore, ~, ~, ~, ~, fluo385_each, fluo405_each] = ...
        extract_fluo_no_correction_cells(path_selected, file, min_wl, max_wl, fluorophores_385, fluorophores_405);

    Nspec = size(fluo385_each,1);

    % Compute total intensity for normalization
    total_385 = sum(fluo385_each, 2);  
    total_405 = sum(fluo405_each, 2);

    integrals_385 = zeros(Nspec, num_fluor);
    integrals_405 = zeros(Nspec, num_fluor);
   

    for s = 1:Nspec

        % Normalization factors
        norm385 = total_385(s); 
        norm405 = total_405(s);

        if ~normalize
            norm385 = 1;
            norm405 = 1;
        end

        integrals_385(s,:) = [sum(fluorophore.NADH_bound_385_each{s}), ...
                              sum(fluorophore.NADH_free_385_each{s}), ...
                              sum(fluorophore.FAD_385_each{s}), ...
                              sum(fluorophore.PbFMN_385_each{s}), ...
                              sum(fluorophore.Lipo_385_each{s}), ...
                              sum(fluorophore.PpIX620_385_each{s}), ...
                              sum(fluorophore.PpIX636_385_each{s})] / norm385;

        integrals_405(s,:) = [sum(fluorophore.NADH_bound_405_each{s}), ...
                              sum(fluorophore.NADH_free_405_each{s}), ...
                              sum(fluorophore.FAD_405_each{s}), ...
                              sum(fluorophore.PbFMN_405_each{s}), ...
                              sum(fluorophore.Lipo_405_each{s}), ...
                              sum(fluorophore.PpIX620_405_each{s}), ...
                              sum(fluorophore.PpIX636_405_each{s})] / norm405;
    end

    means_385(i,:) = mean(integrals_385,1); stds_385(i,:) = std(integrals_385,0,1);
    means_405(i,:) = mean(integrals_405,1); stds_405(i,:) = std(integrals_405,0,1);
    optical_redox_ratios(i,:)=[means_385(i,3)/(means_385(i,2)+means_385(i,3)),...
                               means_405(i,3)/(means_405(i,2)+means_405(i,3))];


    % use the propagation of uncertainty for redox=FAD/(FAD+NADH)
    % stdev(redox)=|redox/(FAD+NADH)|*sqrt( NADH^2/FAD^2*stdev(FAD)^2+stdev(NADH)^2-2*NADH/FAD*covariance(FAD,NADH) )
    covariance_NADH_FAD={cov(integrals_385(:,3),integrals_385(:,2)),...
                         cov(integrals_405(:,3),integrals_405(:,2))}; %cell array of covariance matrices for 385 and 406 nm laser
  
    % disp(covariance_NADH_FAD{1}(1,2))
    % disp(covariance_NADH_FAD{2}(1,2))
    % 
    % disp(covariance_NADH_FAD)
    stds_optical_redox_ratios(i,:)=[abs(optical_redox_ratios(i,1)/((means_385(i,2)+means_385(i,3))))*...
                                    sqrt((means_385(i,2)^2)/(means_385(i,3)^2)*stds_385(i,3)^2+...
                                    stds_385(i,2)^2-...
                                    2*means_385(i,2)/means_385(i,3)*covariance_NADH_FAD{1}(1,2)),...
                                    abs(optical_redox_ratios(i,2)/((means_405(i,2)+means_405(i,3))))*...
                                    sqrt((means_405(i,2)^2)/(means_405(i,3)^2)*stds_405(i,3)^2+...
                                    stds_405(i,2)^2-...
                                    2*means_405(i,2)/means_405(i,3)*covariance_NADH_FAD{2}(1,2))];
end


%% Select directory for saving figures
save_dir = uigetdir('/home/asatiani/Desktop/Thesis/Figures', 'Select folder to save the generated figures');
if save_dir == 0
    error('No folder selected. Aborting figure saving.');
end

%% ----------------------------------------------
%  AUTOMATIC COMPARISON PLOTS (385 vs 405)
% ----------------------------------------------

% Names of fluorophores in order (corresponding to integrals)
fluor_names = { ...
    'NADH_bound', ...
    'NADH_free', ...
    'FAD', ...
    'PbFMN', ...
    'Lipopigments', ...
    'PpIX620', ...
    'PpIX636'};

num_fluor = numel(fluor_names);

for f = 1:num_fluor
    % Data for this fluorophore
    data_f = [means_385(:,f), means_405(:,f)];
    stds_f = [stds_385(:,f), stds_405(:,f)];
    
    % Create figure
    figure('Position',[100 100 1000 700]); hold on;
    b = bar(data_f);  

    for k = 1:2
        errorbar(b(k).XEndPoints, data_f(:,k), stds_f(:,k), ...
                 'k', 'linestyle','none');
    end
    
    set(gca,'XTick',1:num_files, ...
            'XTickLabel', file_labels, ...
            'XTickLabelRotation', 45);
    

    ylabel('Fluorescence intensity (a.u.)', 'Interpreter','none'); 
    title([fluor_names{f} ' comparison 385 vs 405' norm_label], ...
          'Interpreter','none');
    
    legend({[fluor_names{f} ' 385'], [fluor_names{f} ' 405']}, ...
           'Interpreter','none', 'Location','bestoutside');
    ax = gca;
    annotation('doublearrow', [ax.Position(1),0.68], [0.91,0.91],Color='r');
    text(0.35,0.95,"Kept at optimal temperature", Units="normalized", HorizontalAlignment='center', Color='r')

    x_norm = ax.Position(1) + (ax.XTick(end)-ax.XLim(1))/diff(ax.XLim)*ax.Position(3);
    y_norm = ax.Position(2);  % at bottom, pointing at XTick
    annotation('textarrow', [x_norm+0.15, x_norm-0.01], [y_norm-0.06, y_norm-0.06], ...
               'String', "Killed by heat", 'Color', 'r');
    grid on;
    
    % Auto filename
    filename = [fluor_names{f} '_385_vs_405' norm_label '.png'];
    filename = strrep(filename, ' ', '_');  % clean names
    
    % Save
    saveas(gcf, fullfile(save_dir, filename));

    pause(0.01);
    
end


%% ----------------------------------------------
%  ADDITIONAL PLOTS: All fluorophores 385 + 405 
% ----------------------------------------------

% 385 nm full plot
figure('Position',[100 100 1000 700]); hold on;
b = bar(means_385);
colors = lines(num_fluor);
for k = 1:num_fluor
    b(k).FaceColor = colors(k,:);
    errorbar(b(k).XEndPoints, means_385(:,k), stds_385(:,k), ...
             'k','linestyle','none');
end
set(gca,'XTick',1:num_files,'XTickLabel',file_labels,'XTickLabelRotation',45);

ylabel('Fluorescence intensity (a.u.)');
title(['All fluorophores - 385 nm' norm_label]);
legend(fluor_names, 'Location','bestoutside');
ax = gca;
x_norm = ax.Position(1) + (ax.XTick(end)-ax.XLim(1))/diff(ax.XLim)*ax.Position(3);
y_norm = ax.Position(2);  % at bottom, pointing at XTick
annotation('textarrow', [x_norm+0.15, x_norm-0.01], [y_norm-0.06, y_norm-0.06], ...
           'String', "Killed by heat", 'Color', 'r');
annotation('doublearrow', [ax.Position(1),0.68], [0.91,0.91],Color='r');
text(0.35,0.95,"Kept at optimal temperature", Units="normalized", HorizontalAlignment='center', Color='r');
grid on;
filename=['AllFluorophores_385' norm_label '.png'];
filename = strrep(filename, ' ', '_');  % clean names
saveas(gcf, fullfile(save_dir, filename));

% 405 nm full plot
figure('Position',[100 100 1000 700]); hold on;
b = bar(means_405);
for k = 1:num_fluor
    b(k).FaceColor = colors(k,:);
    errorbar(b(k).XEndPoints, means_405(:,k), stds_405(:,k), ...
             'k','linestyle','none');
end
set(gca,'XTick',1:num_files,'XTickLabel',file_labels,'XTickLabelRotation',45);

    
ylabel('Fluorescence intensity (a.u.)');
title(['All fluorophores - 405 nm' norm_label]);
legend(fluor_names, 'Location','bestoutside');
ax = gca;


x_norm = ax.Position(1) + (ax.XTick(end)-ax.XLim(1))/diff(ax.XLim)*ax.Position(3);
y_norm = ax.Position(2);  % at bottom, pointing at XTick
annotation('textarrow', [x_norm+0.15, x_norm-0.01], [y_norm-0.06, y_norm-0.06], ...
           'String', "Killed by heat", 'Color', 'r');
annotation('doublearrow', [ax.Position(1),0.68], [0.91,0.91],Color='r');
text(0.35,0.95,"Kept at optimal temperature", Units="normalized", HorizontalAlignment='center', Color='r');

grid on;
filename=['AllFluorophores_405' norm_label '.png'];
filename = strrep(filename, ' ', '_');  % clean names
saveas(gcf, fullfile(save_dir, filename));


%% ----------------------------------------------
%  ADDITIONAL PLOTS: Optical redox ratios
% ----------------------------------------------

% 385 nm optical redox ratio
figure('Position',[100 100 1000 700]); hold on;
b = bar(optical_redox_ratios(:,1));
colors = lines(1);
b(1).FaceColor = colors(1,:);
errorbar(b(1).XEndPoints, optical_redox_ratios(:,1), stds_optical_redox_ratios(:,1), ...
         'k','linestyle','none');
set(gca,'XTick',1:num_files,'XTickLabel',file_labels,'XTickLabelRotation',45);


ylabel('Redox ratio of FAD/(NADH+FAD)');
title(['Optical redox ratio - 385 nm' norm_label]);
ax = gca;
x_norm = ax.Position(1) + (ax.XTick(end)-ax.XLim(1))/diff(ax.XLim)*ax.Position(3);
y_norm = ax.Position(2);  % at bottom, pointing at XTick
annotation('textarrow', [x_norm+0.15, x_norm-0.01], [y_norm-0.06, y_norm-0.06], ...
           'String', "Killed by heat", 'Color', 'r');
annotation('doublearrow', [ax.Position(1),0.75], [0.91,0.91],Color='r');
text(0.35,0.95,"Kept at optimal temperature", Units="normalized", HorizontalAlignment='center', Color='r')
grid on;
filename=['Optical_redox_ratio_385' norm_label '.png'];
filename = strrep(filename, ' ', '_');  % clean names
saveas(gcf, fullfile(save_dir, filename));

% 405 nm full plot
figure('Position',[100 100 1000 700]); hold on;
b = bar(optical_redox_ratios(:,2));
b(1).FaceColor = colors(1,:);
errorbar(b(1).XEndPoints, optical_redox_ratios(:,2), stds_optical_redox_ratios(:,2), ...
         'k','linestyle','none');
set(gca,'XTick',1:num_files,'XTickLabel',file_labels,'XTickLabelRotation',45);

    
ylabel('Redox ratio of FAD/(NADH+FAD)');
title(['Optical redox ratio - 405 nm' norm_label]);
ax = gca;
x_norm = ax.Position(1) + (ax.XTick(end)-ax.XLim(1))/diff(ax.XLim)*ax.Position(3);
y_norm = ax.Position(2);  % at bottom, pointing at XTick


annotation('doublearrow', [ax.Position(1),0.75], [0.91,0.91],Color='r');
annotation('textarrow', [x_norm+0.15, x_norm-0.01], [y_norm-0.06, y_norm-0.06], ...
           'String', "Killed by heat", 'Color', 'r');
text(0.35,0.95,"Kept at optimal temperature", Units="normalized", HorizontalAlignment='center', Color='r')

grid on;
filename=['Optical_redox_ratio_405' norm_label '.png'];
filename = strrep(filename, ' ', '_');  % clean names
saveas(gcf, fullfile(save_dir, filename));
