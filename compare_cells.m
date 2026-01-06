%% compare_cells.m
clear all; close all; warning('off','all');

dataDir = '/home/asatiani/Desktop/Thesis/Experiment/optics/';
min_wl = 430; max_wl = 645;

% =======================
% NORMALIZATION SWITCH
% =======================
normalize = 1;     % 1 = normalize by total spectrum integral, 0 = no normalization

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
file_labels = cell(num_files,1);

% Add label depending on normalization
if normalize
    norm_label = ' (normalized)';
else
    norm_label = '';
end

for i = 1:num_files
    file = files{i};
    [~, fname, ~] = fileparts(file);
    parts = strsplit(fname, '_');
    main_name = strjoin(parts([1 2 3]), ':');
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
end

% %% Plot all fluorophores - 385 nm
% figure('Position',[100 100 1000 700]); hold on;
% b = bar(means_385);
% colors = lines(num_fluor);
% for k = 1:num_fluor
%     b(k).FaceColor = colors(k,:);
%     errorbar(b(k).XEndPoints, means_385(:,k), stds_385(:,k), 'k','linestyle','none');
% end
% set(gca,'XTick',1:num_files,'XTickLabel',file_labels,'XTickLabelRotation',45);
% ylabel('Fluorescence intensity (a.u.)','Interpreter','none'); 
% title(['All fluorophores - Excitation 385 nm' norm_label],'Interpreter','none');
% legend(fields,'Location','bestoutside','Interpreter','none'); grid on;
% 
% %% Plot all fluorophores - 405 nm
% figure('Position',[100 100 1000 700]); hold on;
% b = bar(means_405);
% for k = 1:num_fluor
%     b(k).FaceColor = colors(k,:);
%     errorbar(b(k).XEndPoints, means_405(:,k), stds_405(:,k), 'k','linestyle','none');
% end
% set(gca,'XTick',1:num_files,'XTickLabel',file_labels,'XTickLabelRotation',45);
% ylabel('Fluorescence intensity (a.u.)','Interpreter','none'); 
% title(['All fluorophores - Excitation 405 nm' norm_label],'Interpreter','none');
% legend(fields,'Location','bestoutside','Interpreter','none'); grid on;
% 
% 
% %% NADH bound comparison (both lasers)
% figure('Position',[100 100 1000 700]); hold on;
% data_nadh_bound = [means_385(:,1), means_405(:,1)];
% stds_nadh_bound = [stds_385(:,1), stds_405(:,1)];
% b = bar(data_nadh_bound);
% for k = 1:2
%     errorbar(b(k).XEndPoints, data_nadh_bound(:,k), stds_nadh_bound(:,k), 'k','linestyle','none');
% end
% set(gca,'XTick',1:num_files,'XTickLabel',file_labels,'XTickLabelRotation',45);
% ylabel('Fluorescence intensity (a.u.)','Interpreter','none'); 
% title(['NADH bound comparison 385 vs 405' norm_label],'Interpreter','none'); 
% legend({'NADH bound 385','NADH bound 405'},'Interpreter','none'); grid on;
% 
% 
% %% NADH comparison (both lasers)
% figure('Position',[100 100 1000 700]); hold on;
% data_nadh = [means_385(:,2), means_405(:,2)];
% stds_nadh = [stds_385(:,2), stds_405(:,2)];
% b = bar(data_nadh);
% for k = 1:2
%     errorbar(b(k).XEndPoints, data_nadh(:,k), stds_nadh(:,k), 'k','linestyle','none');
% end
% set(gca,'XTick',1:num_files,'XTickLabel',file_labels,'XTickLabelRotation',45);
% ylabel('Fluorescence intensity (a.u.)','Interpreter','none'); 
% title(['NADH free comparison 385 vs 405' norm_label],'Interpreter','none'); 
% legend({'NADH 385','NADH 405'},'Interpreter','none'); grid on;
% 
% %% FAD comparison (both lasers)
% figure('Position',[100 100 1000 700]); hold on;
% data_fad = [means_385(:,3), means_405(:,3)];
% stds_fad = [stds_385(:,3), stds_405(:,3)];
% b = bar(data_fad);
% for k = 1:2
%     errorbar(b(k).XEndPoints, data_fad(:,k), stds_fad(:,k), 'k','linestyle','none');
% end
% set(gca,'XTick',1:num_files,'XTickLabel',file_labels,'XTickLabelRotation',45);
% ylabel('Fluorescence intensity (a.u.)','Interpreter','none'); 
% title(['FAD comparison 385 vs 405' norm_label],'Interpreter','none'); 
% legend({'FAD 385','FAD 405'},'Interpreter','none'); grid on;
% 
% 
% %% Lipopigment comparison (both lasers)
% 
% figure('Position',[100 100 1000 700]); hold on;
% data_lipo = [means_385(:,5), means_405(:,5)];
% stds_lipo = [stds_385(:,5), stds_405(:,5)];
% b = bar(data_lipo);
% for k = 1:2
%     errorbar(b(k).XEndPoints, data_lipo(:,k), stds_lipo(:,k), 'k','linestyle','none');
% end
% set(gca,'XTick',1:num_files,'XTickLabel',file_labels,'XTickLabelRotation',45);
% ylabel('Fluorescence intensity (a.u.)','Interpreter','none'); 
% title(['Lipopigments comparison 385 vs 405' norm_label],'Interpreter','none'); 
% legend({'Lipopigments 385','Lipopigments 405'},'Interpreter','none'); grid on;
% 
% %% PbFMN (bound FMN) comparison (both lasers)
% figure('Position',[100 100 1000 700]); hold on;
% data_fmn = [means_385(:,4), means_405(:,4)];
% stds_fmn = [stds_385(:,4), stds_405(:,4)];
% b = bar(data_fmn);
% for k = 1:2
%     errorbar(b(k).XEndPoints, data_fmn(:,k), stds_fmn(:,k), 'k','linestyle','none');
% end
% set(gca,'XTick',1:num_files,'XTickLabel',file_labels,'XTickLabelRotation',45);
% ylabel('Fluorescence intensity (a.u.)','Interpreter','none'); 
% title(['PbFMN comparison 385 vs 405' norm_label],'Interpreter','none'); 
% legend({'PbFMN 385','PbFMN 405'},'Interpreter','none'); grid on;
% 
% 
% %% PpIX 620 comparison (both lasers)
% figure('Position',[100 100 1000 700]); hold on;
% data_ppix620 = [means_385(:,6), means_405(:,6)];
% stds_ppix620 = [stds_385(:,6), stds_405(:,6)];
% b = bar(data_ppix620);
% for k = 1:2
%     errorbar(b(k).XEndPoints, data_ppix620(:,k), stds_ppix620(:,k), 'k','linestyle','none');
% end
% set(gca,'XTick',1:num_files,'XTickLabel',file_labels,'XTickLabelRotation',45);
% ylabel('Fluorescence intensity (a.u.)','Interpreter','none'); 
% title(['PpIX 620 comparison 385 vs 405' norm_label],'Interpreter','none'); 
% legend({'PpIX 620 385','PpIX 620 405'},'Interpreter','none'); grid on;
% 
% 
% 
% 
% 
% %% PpIX 636 comparison (both lasers)
% figure('Position',[100 100 1000 700]); hold on;
% data_ppix636 = [means_385(:,7), means_405(:,7)];
% stds_ppix636 = [stds_385(:,7), stds_405(:,7)];
% b = bar(data_ppix636);
% for k = 1:2
%     errorbar(b(k).XEndPoints, data_ppix636(:,k), stds_ppix636(:,k), 'k','linestyle','none');
% end
% set(gca,'XTick',1:num_files,'XTickLabel',file_labels,'XTickLabelRotation',45);
% ylabel('Fluorescence intensity (a.u.)','Interpreter','none'); 
% title(['PpIX 636 comparison 385 vs 405' norm_label],'Interpreter','none'); 
% legend({'PpIX 636 385','PpIX 636 405'},'Interpreter','none'); grid on;
% 
%% Select directory for saving figures
save_dir = uigetdir(pwd, 'Select folder to save the generated figures');
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
    
    grid on;
    
    % Auto filename
    filename = [fluor_names{f} '_385_vs_405' norm_label '.png'];
    filename = strrep(filename, ' ', '_');  % clean names
    
    % Save
    saveas(gcf, fullfile(save_dir, filename));
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
grid on;
saveas(gcf, fullfile(save_dir, ['AllFluorophores_385' norm_label '.png']));

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
grid on;
saveas(gcf, fullfile(save_dir, ['AllFluorophores_405' norm_label '.png']));


