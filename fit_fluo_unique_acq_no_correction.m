%% Fit and plot one acquisition (stacked 385 + 405)
clear all; close all; warning('off','all');

dataDir = '/home/asatiani/Desktop/Thesis/Experiment/optics/';
min_wl = 430; max_wl = 645;

[file, path_selected] = uigetfile('*fluo.mat','Select a .fluo.mat file', dataDir);
if isequal(file,0)
    error('No file selected');
end

% Load reference fluorophores
[FAD385total, FAD405total, lambda] = load_fluo("Reference/flavine_fluo.mat", min_wl, max_wl);
[NADH385total, NADH405total, ~] = load_fluo("Reference/NADH_fluo_4.mat", min_wl, max_wl);
water = '/home/asatiani/Desktop/Thesis/Experiment/optics/cell_cultures_37Degree/PBS_000_fluo.mat';
[water385total, water405total, ~] = load_fluo(water, min_wl, max_wl);

fluorophores_385 = [FAD385total; NADH385total; water385total];
fluorophores_405 = [FAD405total; NADH405total; water405total];

% Run extraction (fits each spectrum)
[S385total, S405total, fluorophore, lambda, ~, ~, ~, fluo385_each, fluo405_each] = ...
    extract_fluo_no_correction_cells(path_selected, file, min_wl, max_wl, fluorophores_385, fluorophores_405);

fields = {'NADH_bound','NADH_free','FAD','PbFMN','Lipo','PpIX620','PpIX636'};
colors = lines(length(fields));
Nspec = size(fluo385_each,1);

for s = 1:Nspec
    figure('Name',sprintf('%s - Spectrum %d', file, s),'NumberTitle','off');
    
    % Stacked plot
    subplot(2,1,1); % 385 nm
    hold on
    plot(lambda, fluo385_each(s,:), 'k', 'LineWidth', 1.2, 'DisplayName', 'Experimental');
    plot(lambda, fluorophore.fit_385_each{s}, 'r', 'LineWidth', 1.4, 'DisplayName', 'Total fit');
    for f = 1:length(fields)
        compField = sprintf('%s_385_each', fields{f});
        plot(lambda, fluorophore.(compField){s}, '--', 'Color', colors(f,:), 'LineWidth', 1.0, 'DisplayName', fields{f});
    end
    xlabel('Wavelength (nm)','Interpreter','none'); 
    ylabel('Fluorescence intensity (a.u.)','Interpreter','none');
    title(sprintf('laser 385 - %s - Spectrum %d', file, s),'Interpreter','none');
    grid on; legend('show','Interpreter','none'); hold off

    subplot(2,1,2); % 405 nm
    hold on
    plot(lambda, fluo405_each(s,:), 'k', 'LineWidth', 1.2, 'DisplayName', 'Experimental');
    plot(lambda, fluorophore.fit_405_each{s}, 'r', 'LineWidth', 1.4, 'DisplayName', 'Total fit');
    for f = 1:length(fields)
        compField = sprintf('%s_405_each', fields{f});
        plot(lambda, fluorophore.(compField){s}, '--', 'Color', colors(f,:), 'LineWidth', 1.0, 'DisplayName', fields{f});
    end
    xlabel('Wavelength (nm)','Interpreter','none'); 
    ylabel('Fluorescence intensity (a.u.)','Interpreter','none');
    title(sprintf('laser 405 - %s - Spectrum %d', file, s),'Interpreter','none');
    grid on; legend('show','Interpreter','none'); hold off
end

%% Compute and plot mean spectra with contributions

% Preallocate mean contribution matrices
mean_contrib_385 = zeros(length(fields), length(lambda));
mean_contrib_405 = zeros(length(fields), length(lambda));

for f = 1:length(fields)
    compField385 = sprintf('%s_385_each', fields{f});
    compField405 = sprintf('%s_405_each', fields{f});
    
    temp385 = zeros(Nspec, length(lambda));
    temp405 = zeros(Nspec, length(lambda));
    
    for s = 1:Nspec
        % Interpolate each spectrum to common lambda grid
        temp385(s,:) = interp1(lambda, fluorophore.(compField385){s}, lambda, 'linear', 0);
        temp405(s,:) = interp1(lambda, fluorophore.(compField405){s}, lambda, 'linear', 0);
    end
    
    mean_contrib_385(f,:) = mean(temp385, 1);
    mean_contrib_405(f,:) = mean(temp405, 1);
end

% Mean of total fit and experimental spectra
mean_fluo385 = mean(fluo385_each, 1);
mean_fluo405 = mean(fluo405_each, 1);

% Interpolate and compute mean of total fits
temp_fit385 = zeros(Nspec, length(lambda));
temp_fit405 = zeros(Nspec, length(lambda));

for s = 1:Nspec
    temp_fit385(s,:) = interp1(lambda, fluorophore.fit_385_each{s}, lambda, 'linear', 0);
    temp_fit405(s,:) = interp1(lambda, fluorophore.fit_405_each{s}, lambda, 'linear', 0);
end

mean_fit385 = mean(temp_fit385, 1);
mean_fit405 = mean(temp_fit405, 1);



%% Plot mean spectra (385 + 405) in a single figure with proper labels
figure('Name','Mean Spectra 385 + 405','NumberTitle','off');

% 385 nm
subplot(2,1,1); hold on
plot(lambda, mean_fluo385, 'k', 'LineWidth', 1.5, 'DisplayName', 'Experimental Mean');
plot(lambda, mean_fit385, 'r', 'LineWidth', 1.5, 'DisplayName', 'Total Fit Mean');
for f = 1:length(fields)
    plot(lambda, mean_contrib_385(f,:), '--', 'Color', colors(f,:), 'LineWidth', 1.2, 'DisplayName', fields{f});
end
xlabel('Wavelength (nm)','Interpreter','none'); 
ylabel('Fluorescence intensity (a.u.)','Interpreter','none');
title('Mean Spectrum 385 nm','Interpreter','none');
grid on; 
legend('show','Interpreter','none'); 
hold off

% 405 nm
subplot(2,1,2); hold on
plot(lambda, mean_fluo405, 'k', 'LineWidth', 1.5, 'DisplayName', 'Experimental Mean');
plot(lambda, mean_fit405, 'r', 'LineWidth', 1.5, 'DisplayName', 'Total Fit Mean');
for f = 1:length(fields)
    plot(lambda, mean_contrib_405(f,:), '--', 'Color', colors(f,:), 'LineWidth', 1.2, 'DisplayName', fields{f});
end
xlabel('Wavelength (nm)','Interpreter','none'); 
ylabel('Fluorescence intensity (a.u.)','Interpreter','none');
title('Mean Spectrum 405 nm','Interpreter','none');
grid on; 
legend('show','Interpreter','none'); 
hold off