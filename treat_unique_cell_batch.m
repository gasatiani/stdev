clear all
close all
warning('off','all')
%addpath('D:\Lyon thèse\soft\soft manip Arthur\functions')
addpath('experimental_data')
addpath("theoretical_data")
addpath('Functions/') 
dataDir='/home/asatiani/Desktop/Thesis/Experiment/optics';
% dataDir='D:\Lyon thèse\data\Manip fantomes\Phantom_20241108_old_NADH\System_008_fluo';
min_wl=430;
max_wl=650;

[data,path_selected] = uigetfile('*fluo.mat',...
   'Select One or More Files',(dataDir));
file=convertCharsToStrings(data);
[FAD385total,FAD405total,lambda,Power]=load_fluo("Reference/flavine_fluo.mat",min_wl,max_wl);
% [NADH385total,NADH405total,lambda_NADH]=load_fluo("D:\Lyon thèse\soft\soft manip Arthur\correction_optical_properties\experimental_data\NADH_fluo_4.mat",min_wl,max_wl);
[NADH385total,NADH405total,lambda_NADH]=load_fluo("Reference/NADH_fluo_4.mat",min_wl,max_wl);
% water='D:\Lyon thèse\data\Manip fantomes\Phantom_20240320_NADH\fantomes_NADH_20240321\water_probe_001_fluo.mat';
water='/home/asatiani/Desktop/Thesis/Experiment/optics/cell_cultures_37Degree/PBS_000_fluo.mat';
[water385total,water405total,lambda_water]=load_fluo(water,min_wl,max_wl);
fluorophores_385=[FAD385total;NADH385total;water385total];

fluorophores_405=[FAD405total;NADH405total;water405total];
[S385total,S405total,fluorophore,lambda,res_385_correc_Kim_full_exp,res_405_correc_Kim_full_exp,residu] = extract_fluo_no_correction_cells(path_selected,data,min_wl,max_wl,fluorophores_385,fluorophores_405);




figure
subplot(2,1,1)
p=plot(lambda,S385total,'o',lambda,fluorophore.fit_385,lambda,fluorophore.NADH_bound_385,lambda,fluorophore.NADH_385,lambda,fluorophore.flavine_385,lambda,fluorophore.gaussian_385,lambda,fluorophore.lipo_385,lambda,fluorophore.PpIX_620_385,lambda,fluorophore.PpIX_636_385,lambda,fluorophore.water_385);
grid on
yl=ylim;
% xlim([min_wl_fluo max_wl_fluo])
set(gca,'FontSize',20)

legend('Experimental data','Sum of fluorophores','NADH bound','NADH free','FAD','Protein bound FMN','Lipopigments','PpIX 620','PpIX 636','water')
title('laser 385')
xlabel('Wavelength (nm)')
ylabel('Fluorescence intensity (a.u)')
p(2).LineWidth=2;
p(3).LineWidth=2;
p(4).LineWidth=2;
p(5).LineWidth=2;
p(6).LineWidth=2;
p(7).LineWidth=2;
p(8).LineWidth=2;
p(9).LineWidth=2;
p(10).LineWidth=2;


subplot(2,1,2)
p=plot(lambda,S405total,'o',lambda,fluorophore.fit_405,lambda,fluorophore.NADH_bound_405,lambda,fluorophore.NADH_405,lambda,fluorophore.flavine_405,lambda,fluorophore.gaussian_405,lambda,fluorophore.lipo_405,lambda,fluorophore.PpIX_620_405,lambda,fluorophore.PpIX_636_405,lambda,fluorophore.water_405);
yl=ylim;
% xlim([min_wl_fluo max_wl_fluo])
grid on
legend('Experimental data','Sum of fluorophores','NADH bound','NADH free','FAD','Protein bound FMN','Lipopigments','PpIX 620','PpIX 636','water')
title('laser 405')
xlabel('Wavelength (nm)')
ylabel('Fluorescence intensity (a.u)')
p(2).LineWidth=2;
p(3).LineWidth=2;
p(4).LineWidth=2;
p(5).LineWidth=2;
p(6).LineWidth=2;
p(7).LineWidth=2;
p(8).LineWidth=2;
p(9).LineWidth=2;
p(10).LineWidth=2;




figure
subplot(2,1,1)
p=plot(lambda,S385total-fluorophore.water_385,'o',lambda,fluorophore.fit_385-fluorophore.water_385,lambda,fluorophore.NADH_bound_385,lambda,fluorophore.NADH_385,lambda,fluorophore.flavine_385,lambda,fluorophore.gaussian_385,lambda,fluorophore.lipo_385,lambda,fluorophore.PpIX_620_385,lambda,fluorophore.PpIX_636_385);
grid on
yl=ylim;
% xlim([min_wl_fluo max_wl_fluo])
set(gca,'FontSize',20)

legend('Experimental data','Sum of fluorophores','NADH bound','NADH free','FAD','Protein bound FMN','Lipopigments','PpIX 620','PpIX 636')
title('laser 385')
xlabel('Wavelength (nm)')
ylabel('Fluorescence intensity (a.u)')
p(2).LineWidth=2;
p(3).LineWidth=2;
p(4).LineWidth=2;
p(5).LineWidth=2;
p(6).LineWidth=2;
p(7).LineWidth=2;
p(8).LineWidth=2;
p(9).LineWidth=2;

subplot(2,1,2)
p=plot(lambda,S405total-fluorophore.water_405,'o',lambda,fluorophore.fit_405-fluorophore.water_405,lambda,fluorophore.NADH_bound_405,lambda,fluorophore.NADH_405,lambda,fluorophore.flavine_405,lambda,fluorophore.gaussian_405,lambda,fluorophore.lipo_405,lambda,fluorophore.PpIX_620_405,lambda,fluorophore.PpIX_636_405);
yl=ylim;
% xlim([min_wl_fluo max_wl_fluo])
grid on
legend('Experimental data','Sum of fluorophores','NADH bound','NADH free','FAD','Protein bound FMN','Lipopigments','PpIX 620','PpIX 636')
title('laser 405')
xlabel('Wavelength (nm)')
ylabel('Fluorescence intensity (a.u)')
p(2).LineWidth=2;
p(3).LineWidth=2;
p(4).LineWidth=2;
p(5).LineWidth=2;
p(6).LineWidth=2;
p(7).LineWidth=2;
p(8).LineWidth=2;
p(9).LineWidth=2;
