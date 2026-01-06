function [S385total, S405total, fluorophore, lambda, res_385_all, res_405_all, resid, fluo385_each, fluo405_each] = ...
    extract_fluo_no_correction_cells(path, file, min_wl_fluo, max_wl_fluo, fluorophores_385, fluorophores_405)

% Load fluorescence data (expects load_fluo to exist)
[S385total, S405total, lambda, fluo385_each, fluo405_each] = load_fluo(fullfile(path, file), min_wl_fluo, max_wl_fluo);

% Safety: ensure fluo_*_each are oriented [spectra x wavelengths]
if size(fluo385_each,1) < size(fluo385_each,2) && size(fluo385_each,1) == length(lambda)
    fluo385_each = fluo385_each.'; % transpose if needed
end
if size(fluo405_each,1) < size(fluo405_each,2) && size(fluo405_each,1) == length(lambda)
    fluo405_each = fluo405_each.';
end

Nspec = size(fluo385_each,1);

% Bounds, options and initial guess (same as your originals)
lb = [0,0,0,0,0,0,494,14,589,9,636,5.5,618,0,0,440,5,0];
ub = [1,1000,1000,1000,1000,1000,496,16,591,11,638,7.5,620.5,9,1000,460,15,1000];

options = optimoptions('lsqcurvefit', ...
    'MaxFunctionEvaluations', 1e10, ...
    'FunctionTolerance', 1e-6, ...
    'MaxIter', 1e8, ...
    'StepTolerance', 1e-10, ...
    'Display', 'off');

x0 = [0.1,0.1,0.1,0.1,0.1,0.1,455,10,590,10,619,6,636,6,0.1,450,10,0.1];

% Model handles (must exist in path): fluo_exp_cells
fit_385 = @(x, lambda) fluo_exp_cells(x, fluorophores_385, lambda);
fit_405 = @(x, lambda) fluo_exp_cells(x, fluorophores_405, lambda);

% Prepare outputs
fluorophore = struct();

% Total fits per spectrum
fluorophore.fit_385_each = cell(Nspec,1);
fluorophore.fit_405_each = cell(Nspec,1);

% Components per spectrum for 385 nm
fluorophore.NADH_free_385_each      = cell(Nspec,1);
fluorophore.NADH_bound_385_each     = cell(Nspec,1);
fluorophore.FAD_385_each            = cell(Nspec,1);
fluorophore.PbFMN_385_each         = cell(Nspec,1);
fluorophore.Lipo_385_each           = cell(Nspec,1);
fluorophore.PpIX620_385_each        = cell(Nspec,1);
fluorophore.PpIX636_385_each        = cell(Nspec,1);

% Components per spectrum for 405 nm
fluorophore.NADH_free_405_each      = cell(Nspec,1);
fluorophore.NADH_bound_405_each     = cell(Nspec,1);
fluorophore.FAD_405_each            = cell(Nspec,1);
fluorophore.PbFMN_405_each         = cell(Nspec,1);
fluorophore.Lipo_405_each           = cell(Nspec,1);
fluorophore.PpIX620_405_each        = cell(Nspec,1);
fluorophore.PpIX636_405_each        = cell(Nspec,1);

res_385_all = cell(Nspec,1);
res_405_all = cell(Nspec,1);
resid = struct();
resid.fit_385 = cell(Nspec,1);
resid.fit_405 = cell(Nspec,1);

% Loop fitting each spectrum separately
for i = 1:Nspec
    % Fit 385 nm
    y385 = fluo385_each(i,:);
    [res385, resnorm385, residuals385, exitflag385, output385] = lsqcurvefit(fit_385, x0, lambda, y385, lb, ub, options);
    res_385_all{i} = res385;
    resid.fit_385{i} = struct('resnorm', resnorm385, 'exitflag', exitflag385, 'output', output385);

    % Fit 405 nm
    y405 = fluo405_each(i,:);
    [res405, resnorm405, residuals405, exitflag405, output405] = lsqcurvefit(fit_405, x0, lambda, y405, lb, ub, options);
    res_405_all{i} = res405;
    resid.fit_405{i} = struct('resnorm', resnorm405, 'exitflag', exitflag405, 'output', output405);

    % Total fits
    fluorophore.fit_385_each{i} = fit_385(res385, lambda);
    fluorophore.fit_405_each{i} = fit_405(res405, lambda);

    % COMPONENTS — 385 nm
    % FAD (flavine)
    fluorophore.FAD_385_each{i} = fit_385([res385(1),0,0,0,0,0,1,1,1,1,1,1,1,1,0,1,1,0], lambda);

    % NADH free
    fluorophore.NADH_free_385_each{i} = fit_385([0,res385(2),0,0,0,0,1,1,1,1,1,1,1,1,0,1,1,0], lambda);

    % Gaussian / PbFMN (using 3rd param + gaussian position/width)
    fluorophore.PbFMN_385_each{i} = fit_385([0,0,res385(3),0,0,0,res385(7),res385(8),1,1,1,1,1,1,0,1,1,0], lambda);

    % Lipopigments
    fluorophore.Lipo_385_each{i} = fit_385([0,0,0,res385(4),0,0,1,1,res385(9),res385(10),1,1,1,1,0,1,1,0], lambda);

    % PpIX 620 & 636
    fluorophore.PpIX620_385_each{i} = fit_385([0,0,0,0,0,res385(6),1,1,1,1,1,1,res385(13),res385(14),0,1,1,0], lambda);
    fluorophore.PpIX636_385_each{i} = fit_385([0,0,0,0,res385(5),0,1,1,1,1,res385(11),res385(12),1,1,0,1,1,0], lambda);

    % NADH bound (protein-bound)
    fluorophore.NADH_bound_385_each{i} = fit_385([0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,res385(16),res385(17),res385(18)], lambda);

    % Optional: water/background if needed (kept out unless in model)
    % fluorophore.Water_385_each{i} = ... (not included unless in fluorophores arrays)

    % COMPONENTS — 405 nm (same parameter mapping)
    fluorophore.FAD_405_each{i} = fit_405([res405(1),0,0,0,0,0,1,1,1,1,1,1,1,1,0,1,1,0], lambda);
    fluorophore.NADH_free_405_each{i} = fit_405([0,res405(2),0,0,0,0,1,1,1,1,1,1,1,1,0,1,1,0], lambda);
    fluorophore.PbFMN_405_each{i} = fit_405([0,0,res405(3),0,0,0,res405(7),res405(8),1,1,1,1,1,1,0,1,1,0], lambda);
    fluorophore.Lipo_405_each{i} = fit_405([0,0,0,res405(4),0,0,1,1,res405(9),res405(10),1,1,1,1,0,1,1,0], lambda);
    fluorophore.PpIX620_405_each{i} = fit_405([0,0,0,0,0,res405(6),1,1,1,1,1,1,res405(13),res405(14),0,1,1,0], lambda);
    fluorophore.PpIX636_405_each{i} = fit_405([0,0,0,0,res405(5),0,1,1,1,1,res405(11),res405(12),1,1,0,1,1,0], lambda);
    fluorophore.NADH_bound_405_each{i} = fit_405([0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,res405(16),res405(17),res405(18)], lambda);
end

end
