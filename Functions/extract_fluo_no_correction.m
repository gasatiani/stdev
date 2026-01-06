function [S385total,S405total,fluorophore,lambda,res_385,res_405,residu] = extract_fluo_no_correction(path,file,min_wl_fluo,max_wl_fluo,fluorophores_385,fluorophores_405)
    %% obtain fluo data
    [S385total,S405total,lambda]=load_fluo(fullfile(path,file),min_wl_fluo,max_wl_fluo);
    lb = [0,0,0,480,10,0,590,10,0];
    ub = [1000,1000,1000,500,20,1000,590,10,1000];
    
    options=optimoptions('lsqcurvefit',...
        'MaxFunctionEvaluations',1e10,...
        'FunctionTolerance',1e-10, ...
        'MaxIter',1e8, ...
        'StepTolerance',1e-10); 
    fit_405=@(x,lambda)fluo_exp_405(x,fluorophores_405,lambda);
    fit_385=@(x,lambda)fluo_exp_385(x,fluorophores_385,lambda);
    
        lb = [0,0,0,0,0,0,494,14,589,9,636,5.5,618,7.5,0];
    ub = [1000,1000,1000,1000,1000,1000,496,16,591,11,638,7.5,620.5,9,1000];
    x0=[0.1,0.1,0.1,0.1,0.1,0.1,495,10,590,10,619,6,636,6];

    [res_385,residu.fit_385] = lsqcurvefit(fit_385,x0,lambda,S385total,lb,ub,options);
    [res_405,residu.fit_405] = lsqcurvefit(fit_405,x0,lambda,S405total,lb,ub,options);
 
     
    fluorophore.flavine_385=fit_385([res_385(1),0,0,0,0,0,1,1,1,1,1,1,1,1],lambda);
    fluorophore.flavine_405=fit_405([res_405(1),0,0,0,0,0,1,1,1,1,1,1,1,1],lambda);
    fluorophore.NADH_385=fit_385([0,res_385(2),0,0,0,0,1,1,1,1,1,1,1,1],lambda);
    fluorophore.NADH_405=fit_405([0,res_405(2),0,0,0,0,1,1,1,1,1,1,1,1],lambda);
    fluorophore.gaussian_385=fit_385([0,0,res_385(3),0,0,0,res_385(7),res_385(8),1,1,1,1,1,1],lambda);
    fluorophore.gaussian_405=fit_405([0,0,res_405(3),0,0,0,res_385(7),res_385(8),1,1,1,1,1,1],lambda);
    fluorophore.lipo_385=fit_385([0,0,0,res_385(4),0,0,1,1,res_385(9),res_385(10),1,1,1,1],lambda);
    fluorophore.lipo_405=fit_405([0,0,0,res_405(4),0,0,1,1,res_385(9),res_385(10),1,1,1,1],lambda);
    fluorophore.PpIX_636_385=fit_385([0,0,0,0,res_385(5),0,1,1,1,1,res_385(11),res_385(12),1,1],lambda);
    fluorophore.PpIX_636_405=fit_405([0,0,0,0,res_405(5),0,1,1,1,1,res_385(11),res_385(12),1,1],lambda);
    fluorophore.PpIX_620_385=fit_385([0,0,0,0,0,res_385(6),1,1,1,1,1,1,res_385(13),res_385(14)],lambda);
    fluorophore.PpIX_620_405=fit_405([0,0,0,0,0,res_405(6),1,1,1,1,1,1,res_385(13),res_385(14)],lambda);
    fluorophore.fit_385=fit_385(res_385,lambda);
    fluorophore.fit_405=fit_405(res_405,lambda);
end