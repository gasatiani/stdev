function fluo=fluo_exp_405(param,fluorophores,lambda)   
%     FAD="D:\Lyon thèse\soft\soft manip Arthur\correction_optical_properties\experimental_data\flavine_fluo.mat";
%     NADH="D:\Lyon thèse\soft\soft manip Arthur\correction_optical_properties\experimental_data\NADH_fluo_4.mat";
%     bile=load("D:\Lyon thèse\soft\soft manip Arthur\bile_modif.mat");
%     [NADH385total,NADH405total,lambda_NADH]=load_fluo(NADH,wl_min,wl_max);
%     [FAD385total,FAD405total,lambda_FAD]=load_fluo(FAD,wl_min,wl_max);
%     lambda=bile.lambda;
%     idx_min=find(min((lambda-wl_min).^2)==(lambda-wl_min).^2);
%     idx_max=find(min((lambda-wl_max).^2)==(lambda-wl_max).^2);
    pd=makedist('Normal',param(7),param(8));
    pd2=makedist('Normal',param(9),param(10));
    pd3=makedist('Normal',param(11),param(12));
    pd4=makedist('Normal',param(13),param(14));
    FMN=1000*pdf(pd,lambda);
    lipo=1000*pdf(pd2,lambda);
    PpIX_636=1000*pdf(pd3,lambda);
    PpIX_620=1000*pdf(pd4,lambda);
    fluo=param(1)*fluorophores(1,:)+param(2)*fluorophores(2,:)+param(3)*FMN+param(4)*lipo+param(5)*PpIX_636+param(6)*PpIX_620;
end