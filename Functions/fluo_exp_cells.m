function fluo=fluo_exp_cells(param,fluorophores,lambda)   
    pd=makedist('Normal',param(7),param(8));
    pd2=makedist('Normal',param(9),param(10));
    pd3=makedist('Normal',param(11),param(12));
    pd4=makedist('Normal',param(13),param(14));
    pd5=makedist('Normal',param(16),param(17));
    FMN=1000*pdf(pd,lambda);
    lipo=1000*pdf(pd2,lambda);
    PpIX_636=1000*pdf(pd3,lambda);
    PpIX_620=1000*pdf(pd4,lambda);
    NADH_bound=1000*pdf(pd5,lambda);
    fluo=param(1)*fluorophores(1,:)+param(2)*fluorophores(2,:)+param(3)*FMN+param(4)*lipo+param(5)*PpIX_636+param(6)*PpIX_620+param(15)*fluorophores(3,:)+param(18)*NADH_bound;
end