function fluo=fluo_exp_385_vitamin_A(param,fluorophores,lambda)  
    fluo=fluo_exp_385(param(1:14),fluorophores,lambda)+param(15)*fluorophores(3,:);
end