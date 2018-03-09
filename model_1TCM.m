function yCt = model_1TCM(par,info,t)

if info.modelFormat == 'K'
    K1 = par(1);
    k2 = par(2);
    Vb = par(3);

    a1 = k2;
    A1 = K1;

else
    A1 = par(1);
    a1 = par(2);
    Vb = par(3);
end

A = [(1 - Vb)*A1; Vb];
a = a1;

if strcmp(info.inputFormat,'ana')
    yCt = modelCt_ana(A,a,info,t);
    
elseif strcmp(info.inputFormat,'num')
    yCt = modelCt_num(A,a,info,t);
    
else
    error('input format not recognized')
    
end