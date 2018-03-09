function yCt = model_2TCM1K(par,info,t)

if info.modelFormat == 'K'
    K1 = par(1);
    k2 = par(2);
    k3 = par(3);
    k4 = par(4);
    Ki = par(5);
    Vb = par(6);
    
    a1 = ((k2 + k3 + k4) - sqrt((k2 + k3 + k4)^2 - 4*k2*k4))/2;
    a2 = ((k2 + k3 + k4) + sqrt((k2 + k3 + k4)^2 - 4*k2*k4))/2;
    A1 = K1 * (a1 - k3 - k4)/(a1 - a2);
    A2 = K1 * (k3 + k4 - a2)/(a1 - a2);    
else
    A1 = par(1);
    a1 = par(2);
    A2 = par(3);
    a2 = par(4);
    Ki = par(5);
    Vb = par(6);
end

A = [(1 - Vb)*A1; (1 - Vb)*A2; Vb*Ki; Vb];
a = [a1; a2; 0];

if strcmp(info.inputFormat,'ana')
    yCt = modelCt_ana(A,a,info,t);
    
elseif strcmp(info.inputFormat,'num')
    yCt = modelCt_num(A,a,info,t);
    
else
    error('input format not recognized')
    
end