 function C = compart_2TCM(par,info,t)

if info.modelFormat == 'K'
    K1 = par(1);
    k2 = par(2);
    k3 = par(3);
    k4 = par(4);
    Vb = par(5);

    s1  = (-(k2 + k3 + k4) + sqrt((k2 + k3 + k4)^2 - 4*k2*k4))/2;
    s2  = (-(k2 + k3 + k4) - sqrt((k2 + k3 + k4)^2 - 4*k2*k4))/2;
    
    A11 = K1*k3/(s1 - s2);
    A12 = -A11;
    
    A21 = K1*(k4 + s1)/(s1 - s2);
    A22 = -K1*(k4 + s2)/(s1 - s2);

else
    error('Exponential format not supported yet')
end

A1 = [(1 - Vb)*A11; (1 - Vb)*A12];
A2 = [(1 - Vb)*A21; (1 - Vb)*A22];
a  = -[s1; s2];

if strcmp(info.inputFormat,'ana')
    C1 = modelCt_ana(A1,a,info,t);
    C2 = modelCt_ana(A2,a,info,t);
    
elseif strcmp(info.inputFormat,'num')
    C1 = modelCt_num(A1,a,info,t);
    C2 = modelCt_num(A2,a,info,t);
    
else
    error('input format not recognized')
    
end

C = [C1,C2];