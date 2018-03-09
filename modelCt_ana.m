function yCt = modelCt_ana(A,a,info,t)

l0 = info.Blood.MB_Cp.par(1);
t0 = info.Blood.MB_Cp.par(2) - info.delay;
T  = info.Blood.MB_Cp.par(3);

A_good = info.Blood.MB_Cp.info.A_good;
l_good = info.Blood.MB_Cp.info.l_good;

% Memory allocation
Gt = zeros(length(t),length(A));

% Model construction
%--------------------------------------------------------------------------
for j = 1 : length(a)
    l0_Pole = SimplePoleDelayConvRectConvExp([l0,t0,T,a(j)],t);
    for i = 1 : length(l_good)
        Gt(:,j) = Gt(:,j) + A_good(i)*(SimplePoleDelayConvRectConvExp([l_good(i),t0,T,a(j)],t) - l0_Pole);
    end
end

if length(A) > length(a)
    Gt(:,end) = repmat(info.Cb_int,length(t)/length(info.Cb_int),1);
end

% Output
%--------------------------------------------------------------------------
yCt  = Gt * A;

