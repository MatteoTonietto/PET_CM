function yCt = modelCt_num(A,a,info,t)

% Unpack Cp
tCp = info.Blood.tCp - info.delay;
Cp  = info.Blood.Cp;

% Memory allocation
Gt = zeros(length(t),length(A));

% Model construction
%--------------------------------------------------------------------------
for j = 1 : length(a)
    Gt(:,j) = interp1(tCp,(tCp(2)-tCp(1)) * filter(exp(-a(j)*tCp).*(tCp>=0),1,Cp),t);
end

if length(A) > length(a)
    Gt(:,end) = info.Cb_int;
end

% Output
%--------------------------------------------------------------------------
yCt  = Gt * A;