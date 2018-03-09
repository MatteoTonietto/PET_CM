% Returns the integral of Cb for each PET frame, eventually correcting for
% the delay.

function Cb_int = calcCbint(tCb,Cb,frameTimes,delay)

if nargin < 4
    delay = 0;
end

Cb_int       = zeros(size(frameTimes,1),1);
t_Cb         = tCb - delay;
doub_t       = find(diff(t_Cb)==0);
t_Cb(doub_t) = t_Cb(doub_t) + 1e-6;

for i = 1 : size(frameTimes)
    
    idx     = find(t_Cb >= frameTimes(i,1) & t_Cb < frameTimes(i,2));    
    t_frame = t_Cb(idx);
    c_frame = Cb(idx);
    
    if isempty(idx) || t_Cb(idx(1)) ~= frameTimes(i,1)
        val_1   = interp1(t_Cb,Cb,frameTimes(i,1),'linear','extrap');
        t_frame = [frameTimes(i,1);t_frame];
        c_frame = [val_1;c_frame];
    end
    
    if isempty(idx) || t_Cb(idx(end)) ~= frameTimes(i,2)
        val_2   = interp1(t_Cb,Cb,frameTimes(i,2),'linear','extrap');
        t_frame = [t_frame;frameTimes(i,2)];
        c_frame = [c_frame;val_2];
    end
    
    Cb_int(i) = trapz(t_frame,c_frame)./(frameTimes(i,2)-frameTimes(i,1));
end
