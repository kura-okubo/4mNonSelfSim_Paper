function [vpps, delt, loop, T] = get_QTparam(ds, fs, l)
% return program parameter for QT edit
% ds: grid space (assuming uniform in x and y)
% fs: sampling frequency of TEAC
% l: length of scan space
% OUTPUT
% vpps: speed of motion of stage [pulse per s]
% delt: time to complete one subroutine [s]
% loop: count of loop for subroutine
% T: total time to complete the scan [s]
%
% 500 pulses = 1mm

vpps = ds*fs*500; %[pulse/sec]
delt = l * 500 / vpps;

loop = ceil(l/(2*ds));

T = l*delt/ds;

if vpps>3400
    warning("vpps is more than the limitation of QT stage");
end

end