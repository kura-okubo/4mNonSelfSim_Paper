function [R, pvel, ruptype, event_type] = rupvel_linearfit(event_id, Tstart,...
    SGB_x, SGT_x, Taumat2_removeoffset_smoothed, Taumat3_removeoffset_smoothed,...
    tmat_strain_event, x_bound_typecheck, xrange_east, xrange_west, SNratio_threshold)

% SNratio_threshold = 50;
noisewin = 60; % this noise window is different from p02_evalevent_strain.m due to the data trimming. We tuned this noisewin to synchronize the result to the p02.

NSGB = size(Taumat2_removeoffset_smoothed, 2);
NSGT = size(Taumat3_removeoffset_smoothed, 2);

R = zeros(NSGB+NSGT, 3); %1. coordinate[mm], 2. timing [s] 3. max strain value

% search SGB
for ii = 1:NSGB
    [maxeps, imax] = max(Taumat2_removeoffset_smoothed(:, ii));
    noiselevel = mean(abs(Taumat2_removeoffset_smoothed(1:noisewin, ii)));
    SNratio = maxeps/noiselevel
    if SNratio > SNratio_threshold
        R(ii, :) = [SGB_x(ii)/1e3, tmat_strain_event(imax), maxeps];
    else
        R(ii, :) = [NaN, NaN, NaN];
    end

end

% search SGT
for ii = 1:NSGT
    [maxeps, imax] = max(Taumat3_removeoffset_smoothed(:, ii));
    noiselevel = mean(abs(Taumat3_removeoffset_smoothed(1:noisewin, ii)));
    SNratio = maxeps/noiselevel;

    if SNratio > SNratio_threshold
        R(ii+NSGB, :) = [SGT_x(ii)/1e3, tmat_strain_event(imax), maxeps];
    else
        R(ii+NSGB, :) = [NaN, NaN, NaN];
    end

end

disp(R)
%% Compute linear regression

% To identify the rupture from both sides, we compare the slopes associated
% with easetern and western of x = 1250[mm], and check the difference of
% them.
% We also compute the rupture velocity of nucleation with different range
% of location for east and west nucleations.

R(any(isnan(R), 2), :) = [];

R_east = R(R(:, 1) > x_bound_typecheck(2), :);
R_west = R(R(:, 1) <= x_bound_typecheck(1), :);

pall = polyfit(R(:, 2), R(:, 1), 1);

peast = polyfit(R_east(:, 2), R_east(:, 1), 1);
pwest = polyfit(R_west(:, 2), R_west(:, 1), 1);

rupvel = pall(1); % [m/s]

% Evaluate the rupture type
% Case2 rupture propagate from west; Note that some events from west are very slow, and no data point within the
% range of polyfit. In this case, they are categorized in type 2.
if (pwest(1) == 0) && (peast(1) > 0)
    ruptype = 2;

    disp("test0")

    R_rupvel_nuc = R((xrange_west(1) <= R(:, 1)) & (R(:, 1) <= xrange_west(2)) , :);
    pvel = polyfit(R_rupvel_nuc(:, 2), R_rupvel_nuc(:, 1), 1);
    rupvel = pvel(1);

    % case 1. rupture propagate from east
elseif (pwest(1) == 0) && (peast(1) < 0)
    disp("test1")

    ruptype = 1;
    R_rupvel_nuc = R((xrange_east(1) <= R(:, 1)) & (R(:, 1) <= xrange_east(2)) , :);
    pvel = polyfit(R_rupvel_nuc(:, 2), R_rupvel_nuc(:, 1), 1);
    rupvel = pvel(1);

    % case 3. if the signs of slope in east and west are different, rupture nucleated
    % from both side
elseif sign(peast(1)) ~= sign(pwest(1))
    disp("test2")

    ruptype = 3;
    pvel = NaN;
    rupvel = NaN;

    % case 1. rupture propagate from east
elseif sign(rupvel) == -1
    disp("test3")
    ruptype = 1;
    R_rupvel_nuc = R((xrange_east(1) <= R(:, 1)) & (R(:, 1) <= xrange_east(2)) , :);
    pvel = polyfit(R_rupvel_nuc(:, 2), R_rupvel_nuc(:, 1), 1);
    rupvel = pvel(1);
    % case 2. rupture propagate from west

elseif sign(rupvel) == 1
    disp("test4")
    ruptype = 2;

    R_rupvel_nuc = R((xrange_west(1) <= R(:, 1)) & (R(:, 1) <= xrange_west(2)) , :);
    pvel = polyfit(R_rupvel_nuc(:, 2), R_rupvel_nuc(:, 1), 1);
    rupvel = pvel(1);
end

event_type(event_id, 1) = ruptype;
event_type(event_id, 2) = rupvel;
event_type(event_id, 3) = Tstart;

end