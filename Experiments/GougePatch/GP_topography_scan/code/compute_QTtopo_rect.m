function [Ztopo, x1vec, x2vec, I] = compute_QTtopo_rect(Zmat, tmat, Trigmat, param)
% Compute topography from the laser displacement sensor with QT stage at NIED
%
% [Input]
% Zmat: [um] time series of laser displacement sensor
% tmat: [s] time vector of data
% Trigmat: [V] V output of QT stage used for triggering
% param.l_space: [mm] size of scanning space
% param.ds: [mm] prescribed grid size  
% param.fs: [Hz] sampling frequency of data
% param.vpps: [pulse/s] speed of motion of stage
% param.R_sample: [mm] length of sample
% param.center_x1: [mm] central location of the sample in x1
% param.center_x2: [mm] central location of the sample in x2
% param.Vout_thresh: [V] threshold to pick the motion of stage
%
% param.low_clipval: [um] lower value of clip and interpolation for missing data point
% param.gsize: [npts] number of point of gaussian filter
% param.sig: [um] std of gaussian filter

% [Output]
% x1vec, x2vec: coordinates of the image of topography 
% Ztopo : processed topography 
% I : other data for debugging
%
% [Workflow]
% 
% 1. Pick the separation of stage motion using Trigmat 
% 2. Reshape the trace
% 3. Trim the target sample
% 4. Correct the offset due to the direction of scanning   
% 5. Fit the plane to remove slope of rock surface
% 6. remove the offset on the hostrock
% 7. interpolate the missing data in the image
%
% [Note]
% Some of the parameters are fixed in the function for the sake of
% simplicity. You can change them if needed.
%
% 2023.02.28 Kurama Okubo
% 2023.12.21 update to use the cubic sample of base rock

%% Constant parameters
PPM = 500; % pulse per mm; 1mm for 500 pulse.
Planefitfactor = 0.8; % Factor of radius to fit the plane with avoiding the sample edges.
Offset_prctile = 0.1; % offset of hostrock is taken from this percentile

%% Pre-processing %%
% Shift Coordinate
x1vec = 0:param.ds:param.l_space;
x2vec = 0:param.ds:param.l_space;
x1vec = x1vec - param.center_x1;
x2vec = x2vec - param.center_x2;

% compute npts
I.npulse_l = param.l_space * PPM; % number of pulse to send per scan trace
I.npts_l = round(param.fs * I.npulse_l/param.vpps)+1; % number of data points per scan trace; note that the coordinate starts from 0.
I.Ntrace = param.l_space/param.ds + 1; % number of trace in x2 direction

% find Tstart
I.Tstart_ind = find(diff(Trigmat)>param.Vout_thresh, 1);
I.Tstart = tmat(I.Tstart_ind - round(1.0*param.fs)); % 1 second before the state of motion

% trim the time series
Zmat_trim = Zmat(I.Tstart_ind:end);
tmat_trim = tmat(I.Tstart_ind:end);
Trig_trim = Trigmat(I.Tstart_ind:end);

Ndata = length(tmat_trim);

% search initiation and termination of stage motion
Voutdiff = diff(Trig_trim);

i = 1;
iskip = 5; % skip npts after detection of onset or stop of motion
itrace = 0;
sep_init = zeros(I.Ntrace, 1);
sep_end = zeros(I.Ntrace, 1);

while (i<Ndata) && (itrace<=I.Ntrace)
    if Voutdiff(i) > +param.Vout_thresh
        itrace = itrace + 1;
        sep_init(itrace) = i
        i = i + iskip;
    elseif Voutdiff(i) < -param.Vout_thresh
        sep_end(itrace) = i;
        i = i + iskip;
    else
        i = i + 1;
    end
end

I.Ntrace
% 
% disp(size(sep_init));
% disp(size(sep_end));

% sep_init
% figure(11); clf; hold on;
% plot(tmat_trim(1:end-1), Voutdiff, "-o");
% % plot(tmat_trim(nonzeros(sep_init)), Voutdiff(nonzeros(sep_init)), "s", "Color", "r");
% plot(tmat_trim(nonzeros(sep_init)), zeros(length(nonzeros(sep_init)), 1), "s", "Color", "r");
% 
% figure(12); clf;
% plot(diff(nonzeros(sep_init)), "o-");

%% Reshape with fixed sampling frequency
sep_len = sep_end - sep_init;
sampling_pts_tolerance = 5; % check if the error of trace npts is less than this number

% abs(sep_len - I.npts_l)
% itrace

assert(all(abs(sep_len - I.npts_l) < sampling_pts_tolerance), "sampling error exceedes tolerance.");

Ntrace_trig = length(sep_init);

I.Zmat_topo = zeros(Ntrace_trig, I.npts_l);

for i = 1:Ntrace_trig
    stind = sep_init(i);
    tr = Zmat_trim(stind:stind+I.npts_l-1)';
    if mod(i, 2) == 1
%         tr = fliplr(tr); % flip the scan data on the odd traces
    end
    I.Zmat_topo(i, :) = tr;
end

%% Trim the target sample
I.Zmat_topo_trim = I.Zmat_topo;

for i = 1:I.npts_l
    for j = 1:I.Ntrace
       x1 = x1vec(i);
       x2 = x2vec(j);
       if (abs(x1) > param.R_sample) || (abs(x2) > param.R_sample) 
           I.Zmat_topo_trim(j, i) = NaN; % (x1, x2) is with (j, i)
       end
    end
end

% %% Correct the mean between the RL and LR scanning
% % There is a offset between the scanning from RL and from LR traces.
% % It could be due to the motion of laser sensor on the inclined surface.
% % To correct that, we subtract the mean of RL trace from the LR so that
% % they have the same mean.
% 
% % compute difference in mean
% I.A1 = median(I.Zmat_topo_trim(1:2:end), 'omitnan'); % mean of RL trace
% I.A2 = median(I.Zmat_topo_trim(2:2:end), 'omitnan'); % mean of LR trace
% I.offsetdiff = I.A2-I.A1;
% I.Zmat_topo_trim(2:2:end) = I.Zmat_topo_trim(2:2:end) - I.offsetdiff;

%% Fit the surface
Nzall = nnz(~isnan(I.Zmat_topo_trim));
x1all = zeros(Nzall, 1);
x2all = zeros(Nzall, 1);
zall = zeros(Nzall, 1);

k = 0;
for i = 1:I.npts_l
    for j = 1:I.Ntrace
       x1 = x1vec(i);
       x2 = x2vec(j);
       if ~isnan(I.Zmat_topo_trim(j, i))
           % % avoid edges
           % if sqrt((x1)^2 + (x2)^2) < Planefitfactor*param.R_sample; we
           % switched to the rectangular host rock. We estimate the plane
           % on the entire area of the host rock.
               k = k + 1;
               x1all(k) = x1vec(i);
               x2all(k) = x2vec(j);
               zall(k) = I.Zmat_topo_trim(j, i);
           % end
       end
    end
end

%--- reference: https://www.mathworks.com/matlabcentral/answers/357947-how-to-fit-plane-z-ax-by-c-to-3d-point-cloud-data#answer_282719
A = [x1all, x2all, ones(Nzall, 1)];
B = A \ zall;
%---%

I.Zmat_topo_planeremoval = I.Zmat_topo_trim;
for i = 1:I.npts_l
    for j = 1:I.Ntrace
        x1 = x1vec(i);
        x2 = x2vec(j);
        zz = B(1) * x1 + B(2) * x2 + B(3) * 1;
        I.Zmat_topo_planeremoval(j, i) = I.Zmat_topo_trim(j, i) - zz;
    end
end 

%% Remove the offset from host rock
Zmat_sort = sort(I.Zmat_topo_planeremoval(:));
Zmat_sort = Zmat_sort(~isnan(Zmat_sort));
offset_host = Zmat_sort(round(Offset_prctile*Nzall));

I.Zmat_topo_planeremoval_offsetremoved = I.Zmat_topo_planeremoval - offset_host;

%% interpolate the missing data

%1. clip the low value, which is a missing data due to e.g. low light strength
Zmat_topo_fillmissingval = I.Zmat_topo_planeremoval_offsetremoved;
Zmat_topo_fillmissingval(Zmat_topo_fillmissingval < param.low_clipval) = NaN;

% 2. gaussian filter with nan
% reference: https://stackoverflow.com/questions/18697532/gaussian-filtering-a-image-with-nan-in-python
% We apply conv2 two times; first one is with zero-replaced array, which
% under estimate the gaussian summation. then, we prepare the auxiliary
% matrix to rescale the value so that the sum of weight to be unity.

%1. auxiliary matrix with replacing nan to zero
Zmat_0 = Zmat_topo_fillmissingval;
Zmat_0(isnan(Zmat_0)) = 0.0;

filt_kernel = gauss_filterkernel(param.gsize, param.sig);

V = conv2(Zmat_0, filt_kernel, 'same');

%2. second auxiliary matrix with one and zero at nan
pressure_W = ones(size(Zmat_0));
pressure_W(isnan(Zmat_topo_fillmissingval)) = 0;
W = conv2(pressure_W, filt_kernel, 'same');

I.Zmat_gaussfiltered=V./W;

% Replace the clipped value with gaussian filtered distribution
for i = 1:I.npts_l
    for j = 1:I.Ntrace
        x1 = x1vec(i);
        x2 = x2vec(j);
       
        if ((abs(x1) < param.R_sample) && (abs(x2) < param.R_sample)) && (isnan(Zmat_topo_fillmissingval(j, i)))       
            Zmat_topo_fillmissingval(j, i) = I.Zmat_gaussfiltered(j, i);
        end
    end
end 

Ztopo = Zmat_topo_fillmissingval;

end

%----------------------------%
%---Gaussian filter kernel---%
%----------------------------%
function filt_kernel = gauss_filterkernel(gsize, sig)

assert(gsize>1, 'Size of gaussian filter must be lager than 1.');
filt_kernel=zeros(gsize, gsize);

for i = 1:gsize
    for j = 1:gsize
        x = i-gsize/2-0.5;
        y = j-gsize/2-0.5;
        if sqrt(x^2 + y^2) <= gsize/2 % crop the circular area to avoid the asymmetry of edges outside the circle.
            filt_kernel(i, j) = 1/(2*pi*sig^2)* exp(-(x^2 + y^2)/(2*sig^2));
        end
    end
end

filt_kernel = filt_kernel/sum(filt_kernel, 'all');
end

