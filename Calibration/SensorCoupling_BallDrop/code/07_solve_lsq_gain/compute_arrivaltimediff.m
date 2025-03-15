clear all
theta = linspace(0, 85,101);
h = 0.05; %m
cp = 6919; %m/s
d0 = h ./ cosd(theta);
% figure(1); clf;
% plot(d0);

d1 = sqrt(d0 .^ 2 + 16 * h^2 - 8 .* d0 * h .* cosd(theta)); 
del_t = (d1 - d0) / cp;

fig = figure(2);
fig.Units = 'normalized';
fig.Position = [0 1 0.4 0.4];
clf(fig,'reset'); cla(fig,'reset'); hold on;clf

plot(theta, del_t*1e6, 'k-');
xlabel("Incident angle");
ylabel("Arrival time difference [μs]");
axis square
grid on;
box on;
set(gcf, 'Color', 'w');
export_fig("../../figure/directionality/reflection_arrivaltimediff.png");