clear all;
theta = linspace(1, 85,101);
a = 3.2725e-09 %1e-4; %1e-4;
b = 3.8678 %2.5;
y = exp(-a.*theta .^b);

fig = figure(1);
fig.Units = 'normalized';
fig.Position = [0 1 0.4 0.6];
clf(fig,'reset'); cla(fig,'reset'); clf;
hold on;
plot(theta, y, "+-");
% plot(theta, y .^ 2, "r+-");
ylim([0, 1])