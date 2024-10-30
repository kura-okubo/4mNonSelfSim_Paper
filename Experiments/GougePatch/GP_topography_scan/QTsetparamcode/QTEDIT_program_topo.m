% Script to compute program parameters for QT-EDIT
% 2023/02/23 Kurama Okubo
clear all;
set(0,'DefaultTextFontsize',14, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Arial', ...
    'DefaultAxesFontsize',14, ... 
    'DefaultAxesFontname','Arial', ...
    'DefaultLineLineWidth', 1.5)

%% Default Unit: [pulse], [millimeter], [second], [Hz]
l_space = 27; % [mm] total length of measurement line

%% ds=0.1
Q100.ds = 0.1;
Q100.fs = 20;

[Q100.vpps, Q100.delt, Q100.loop, Q100.T] = get_QTparam(Q100.ds, Q100.fs, l_space);

fprintf("%f, %f, %d, %f\n", Q100.vpps, Q100.delt, Q100.loop, Q100.T/60);

%% ds = 0.05
Q50.ds = 0.05;
Q50.fs = 100;

[Q50.vpps, Q50.delt, Q50.loop, Q50.T] = get_QTparam(Q50.ds, Q50.fs, l_space);

fprintf("%f, %f, %d, %f\n", Q50.vpps, Q50.delt, Q50.loop, Q50.T/60);

%% ds = 0.01
Q10.ds = 0.01;
Q10.fs = 200;

[Q10.vpps, Q10.delt, Q10.loop, Q10.T] = get_QTparam(Q10.ds, Q10.fs, l_space);

fprintf("%f, %f, %d, %f\n", Q10.vpps, Q10.delt, Q10.loop, Q10.T/60);

%% Summary

fo = fopen("../data/QT-EDIT_programtable.txt", "w");
fprintf(fo, "ds[mm], fs[Hz], vpps[pulse/s], t0[s], loop[count], T[min]\n");
for Q =[Q100, Q50, Q10]
    fprintf(fo, "%4.3f, %4.1f, %4.1f, %4.3f, %d, %4.3f\n",...
        Q.ds, Q.fs, Q.vpps, Q.delt,  Q.loop, Q.T/60);
end
fclose(fo);
