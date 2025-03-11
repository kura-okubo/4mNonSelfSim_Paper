% Script to stack the signal of LDV from transducer on metal block.
% read data from the triggered traces
% 2022/07/01 Kurama Okubo
% 2022/08/10 Update file name for the validation
% 2022/09/05 Update for master data

clear all;
set(0,'DefaultTextFontsize',14, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Arial', ...
    'DefaultAxesFontsize',14, ... 
    'DefaultAxesFontname','Arial', ...
    'DefaultLineLineWidth', 1.5)
%% Parameters
sectionlength = 2^14; % length of trace set in Sbench configuration.
Nevent_LDV = 10000;
Nevent_AE = 2000;
input_ohm = 50; % input inpedence
fs = 1e7;

fodir = "/Volumes/Okuboetal2025_masterHDD/AEsensorCalibration/AE_FreqResp_Mat_Masterdata";
if ~isfolder(fodir); mkdir(fodir); end

figdir = "../figure/debug_convertdata";
if ~isfolder(figdir); mkdir(figdir); end


% casename_list=["master_S25V", "master_S100V", "master_S150V", "master_S200V",...
%     "termination_1Mohm_S200V", "testdata_S200V", "testdata_side_S200V", "granite_S200V"];

casename_list=["master_S200V", "testdata_S200V", "testdata_side_S200V"]; % "master_S25V", "master_S100V",  "master_S150V"]; %dump only the test cases

casename_list_out=["frontcenter_S200V", "fronttop_S200V", "sidecenter_S200V"]; % "frontcenter_S25V", "frontcenterS100V",  "frontcenter_S150V"];

for ic =1:length(casename_list)
% ic = 4;
casename = casename_list(ic);
finame = sprintf("%s/LDV_%s.mat", fodir, casename);

if casename == "termination_1Mohm_S200V"
    gain_LDV = 0.01; % [(m/s)/V] with 1M ohm
else
    gain_LDV = 0.02; % [(m/s)/V] with 50 ohm
end

A = load(finame);
tr_stacked = mean(A.datmat);

S.t = A.tmat;
S.LDV =  mean(A.datmat)' * gain_LDV;
Stable = struct2table(S);
foname = sprintf("../data/LDV_%s.csv", casename_list_out(ic));
writetable(Stable,foname);

end

%% linear stack


Nevent_linprog = 18000; % 2000 x 9 cases
input_ohm = 50; % input inpedence
gain_LDV = 0.02;

V = [10, 25, 50, 75, 100, 125, 150, 175, 200]; %[V]

casename="linearprog";
Nevent_linprog = 18000; % 2000 x 9 cases

finame = sprintf("%s/LDV_%s.mat", fodir, casename);
A = load(finame);

S1.t = A.tmat;

for j = 1:length(V)
    sind = 2000*(j-1)+1;
    eind = 2000*j; 
    tr_LP = mean(A.datmat(sind:eind, :), 1) * gain_LDV;
    S1.(sprintf("V%d", V(j))) =  tr_LP';
end

Stable = struct2table(S1);
foname = sprintf("../data/LDV_%s.csv", casename);
writetable(Stable,foname);