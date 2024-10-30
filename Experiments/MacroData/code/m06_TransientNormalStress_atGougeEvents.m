% Evaluate the transient normal stress at the gouge events
% to investigate the effect in the non-self-similarity of those events 

clear all;
set(0,'DefaultTextFontsize',16, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Arial', ...
    'DefaultAxesFontsize',16, ... 
    'DefaultAxesFontname','Arial', ...
    'DefaultLineLineWidth', 1.0)

figdir="../figure";
if ~exist(figdir) mkdir(figdir); end
addpath("../../../utils/matlabcode_biax_v03");


%% Load event High-frequency sample data

ifSavefig = true;
ifPlotAE = true; % true if plotting AE

% ifBandpass = 1; % Apply band-pass to find ordinary events.
% ifLFEpass = 0; % Set low-frequency range to detect LFEs.
event_id = 26;
% for event_id = 1:56

    expr_id = 87;
    runID = sprintf('FB03_%03d', expr_id);
    T = readtable(sprintf("../../DetectEvent/data/p02_eventtype_fb03-%03d.csv", expr_id));

    event_type = zeros(length(size(T, 1)), 3); %1. event type %2. rupture vel %3. event start time from strain

    % Read the picked csv file to superimposed the detected events.
    A = readtable("../../DetectEvent/data/p06_visual_pick_gougeevents_merged.csv", 'NumHeaderLines',5);

    % filter by expr_id
    idx_expr = string(A.Var1) == sprintf("fb03-%03d", expr_id);
    A_expr = A(idx_expr,:);


    %%
    % for event_id = 1:size(T, 1)

    % filter by the event
    idx_ev = A.Var2 == event_id;
    if ~isempty(A_expr)
        A_event = A_expr(idx_ev,:);
    else
        A_event = [];
    end

    foreshock_pt = 0; %[ms] visually pick from the plot

    event_datdir=sprintf("/Volumes/4mGouge_WorkHDD/FB03data/4mBIAX_paper_tmp/p03_eventdata_FB03_%03d", expr_id);

    load(event_datdir + sprintf("/eventdata_FB03_%03d_event%02d.mat", expr_id, event_id));

    %%
    % obtain the tri-axial strain gouge near the target patch
    gouge_x = 1750;
    search_x = 200; %[mm] search the SGTs within this distance

    SGT_selected = find(abs(SGT_x-gouge_x) < search_x);

    %%

    plot(tmat_strain_event, Snmat(:, 1:16));
    legend()