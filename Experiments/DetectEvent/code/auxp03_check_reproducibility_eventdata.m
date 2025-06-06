% We check if the event data is reproducible from the raw continuous data.

% update 2025/3/7 We generated the master MAT data with modifications to 
% the elastic modulus to be consistent with the 4mNonSelfSim analysis.

clear all;

rootdir = "/Volumes/4mGouge_WorkHDD/_FB03data/4mBIAX_paper_tmp/";
masterdir = "/Volumes/Okuboetal2025_masterHDD/4mBIAX_eventdata_master/";
for eventid = 1:56
% eventid = 1;

Orig = load(rootdir+sprintf("p03_eventdata_FB03_087/eventdata_FB03_087_event%02d.mat", eventid)); % original data
Rep = load(masterdir+sprintf("p03_eventdata_FB03_087/eventdata_FB03_087_event%02d.mat", eventid)); % reproduced master data by rerunning the script to be uploaded to the archive

%%
assert(all(Orig.SGB_x == Rep.SGB_x));
assert(all(Orig.SGT_x == Rep.SGT_x));
assert(all(Orig.Disp_x == Rep.Disp_x));
assert(Orig.Tstart == Rep.Tstart);

%%
eps = 1e-15;
assert(norm(Orig.AEdatmat - Rep.AEdatmat) < eps);
assert(norm(Orig.DXeast - Rep.DXeast) < eps);
assert(norm(Orig.DXwest - Rep.DXwest) < eps);
assert(norm(Orig.Dmat_event - Rep.Dmat_event) < eps);
assert(norm(Orig.NPmacro - Rep.NPmacro) < eps);
assert(norm(Orig.SSmacro - Rep.SSmacro) < eps);
% Skip the stress due to the change of elastic modulus
% assert(norm(Orig.Snmat(~isnan(Orig.Snmat)) - Rep.Snmat(~isnan(Rep.Snmat))) < eps);
% assert(norm(Orig.Spmat(~isnan(Orig.Spmat)) - Rep.Spmat(~isnan(Rep.Spmat))) < eps);
% assert(norm(Orig.Taumat2(~isnan(Orig.Taumat2)) - Rep.Taumat2(~isnan(Rep.Taumat2))) < eps);
% assert(norm(Orig.Taumat3(~isnan(Orig.Taumat3)) - Rep.Taumat3(~isnan(Rep.Taumat3))) < eps);
assert(norm(Orig.tmat_AE_event - Rep.tmat_AE_event) < eps);
assert(norm(Orig.tmat_macro - Rep.tmat_macro) < eps);

fprintf("eventdata %02d is reproducible from the raw continuous data!\n", eventid);
end
