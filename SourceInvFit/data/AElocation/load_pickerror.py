# Error estimation in the determination of source location
# We conducted the aux01_AE_locateevent_evalerror_GUI.py, where we mannualy pick the onset of P wave (not updated by the theoretical arrival time).
# The pick times are stored in `arrivalpick_evallocerror`.
# Then, we evaluated the std of error between the picked and theoretically calculated arrival times. The source location was optimized such that the error was minimized.
# The averaged error is printed with this script.
# 2024.2.4 Kurama Okubo 
# 2025.2.21 updated for the merged catalog.

import numpy as np
import pickle

# gougeevent_id = 16
eventid_list = [4,   9,  18,  19,  20,  21,  24,  27,  30,  31,  37,  38,  40,
        43,  44,  49,  50,  52,  55,  59,  61,  62,  69,  72,  75,
        76,  77,  81,  85,  88,  89,  95,  99, 100, 102, 109, 110, 111,
       118, 120, 126, 128, 129, 131]
std_err_all = []
err_loc_all = []
print(f"total event number: {len(eventid_list)}")
for gougeevent_id in eventid_list:
	finame_loc = f"./arrivalpick_evallocerror/{gougeevent_id:04d}/eventloc__fb03-087__{gougeevent_id:04d}.pickle"
	with open(finame_loc, "rb") as fi:
		event_loc = pickle.load(fi)

	vp = 6200
	std_err = np.sqrt(event_loc["R"]) * 1e-3 #[ms] - > [s]
	err_loc = vp*std_err * 1e3 #[mm]
	std_err_all.append(std_err)
	err_loc_all.append(err_loc)
	print(f"error location of {gougeevent_id:04d} is {std_err*1e6:.4f}us: {err_loc:.4f}mm.")

print(f"average err is {np.mean(std_err_all)*1e6:.4f}us: {np.mean(err_loc_all):.4f}mm")

with  open("./arrivalpick_evallocerror/arrivalpick_evallocerror_summary.txt", "w") as fo:
	fo.write(f"average err is {np.mean(std_err_all)*1e6:.4f}us: {np.mean(err_loc_all):.4f}mm with Cp={vp}m/s.")