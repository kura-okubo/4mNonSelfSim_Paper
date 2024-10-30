# Detect the stick-slip events

This repository performs the detection of events using the slip and strain, and the evaluation of rupture type such as the location of initiation and the rupture velocity. The goal of this repository is to make the event catalog of gouge-mediated seismic events.

1. Detect the event using local slip sensors
`p01_detectevent_slipevolution.m` computes the evolution of slip during the stick-slip experiments. We detect the event by the transient jump of slip caused by the stick-slips. This provides the first-order estimation of the event time, whereas we need to re-evaluate the initiation of the nucleation using the strains to estimate the initiation of slow nucleation, which is insensitive to the slip evolution.

2. `p02_evalevent_strain.m` evaluates the event using the first-order estimation of the event time. We load the time window around the event time and estimate the following:
- Time for the initiation of nucleation
- Location of rupture nucleation
- Rupture velocity

3. `p03_save_eventdata.m` compiles the data of slip, strain and AE traces with the event time window and dump it to `mat` file.

4. `p04_plot_slipevents.m` plots the slip, strain and AE of the stick-slip event. This is used to evaluate the foreshocks to be analyzed. `p04_2_plot_slownucleation.m` plots some events with very slow nucleation so that we need to increase the range of plot.

5. `p05_expr_overview.m` computes the overall stats of the experiments such as the shear loading history, slip, stress drop, rupture type and the rupture velocity of the nucleation. We also use `p05_2_expr_overview_presentation.m` to reduce the plot tile for the use of presentation.

6. `p06_search_gougeevents.m` plot the AE traces to visually search the gouge-mediated seismic events. We get the Data tip from Matlab figure and note the event time in `p06_visual_pick_foreshocks.csv` by hand. We read the picked event and superimpose them on the AE traces to avoid the missing of the picking.

`p06_2_search_gougeevents_tinyevents.m` is tuned to find the tiny events, where we increase the freqmax and plot amplification.

`p07_plot_pickedgougeevents.m` plots all the events with the high and low frequency band-pass filters.

`data/p06_visual_pick_gougeevents_merged.csv` contains the pick time of the gouge-mediated seismic events with their id. Note that we have two event ids; one is for the stick-slip events, and another is the gouge-mediated seismic events, which are named as "event_id" and "gougeevent_id", respectively.

**NOTE:** We first analyzed the ordinary events, and we added the tiny events later to complete the gouge events. We reordered the events in the chronological order, which is documented in the `p06_visual_pick_gougeevents_merged.csv`.

In this study, the low-frequency events are out of scope. We picked them just for the sake of future works.

