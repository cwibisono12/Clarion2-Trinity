# Clarion2-Trinity
**Clarion2-Trinity Data Processing Software and Visualization**


**Clarion2-Trinity Data Processing:**

#How to Compile:
./build.sh

After all of the programs and their dependencies have been compiled, link the object by typing:
./link.sh

#How to Run:
./clarion -up [time-order-file] [calibration file] [detector mapping file] [gagg calibration file] [banana cut file] [scaling factor for doppler correction] [option of doppler correction;1 for kinematic correction;everything else can be made 2] [reaction channel information]

#Example of how to run the program:
./clarion -up data.evt.to cal_fsu.ca3 id_fsu_Jun9.map gagg_calib_proton.txt 2d_all_p.banx 0.8 1 particle.txt

The output for the program will generate multiple files which are basically a histogram represented as 2D matrix. 

#Note for the PID:
PIDs for Clarion2-Trinity are made separately for each GAGG. In order to produce the PIDs, one needs to compile differently with the one previously mentioned. To get the PIDs do the following:

#Compile:
./build_pid.sh

#Linking:
./link_pid.sh

The procedure for running the program is the same except that now the executable name is clarion_pid.

**Visualization[Very Basic version]** 

In order to visualize the output files, one can use the program called vis.py

#How to Run:
./vis.py [outputfile] [ydim] [xdim] [xlowerlimit] [xupperlimit] [ylowerlimit] [yupperlimit]

#Example:
./vis.py pid_qdc21.spn 4096 4096 0 4095 0 4095
Follow the instruction on the screen to project on X or Y axis or to set the limit of the 2D histograms.

Capabilities of visualization:

1). Enabling the user to project on X or Y axis and set the width for the projection.

2). User can also perform background subtraction if needed. (This feature is added to accomodate gamma-gamma coincidences analysis).

3). User can expand the projected region and set the lower and upper x limits to see the region of interest.

4). User can perform a Gauss fit for a peak, and get the mean, standard deviation, and area.

5). User can also draw banana gates on the 2D histogram and get the coordinates.

Note that visualization program is far from perfection and is currently still under development.
