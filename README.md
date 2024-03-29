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

./clarion -up data.evt.to cal_fsu.ca3 id_fsu_Jun9.map gagg_calib_proton.txt 2d_all_p.banx 0.8 1 proton.txt

Note about scaling factor flag:
It is used as an empirical factor needed for kinematic corrections. For the use of simple Doppler correction replace the scaling factor with an empirical beta value.

The output for the programs will generate multiple files which are basically a histogram represented as a 2D matrix. Please see sort/Clarion2_SpectraFiles.xlsx to see the lists of the files along with the notes about dimensions and descriptions.

#Note for the PID:
PIDs for Clarion2-Trinity are made separately for each GAGG. Gamma gate condition can also be made to give more restrictive requirement so that the separation between charged particles are more distinct. In order to produce the PIDs, one needs to compile differently with the one previously mentioned. To get the PIDs do the following:

#Compile:
./build_pid.sh

#Linking:
./link_pid.sh

The procedure for running the program is the same except that now the executable name is clarion_pid. Also, there is an additional argument next to reaction channel parameter that is the information about the gamma gate energy window that we want to put to create gamma-gated pid. 

#Example of how to generate PID:

./clarion_pid -up data.evt.to cal_fsu.ca3 id_fsu_Jun9.map gagg_calib_proton.txt 2d_all_p.banx 0.8 1 proton.txt gamgatefile.txt

The default dimensions for each PID is 4096 x 4096.

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

6). User can rebin the histogram. 

Note that visualization program is far from perfection and is currently still under development.

**Visualization[Scripting Version]**

Another way to visualize the output files is by using several classes that are provided in the directory visualization/custom. In order to use the feature of classes listed there, knowledge about Python Object Oriented Programming is required. A user is also expected to know the Object oriented feature of Matplotlib such as figure and axes objects that are heavily used for plotting.
