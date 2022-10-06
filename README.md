# Clarion2-Trinity
**Clarion2-Trinity Data Processing Software and Visualization**


**#Clarion2-Trinity Data Processing:**
#How to Compile:
./build.sh

After all of the programs and their dependencies have been compiled, link the object by typing:
./link.sh

#How to Run:
./clarion -up [time-order-file] [calibration file] [detector mapping file] [gagg calibration file] [banana cut file] [scaling factor for doppler correction] [option of doppler correction;1 for kinematic correction;everything else can be made 2]

#Example of how to run the program:
./clarion -up data.evt.to cal_fsu.ca3 id_fsu_Jun9.map gagg_calib_proton.txt 2d_all_p.banx 0.8 1

The output for the program will generate multiple files which are basically a histogram represented as 2D matrix. 

**#Visualization[Very Beta version]** 
In order to visualize the output files, one can use the program called vis.py

#How to Run:
./vis.py [outputfile] [ydim] [xdim] [xlowerlimit] [xupperlimit] [ylowerlimit] [yupperlimit]

#Example:
./vis.py pid_qdc21.spn 4096 4096 0 4095 0 4095
Follow the instruction on the screen to project on X or Y axis or to set the limit of the 2D histograms.

Note that visualization program is far from perfection and is currently still under development. 
