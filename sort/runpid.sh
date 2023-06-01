#!/bin/bash
#===================================================================================================================
#Program to run Clarion2-trinity Data Processing
#v1.0 01/04/2022
#C.Wibisono

#How to Run:
#./runall.sh [initial_run_number] [final_run_number] [dir_initial] [date_run] [dateseg] [bangatefile]
#====================================================================================================================

#defining the run number to be filtered:
initial=$1
final=$2

#defining the directory of initial evt file:
directoryinitial=$3

#Date
daterun=$4

#SegmentDate
dateseg=$5

#Background:
bangate=$6

if ([ $initial -ge 10  ] && [ $initial -lt 100 ])
then 
for ((j=initial;j<=final;j++));
	do
	if [ -f  "$directoryinitial/Dec$daterun-$dateseg-0$j.evt.to" ]
	then
	#ls -l $directoryinitial/Dec$daterun-$dateseg-$j.evt.to;
	../src/./clarion_pid -up $directoryinitial/Dec$daterun-$dateseg-0$j.evt.to ../../cal_fsu.ca3 ../../id_fsu_Jun9.map gagg_calib_proton.txt $bangate 0.8 1 proton.txt gamgate.txt;  
  fi
	done
elif [ $initial -lt 10 ]
then
for ((j=initial;j<=final;j++));
  do
  if [ -f "$directoryinitial/Dec$daterun-$dateseg-00$j.evt.to" ]
  then
  ../src/./clarion_pid -up $directoryinitial/Dec$daterun-$dateseg-00$j.evt.to ../../cal_fsu.ca3 ../../id_fsu_Jun9.map gagg_calib_proton.txt $bangate 0.8 1 proton.txt gamgate.txt;
  fi
  done
else
for ((j=initial;j<=final;j++));
  do
  if [ -f "$directoryinitial/Dec$daterun-$dateseg-$j.evt.to" ]
  then
  ../src/./clarion_pid -up $directoryinitial/Dec$daterun-$dateseg-$j.evt.to ../../cal_fsu.ca3 ../../id_fsu_Jun9.map gagg_calib_proton.txt $bangate 0.8 1 proton.txt gamgate.txt;
  fi
  done
fi


