#!/bin/bash
#! /usr/local/bin/bash -l 
#


for((k=1000;k<1001;k++))
do
    
    cat Energy_studies_recoquant_iNUM.C | sed -e "s/iNUM/${k}/g" >  Energy_studies_recoquant_${k}.C;
     ###if these steps have been completed then comment out these lines and comment in the lines above.. 
    cat loadlib_iNUM.C | sed -e "s/iNUM/${k}/g" > loadlib_${k}.C;
    #____ for numu ____
    cat run_iNUM.sh  | sed -e "s/iNUM/${k}/g" > run_${k}.sh; 
    qsub run_${k}.sh
    #____ for nue ____
    #cat run_nue_iNUM.sh  | sed -e "s/iNUM/${k}/g" > run_nue_${k}.sh; 
    #qsub run_nue_${k}.sh
    #____ for antinumu ____
    #cat run_antinumu_iNUM.sh  | sed -e "s/iNUM/${k}/g" > run_antinumu_${k}.sh; 
    #qsub run_antinumu_${k}.sh
    #____ for antinue ____                                                                                                         
    #cat run_antinue_iNUM.sh  | sed -e "s/iNUM/${k}/g" > run_antinue_${k}.sh;
    #qsub run_antinue_${k}.sh 

done
