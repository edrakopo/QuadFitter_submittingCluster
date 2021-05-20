#!/bin/sh
#$ -l h_vmem=8G #this make scripts aborted?!
#$ -N ener_edit
#$ -e /home/edrakopo/using_titus_ppe/running_ener/logs/
#$ -o /home/edrakopo/using_titus_ppe/running_ener/logs/

# Setup work directory
#myID=`echo $PBS_JOBID  | cut -f 1 -d .`
#workdir=/scratch/$USER/$myID
#mkdir -p $workdir
#echo "Workdir " $workdir

#path to input scripts:
Cscript_path=/home/edrakopo/using_titus_ppe/running_ener

# Where do I put the output?
datadir=/exports/csce/datastore/ph/groups/PPE/titus/ts-WChRecoSandBox/scripts/editing_ene/outputs
echo "datadir "$datadir

#cd $workdir
mkdir -p nu_numu_ener_iNUM
cd nu_numu_ener_iNUM
echo "thedir is " `pwd`

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/edrakopo/using_titus_ppe/lib_fromppe
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/exports/csce/datastore/ph/groups/PPE/titus/root/lib/
source /exports/csce/datastore/ph/groups/PPE/titus/ts-titus/Source_At_Start.sh

cp $Cscript_path/loadlib_iNUM.C .
cp $Cscript_path/Energy_studies_recoquant_iNUM.C .
echo "scripts and input files copied"
root -b -q loadlib_iNUM.C Energy_studies_recoquant_iNUM.C
echo "scripts run"

##rm FullEvent.root
##rm generatorcardfile.root

cp *.root $datadir

rm *.C
rm *.root
cd ..
rm -r nu_numu_ener_iNUM
echo "end at dir " `pwd`

#cp -r nu_numu_photmask_iNUM $datadir
#rm -r nu_numu_photmask_iNUM
#rm *.root
#cd $datadir
#cp -f  $workdir/*.root .
#cp -f  $workdir/stdout  ${myID}_stdout
#rm $workdir

