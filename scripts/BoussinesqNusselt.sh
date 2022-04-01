#!/bin/bash -e

# Nusselt number for RB convection

if [ "$#" -ne 3 ]
then
   echo usage: BoussinesqNusselt.sh case startTime startAveraging
   exit
fi

case=$1
startTime=$2
startAveraging=$3

# Calculate vertical heat transport and global mean
multiFluidHeatTransferFoam -case $case -time $startTime':'
writeuvw heatTransfer -case $case -time $startTime':'
globalSum heatTransferz -case $case
boundarySum heatTransferz -case $case ground
boundarySum heatTransferz -case $case top

# Calculate time means
timeMean -time $startAveraging':' heatTransferz
globalSum heatTransferzTimeMean -case $case -time $startAveraging':'
boundarySum heatTransferzTimeMean -case $case ground -time $startAveraging':'
boundarySum heatTransferzTimeMean -case $case top -time $startAveraging':'

# Find b0, bTop, z0 and zTop to non-dimensionalise
meshBounds=`checkMesh -case $case -constant \
    | grep 'Overall domain bounding box'`
z0=`echo $meshBounds | awk '{print $7}' | awk -F')' '{print $1}'`
zTop=`echo $meshBounds | awk '{print $10}' | awk -F')' '{print $1}'`
boundarySum b -case $case ground -time 0
boundarySum b -case $case top -time 0
b0=`tail -1 $case/sum_groundb.dat | awk '{print $5}'`
bTop=`tail -1 $case/sum_topb.dat | awk '{print $5}'`
dBdZ=`awk BEGIN'{print ('$bTop'-'$b0')/('$zTop'-'$z0')}'`
alpha=`grep alpha $case/constant/environmentalProperties | head -1 |
       awk '{print $8}' | awk -F';' '{print $1}'`
heatTransferNorm=`awk 'BEGIN{print -1*'$dBdZ'*'$alpha'}'`
echo b goes from$b0 at $z0 to $bTop at $zTop. dBdZ = $dBdZ. heatTransferNorm = $heatTransferNorm

# Plot global and boundary averages and time mean of global average
gv=0
inputFiles=($case/globalSumheatTransferz.dat
            $case/globalSumheatTransferzTimeMean.dat
            $case/sum_groundheatTransferz.dat
            $case/sum_groundheatTransferzTimeMean.dat
            $case/sum_topheatTransferz.dat
            $case/sum_topheatTransferzTimeMean.dat)
outFile=$case/globalSumheatTransferz.eps
legends=("Global mean" "Global time mean" "Ground mean" "Ground time mean"
         "Top mean" "Top time mean")
pens=("black" "1.5,black,3_3:0" "red" "1,red,3_3:1.5" "blue" "1,blue,3_3:3")
col=5
colx=1
xlabel="Time"
ylabel="Heat transfer ratio (Nusselt number)"
xmin=0
xmax=`tail -1 $case/globalSumheatTransferz.dat | awk '{print int($1+.999)}'`
dx=`echo $xmax | awk '{print int($1/10+.999)}'`
ymin=-20
ymax=100
dy=20
ddy=$dy
dyg=$ymax
xscale=/1
yscale=/$heatTransferNorm
nSkip=1
legPos=x5/6.6
source gmtPlot

awk '{print $1, $5/'$heatTransferNorm'}' $case/globalSumheatTransferzTimeMean.dat \
    > $case/NusseltTimeMean.dat
heatTransfer=`tail -1 $case/globalSumheatTransferzTimeMean.dat | awk '{print $5}'`
Nusselt=`awk 'BEGIN{print '$heatTransfer'/'$heatTransferNorm'}'`
echo Time and space averaged heat transfer = $heatTransfer
echo Nusselt number = $Nusselt
ev $outFile

