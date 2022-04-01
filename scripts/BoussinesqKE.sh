#!/bin/bash -e

# Kinetic energy of multi-fluid Boussinesq equations

if [ "$#" -ne 4 ]
then
   echo usage: BoussinesqKE.sh case startTime startAveraging 'partNames'
   exit
fi

case=$1
startTime=$2
startAveraging=$3
parts=($4)
nparts=${#parts[*]}

# Calculate kinetic energy
for part in ${parts[*]}; do
    postProcess -case $case -func "magSqr(u.${part})" -time $startTime:
done
for time in $case/[0-9]*; do
    time=`filename $time`
    if [ $time -ge $startTime ]; then
        for part in ${parts[*]}; do
            multiplyFields -case $case $time sigmauSqr.${part} \
                $time magSqr\(u.${part}\) $time sigma.${part}
        done
        sumFields -case $case $time KE $time sigmauSqr.${parts[0]} \
            $time sigmauSqr.${parts[0]} -scale0 0.25 -scale1 0.25
        for part in ${parts[1]} ${parts[2]} ${parts[3]} ${parts[4]} ${parts[5]}; do
            sumFields -case $case $time KE $time KE $time sigmauSqr.$part -scale1 0.5
        done
    fi
done

# Calculate time means
timeMean -case $case -time $startAveraging':' KE

# Find nu, b0, bTop, z0 and zTop to non-dimensionalise
nu=`grep nu $case/constant/environmentalProperties | head -1 |
       awk '{print $8}' | awk -F';' '{print $1}'`
meshBounds=`checkMesh -case $case -constant \
    | grep 'Overall domain bounding box'`
z0=`echo $meshBounds | awk '{print $7}' | awk -F')' '{print $1}'`
zTop=`echo $meshBounds | awk '{print $10}' | awk -F')' '{print $1}'`
boundarySum b -case $case ground -time 0
boundarySum b -case $case top -time 0
b0=`tail -1 $case/sum_groundb.dat | awk '{print $5}'`
bTop=`tail -1 $case/sum_topb.dat | awk '{print $5}'`
H=`awk BEGIN'{print '$zTop'-'$z0'}'`
Uscale=`awk BEGIN'{print sqrt(('$b0'-('$bTop'))*'$H')}'`
KEscale=`awk BEGIN'{print '$Uscale'*'$nu'/'$H'}'`

# Globally average the KE and scale
globalSum -case $case -time 0':' KE
globalSum -case $case -time ${startAveraging}':' KETimeMean
awk '{print $1, $5/'$KEscale'}' $case/globalSumKE.dat \
    > $case/globalSumKEnorm.dat
awk '{print $1, $5/'$KEscale'}' $case/globalSumKETimeMean.dat \
    > $case/globalSumKETimeMeannorm.dat

# Plot global averaged and time mean of global average KE
gv=0
inputFiles=($case/globalSumKEnorm.dat $case/globalSumKETimeMeannorm.dat)
outFile=$case/globalSumKEnorm.eps
legends=("Global mean KE" "Time mean")
pens=("black" "1.5,black,3_3:0")
col=(2 2)
colx=1
xlabel="Time"
ylabel="Normalised kinetic energy (RE)"
xmin=0
xmax=`tail -1 $case/globalSumKEnorm.dat | awk '{print int($1+.999)}'`
dx=`echo $xmax | awk '{print int($1/10+.999)}'`
ymin=0
ymax=1400
dy=200
ddy=$dy
dyg=$ymax
xscale=/1
yscale=/1
nSkip=1
legPos=x5/6.6
source gmtPlot
ev $outFile

