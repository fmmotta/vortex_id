cp ../0_template/* 0/
blockMesh
checkMesh
t0=10.0075
dt=0.0120
n0=0
nf=5
for n in `seq $n0 $nf`; 
do
  t=`awk "BEGIN { print $t0+$n*$dt}";` 
  echo "Maping fields for time $t" 
  mkdir $t
  mapFields ../DNS_OPEN_FOAM/ -sourceTime $t
  cp -r 0/* $t 
  echo "Computed and copied the results for time $t"
done