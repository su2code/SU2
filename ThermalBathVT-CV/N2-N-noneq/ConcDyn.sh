
rm -f ConcDyn.dat
touch ConcDyn.dat

for (( i=0; i<=$1; i+=1 )) 
do

  if [[ $i -lt 10 ]]; then
    p='0000'
  elif [[ $i -lt 100 ]]; then
    p='000'
  elif [[ $i -lt 1000 ]]; then
    p='00'
  elif [[ $i -lt 10000 ]]; then
    p='0'
  else
    p=''
  fi

  name='restart_flow_'$p$i'.dat'
  echo $name
  line=`sed -n 5p $name`  
  
  N2=`echo $line | awk '{print $4}'`
  N=`echo $line | awk '{print $5}'`
  Dt=`echo $i $2 | awk '{print $1*$2}'`

  N2n=`echo $N2 60220000000000000000000 28 1.e19 | awk '{print ($1*$2/$3)/$4}'`
  Nn=`echo $N 60220000000000000000000 14 1.e19 | awk '{print ($1*$2/$3)/$4}'`

  N2nn=`echo $N2n $Nn | awk '{print($1)/($1+$2)}'`
  Nnn=`echo $N2n $Nn | awk '{print($2)/($1+$2)}'`

  echo $Dt $N2nn $Nnn >> ConcDyn.dat

done
