
rm -f EDyn.dat
touch EDyn.dat

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
  
  E=`echo $line | awk '{print $8}'`
  Eve=`echo $line | awk '{print $9}'`
  Dt=`echo $i $2 | awk '{print $1*$2}'`

  Etr=`echo $E $Eve | awk '{print $1-$2}'`

  echo $Dt $Etr $Eve >> EDyn.dat

done
