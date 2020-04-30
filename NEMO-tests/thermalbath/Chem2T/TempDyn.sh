
rm -f TempDyn.dat
touch TempDyn.dat

for (( i=0; i<=$1; i+=10 )) 
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

  name='restart_flow_'$p$i'.csv'
  echo $name
  line=`sed -n 5p $name`  
  
  Tt=`echo $line | awk '{print $11}'`
  Tv=`echo $line | awk '{print $12}'`
  Dt=`echo $i $2 | awk '{print $1*$2}'`

  echo $Dt $Tt $Tv >> TempDyn.dat

done
