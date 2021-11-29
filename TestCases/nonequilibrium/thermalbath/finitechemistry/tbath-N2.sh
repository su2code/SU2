rm -f $3
touch $3

echo 'time, rhoN2, rhoN, Tt, Tv' >> $3

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

  name='restart_flow_'$p$i'.csv'
  echo $name
  line=`sed -n 5p $name`  
   
  r1=`echo $line | awk '{print $4}'`
  r2=`echo $line | awk '{print $5}'`
  Tt=`echo $line | awk '{print $12}'`
  Tv=`echo $line | awk '{print $13}'`
  Dt=`echo $i $2 | awk '{print ($1+1)*$2}'`

  echo $Dt', '$r1', '$r2', '$Tt', '$Tv >> $3         

done
