resultsFile=`pwd`/results.out

rm -rf $resultsFile > /dev/null 2>&1

dirs=`find . -maxdepth 1 -type d | grep -v '^\.$' | sed -e 's/^\.//' -e 's/^\///'`

runit()
{
  printf "running  %-25s:" $1 >> $resultsFile
  if test -x ./$1 ; then
    ./$1
    runitresult=$?
    if test $result -ne 0 ; then
      printf "Error\n" >> $resultsFile
    else
      printf "Passed\n" >> $resultsFile
    fi
  else
    printf "Not Run\n" >> $resultsFile
  fi
}


buildStatus=


for dir in $dirs
do
  cd $dir                      > /dev/null 2>&1
  make                         > /dev/null 2>&1
  result=$?
  printf "building  %-24s:" $dir >> $resultsFile
  if test $result -eq 0 ; then
    printf "Passed\n" >> $resultsFile
    runit $dir
    runit ${dir}-f
    runit ${dir}-mpi
    runit ${dir}-mpif90
  else
    printf "Failed\n" >> $resultsFile
  fi
  cd ..                        > /dev/null 2>&1
done

cat $resultsFile
rm $resultsFile
