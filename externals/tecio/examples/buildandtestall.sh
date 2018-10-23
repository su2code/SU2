dirs=`find . -maxdepth 1 -type d | grep -v '^\.$' | sed -e 's/^\.//' -e 's/^\///'`

for dir in $dirs
do
  printf "Building %-20s" "$dir"
  cd $dir                      > /dev/null 2>&1
  make                         > /dev/null 2>&1
  result=$?
  if test $result -eq 0 ; then
    echo "Success!"
  else
    echo "Error!"
  fi
  cd ..                        > /dev/null 2>&1
done
