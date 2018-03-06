#!/bin/bash

CPPFILES=../../*/src/*.cpp
HPPFILES=../../*/include/*.hpp
INLFILES=../../*/include/*.inl

# First fix the CPP files 
for f in $CPPFILES
do
    TESTFILE=$(file "$f")
    if [[ $TESTFILE == *"with CRLF line terminators"*  ]]; then
      echo $TESTFILE
      echo "Fixing "$f"..."
      dos2unix $f
    fi
done

# Now fix the HPP files
for f in $HPPFILES
do
    TESTFILE=$(file "$f")
    if [[ $TESTFILE == *"with CRLF line terminators"*  ]]; then
      echo $TESTFILE
      echo "Fixing "$f"..."
      dos2unix $f
    fi
done

# Finally fix the INL files
for f in $INLFILES
do
    TESTFILE=$(file "$f")
    if [[ $TESTFILE == *"with CRLF line terminators"*  ]]; then
      echo $TESTFILE
      echo "Fixing "$f"..."
      dos2unix $f
    fi
done
