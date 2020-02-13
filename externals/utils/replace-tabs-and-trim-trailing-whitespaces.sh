#! /bin/bash

# has to be called from code-root

for folder in SU2_AD \
              SU2_BASE \
              SU2_CFD \
              SU2_DEF \
              SU2_DIRECTDIFF \
              SU2_DOT \
              SU2_GEO \
              SU2_MSH \
              SU2_SOL
do
  cd $folder
  echo $folder

  # replace all leading tabs by 2 spaces, for some reason trailing tabs and tabs enclosed by characters are only replaced by 2-1=1 spaces.
  # The trailing tabs are irrelevant anyway as they are trimmed and enclosed tabs require new formatting anyway

  #https://stackoverflow.com/questions/11094383/how-can-i-convert-tabs-to-spaces-in-every-file-of-a-directory
  #https://unix.stackexchange.com/questions/15308/how-to-use-find-command-to-search-for-multiple-extensions

  echo -n "Replace all tabs by 2 spaces..."
  find . \( -name 'C*.cpp' -o -name 'C*.hpp' -o -name 'C*.inl' \) ! -type d -exec bash -c 'expand -t 2 "$0" > /tmp/e && mv /tmp/e "$0"' {} \;
  #find . \( -name 'C*.cpp' -o -name 'C*.hpp' -o -name 'C*.inl' \) 
  echo "done"

  # deletes trailing whitespaces in all cpp,hpp,inl,build,py files

  #https://stackoverflow.com/questions/10711051/how-to-remove-trailing-whitespaces-for-multiple-files
  #https://unix.stackexchange.com/questions/15308/how-to-use-find-command-to-search-for-multiple-extensions

  echo -n "Trim all trailing whitespaces..."
  find . -type f \( -name 'C*.cpp' -o -name 'C*.hpp' -o -name 'C*.inl' \) -exec sed --in-place 's/[[:space:]]\+$//' {} \+
  echo "done"

  cd ..
done
