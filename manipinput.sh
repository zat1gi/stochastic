#!/bin/bash
#-------------------------------------------------------------------------------
# author: Aaron Olson, 12-7-15, aolson1@unm.edu, aaronjeffreyolson@gmail.com
# function: bash script to easily exchange one string in a file for another
# original intent: easily and repetatively swap strings in input file to
#                  compare calculations with different parameters
# object of operation: currently 'inputstoc.txt', replace all for any other file
#-------------------------------------------------------------------------------

# test for valid number of inputs
if [ "$#" -ne 2 ] && [ "$#" -ne 3 ]; then
    echo
    echo "|------------------- *** Input Error *** --------------------|"
    echo "|------------------------------------------------------------|"
    echo "|               illegal number of parameters                 |"
    echo "| usage: ./manipinput [old string] [new string]              |"
    echo "|                             or                             |"
    echo "| usage: ./manipinput [linenumber] [old string] [new string] |"
    echo "|------------------------------------------------------------|"
    echo
    exit
fi

if [ "$#" -eq 2 ]; then  # if line number not specified, replace first occurance
  echo 
  echo "Replace first occurance of \"$1\" with \"$2\""
  echo
  echo "Replacement to be made here:"
  cat inputstoc.txt | grep $1 | head -n1
  line="$(grep -n "$1" inputstoc.txt | cut -c-1)"
  sed -i "s/$1/$2/" inputstoc.txt
  echo "Replacement line now reads:"
  sed "$line!d" inputstoc.txt
  echo
else                     # if line number specified, replace first occurance in line
  echo 
  echo "In line number \"$1\" replace \"$2\" with \"$3\""
  echo
  echo "  --Line before operation:"
  sed "$1q;d" inputstoc.txt
  sed -i "$1 s/$2/$3/" inputstoc.txt
  echo "  -- Line after operation:"
  sed "$1q;d" inputstoc.txt
  echo
fi
