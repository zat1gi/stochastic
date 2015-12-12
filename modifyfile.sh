#!/bin/bash
#-------------------------------------------------------------------------------
# author: Aaron Olson, 12-10-15, aolson1@unm.edu, aaronjeffreyolson@gmail.com
# function: bash script to easily exchange one string in a file for another
# input format: see *** Input Error *** below
#-------------------------------------------------------------------------------

# test for valid number of inputs
if [ "$#" -ne 3 ] && [ "$#" -ne 4 ]; then
    echo
    echo "|------------------------- *** Input Error *** ---------------------------|"
    echo "|-------------------------------------------------------------------------|"
    echo "|                     illegal number of parameters                        |"
    echo "| usage: ./manipinput [old string] [new string] [file  path]              |"
    echo "|                             or                                          |"
    echo "| usage: ./manipinput [linenumber] [old string] [new string] [file  path] |"
    echo "|-------------------------------------------------------------------------|"
    echo
    exit
fi

if [ "$#" -eq 3 ]; then  # if line number not specified, replace first occurance
  file=$3
  echo $file
  echo 
  echo "| Replace first occurance of \"$1\" with \"$2\""
  echo "| "
  if [ ! -e $file ]; then  # test that file exists
    echo "| --- File \"$file\" does not exist --- "
    echo
    exit
  fi
  numoccurances="$(cat $file | grep $1 | wc -l)"
  if [ $numoccurances -eq 0 ]; then #test if old string exists in file, throw error if not
    echo "| --- String \"$1\" does not exist in \"$file\" --- "
    echo
    exit
  fi
  echo "| Replacement to be made here:"
  echo "| $(cat $file | grep $1 | head -n1)"
  line="$(grep -n "$1" $file | cut -f -d":")"
  sed -i "s/$1/$2/" $file
  echo "| Replacement line now reads:"
  echo "| $(sed "$line!d" $file)"
  echo
else                     # if line number specified, replace first occurance in line
  file=$4
  echo $file
  echo 
  echo "| In line number \"$1\" replace \"$2\" with \"$3\""
  echo "| "
  if [ ! -e $file ]; then  # test that file exists
    echo "| --- File \"$file\" does not exist --- "
    echo
    exit
  fi
  numlines="$(cat $file | wc -l)"
  if [ $numlines -lt $1 ]; then  # test that line specified exists in file
    echo "|  --- \"$file\" only has $numlines lines, less than user input of \"$1\" ---"
    echo
    exit
  fi
  echo "|   --Line before operation:"
  echo "| $(sed "$1q;d" $file)"
  sed -i "$1 s/$2/$3/" $file
  echo "|   -- Line after operation:"
  echo "| $(sed "$1q;d" $file)"
  echo
fi
