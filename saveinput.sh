#!/bin/bash
# Aaron Olson, 12-12-15
# Zero or one input.  Valid inputs: ls, ws, rm
# If zero inputs, store inputstoc.txt with date suffix
# If ls, list files in texts/inputlog
# If wc, count number of files in texts/inputlog
# If rm, remove files in textx/inputlog

if [ "$#" -eq 1 ]; then
  if [ $1 = "ls" ]; then
    ls texts/inputlog/ -I dummy.txt
  elif [ $1 = "wc" ]; then
    echo "$(($(ls texts/inputlog/ | wc -l) -1))"
  elif [ $1 = "rm" ]; then
    rm texts/inputlog/inputstoc*
  fi
elif [ "$#" -eq 0 ]; then
  suffix="$(date | sed 's/ /_/g')"
  cp inputstoc.txt texts/inputlog/inputstoc.txt_$suffix
fi
