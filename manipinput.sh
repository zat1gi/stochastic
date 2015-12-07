#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo
    echo "|------------------- *** Input Error *** --------------------|"
    echo "|------------------------------------------------------------|"
    echo "|               illegal number of parameters                 |"
    echo "| usage: ./manipinput [linenumber] [old string] [new string] |"
    echo "|------------------------------------------------------------|"
    echo
    exit
fi
echo 
echo "In line number \"$1\" replace \"$2\" with \"$3\""
echo
echo "  --Line before operation:"
sed "$1q;d" inputstoc.txt
sed -i "$1 s/$2/$3/" inputstoc.txt
echo "  -- Line after operation:"
sed "$1q;d" inputstoc.txt
echo
