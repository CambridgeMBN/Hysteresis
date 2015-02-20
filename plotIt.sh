#!/bin/bash
FILES=Arrow_Data/*
for f in $FILES[0]
do
  echo "Processing $f file..."
  # take action on each file. $f store current file name
  python Arrow_Plotter.py $f
  # cat $f
done