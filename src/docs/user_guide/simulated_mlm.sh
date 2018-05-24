#!/usr/bin/env bash

for $filename in $@;
do
   echo "Running $filename now.\n"
   ./run_pipeline.pl -Xmx6g -configFile $filename
done
