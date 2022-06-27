#!/bin/bash

# Script to record memory usage every interval

# For short jobs, set small, e.g 1m, 5, etc. 

# get PID of Rscript job(s)

while true; do
   ps -xm -v -p $(pgrep ^R$)  # Detect processed running "R"
   sleep 600
done
