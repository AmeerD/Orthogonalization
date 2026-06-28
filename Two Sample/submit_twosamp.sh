#!/bin/bash

njobs=500

sbatch -p students-12c128g -t 96:00:00 -e ./iotrash/ts-%A_%a.out -o ./iotrash/ts-%A_%a.out --array=1-$njobs ./twosamp.sh $1 $2
