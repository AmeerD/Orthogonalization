#!/bin/bash

njobs=500

sbatch -p witten-12c128g -t 96:00:00 -e ./iotrash/db-%A_%a.out -o ./iotrash/db-%A_%a.out --array=1-$njobs ./zipclust.sh $1 $2 $3
