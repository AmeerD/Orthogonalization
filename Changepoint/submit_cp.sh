#!/bin/bash

njobs=100

sbatch -p witten-12c128g -t 96:00:00 -e ./iotrash/t-%A_%a.out -o ./iotrash/t-%A_%a.out --array=1-$njobs ./CP.sh $1 $2
