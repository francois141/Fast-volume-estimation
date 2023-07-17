#!/bin/sh

set -e

files=""

for ((i=1; i<=30; i++)); do
    echo $i
    perf stat -o ../perf_output/perf_original_cube_output_${i}.txt -e r5301c7 -e r5304c7 -e r5310c7 -e r5340c7 ../build/PolyVestCpp ../benchmark_infrastructure/cubes/cube_${i} 1600 > /dev/null
    files="${files} ../perf_output/perf_original_cube_output_${i}.txt"
done

# echo $files

python3 ../post-processing/proc_out_process.py $files