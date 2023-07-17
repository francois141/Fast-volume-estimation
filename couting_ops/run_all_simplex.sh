#!/bin/sh

set -e

files=""

for ((i=5; i<=50; i+=5)); do
    echo $i
    # perf stat -o ../perf_output/perf_simplex_output_${i}.txt -e r5301c7 -e r5304c7 -e r5310c7 -e r5340c7 ../build/PolyVestFast ../examples/simplex/simplex_${i} 1600 > /dev/null
    perf stat -o ../perf_output/perf_original_simplex_output_${i}.txt -e r5301c7 -e r5304c7 -e r5310c7 -e r5340c7 ../build/PolyVestCpp ../examples/simplex/simplex_${i} 1600 > /dev/null
    # files="${files} ../perf_output/perf_simplex_output_${i}.txt"
    files="${files} ../perf_output/perf_original_simplex_output_${i}.txt"
done

# echo $files

python3 ../post-processing/proc_out_process.py $files