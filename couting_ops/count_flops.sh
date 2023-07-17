#!/bin/sh

set -e

# for executable in "../build/PolyVestBasicOpt" "../build/PolyVestFastest"; do
# for executable in "../build/PolyVestBaseline"; do
for executable in "../build/PolyVestAlgoOpt"; do
    echo $executable
    # take text from last slash to the end of executable
    executable_name=${executable##*/}
    echo $executable_name

    files=""

    for ((i=3; i<=50; i+=2)); do
        echo $i
        perf stat -o ../perf_output/perf_${executable_name}_cube_output_${i}.txt -e r5301c7 -e r5304c7 -e r5310c7 -e r5340c7 ${executable} ../benchmark_infrastructure/cubes/cube_${i} 1600 > /dev/null
        files="${files} ../perf_output/perf_${executable_name}_cube_output_${i}.txt"
    done

    echo $files

    python3 ../post-processing/proc_out_process.py $files
    mv res.txt ${executable_name}_res.txt
done

# files=""

# for ((i=1; i<=30; i++)); do
#     echo $i
#     perf stat -o ../perf_output/perf_original_cube_output_${i}.txt -e r5301c7 -e r5304c7 -e r5310c7 -e r5340c7 ../build/PolyVestCpp ../benchmark_infrastructure/cubes/cube_${i} 1600 > /dev/null
#     files="${files} ../perf_output/perf_original_cube_output_${i}.txt"
# done

# # echo $files

# python3 ../post-processing/proc_out_process.py $files