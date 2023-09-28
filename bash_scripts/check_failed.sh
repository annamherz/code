#! /bin/bash

slurmlog=$1
slurm_dir=$2

echo checking $slurm_logs in $slurm_dir

file_base="$slurm_dir/prod_$slurmlog\_*"

declare regex="\s+STOP PMEMD Terminated Abnormally!\s+"

for f in $file_base; do

declare file_content=$( cat "${f}" )
if [[ " $file_content " =~ $regex ]] # please note the space before and after the file content
    then
        echo "$f is a failed run"
fi

done
