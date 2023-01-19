#! /bin/bash

echo "running the script..."

export net_file="/backup/anna/benchmark/tyk2/execution_model/network_combined.dat"

# Make a list of the transformations from the network.dat file
trans_array=()
eng_array=()
win_array=()
IFS=' '
while read trans; do
while read -a tra; do tran=${tra[0]}~${tra[1]}; eng=${tra[-1]}; win=${tra[2]};
trans_array+=("$tran"); eng_array+=("$eng"); win_array+=("$win"); done <<< $trans
done < $net_file

repeats_array=( 0 1 2 )

echo ${trans_array[@]}

output_dir="/backup/anna/benchmark/tyk2/outputs_extracted/AMBER"

lamvals=( 0.0000 0.1000 0.2000 0.3000 0.4000 0.5000 0.6000 0.7000 0.8000 0.9000 1.0000 )

for pert in ${trans_array[@]}; do

for dir in 'bound' 'free'; do
for rep in ${repeats_array[@]} ; do

repeat_dir=${dir}_${rep}

echo "$output_dir/$pert/$repeat_dir"

for lam in ${lamvals[@]}; do

correct_folder_name="$output_dir/$pert/$repeat_dir/lambda_$lam"
folder_name="$output_dir/$pert/$repeat_dir/$repeat_dir/lambda_$lam"

if [ -s $correct_folder_name/amber.out ]; then
echo "file is in the correct folder already"
else
echo "copying"
cp $folder_name/amber.out $correct_folder_name/amber.out
fi

done

done
done

done
