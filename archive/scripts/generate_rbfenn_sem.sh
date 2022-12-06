#! /bin/bash 
#

export BSS="/home/anna/anaconda3/bin/activate rbfenn"
source $BSS

prots=( "tyk2" "mcl1" "p38" )

cd RBFENN/ANALYSIS/perturbation_networks

for prot in ${prots[@]}; do
echo $prot
python _04_external_test_all_series.py $prot
echo "done."
done