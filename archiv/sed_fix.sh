#!/bin/bash

current=$(pwd)
trans=( '2v~2w' '2w~2z' '2y~2w' '2z~2y' 'ejm42~ejm31' 'ejm42~ejm55' 'ejm54~ejm42' 'ejm55~ejm54' )
# '67~60' '60~63' '63~61' '61~60' 
for dir in ${trans[@]} ; do
cd $dir
transdir=$(pwd)
bound_reps=( 'bound_0' 'bound_1' 'bound_2' )
for rep in ${bound_reps[@]} ; do
cd $rep
rdir=$(pwd)

lamvals=( 0.0000 0.1000 0.2000 0.3000 0.4000 0.5000 0.6000 0.7000 0.8000 0.9000 1.0000 )

for lam in ${lamvals[@]}; do

cd eq/lambda_$lam
mv somd.cfg old_somd.cfg
string=$(awk '/perturbed residue number/{print}' old_somd.cfg)
IFS=' ' read -a array <<< $string
res_number=${array[4]}
new_res_number="${res_number:1}"
new_string="perturbed residue number = $new_res_number"
sed -e "s/$string/$new_string/g" old_somd.cfg > somd.cfg
cd $rdir

cd lambda_$lam
mv somd.cfg old_somd.cfg
string=$(awk '/perturbed residue number/{print}' old_somd.cfg)
IFS=' ' read -a array <<< $string
res_number=${array[4]}
new_res_number="${res_number:1}"
new_string="perturbed residue number = $new_res_number"
sed -e "s/$string/$new_string/g" old_somd.cfg > somd.cfg
cd $rdir

done

cd $transdir

done

cd $current

done
