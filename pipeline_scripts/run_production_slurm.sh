#!/bin/bash
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=prod
#SBATCH --time=24:00:00
#SBATCH -o ../slurm_logs/prod_%A_%a.out
#SBATCH -e ../slurm_logs/prod_%A_%a.err

# specifying the number of threads (needed for gromacs)
export OMP_NUM_THREADS=8

# sourcing
source $BSS
source $amber/amber.sh
source $gromacs
source $scripts_dir/extract_execution_model_bash.sh

date
start=`date +%s`

echo "Folder for these runs is : $MAINDIRECTORY"
echo "The transformation is $1 using $3 windows and $2 as the MD engine with $4 repeats"

# define no of windows based on the assosciative array generated from reading the net_file in extract execution model
IFS=","
read -r -a lamvals <<< "${wins_array[$3]}"

# define lambda based on slurm array
lam=${lamvals[SLURM_ARRAY_TASK_ID]}

# change to the trans dir, abort and message if not there
cd $MAINDIRECTORY/outputs/$2/$1
if [[ ! -d $MAINDIRECTORY/outputs/$2/$1 ]]; then
    echo "$MAINDIRECTORY/outputs/$2/$1 does not exist. Production run aborted..."
    exit
fi
trans_dir=$(pwd)

# iterate over dir (for each leg) based on no of repeats
for dir in 'bound' 'free'; do
for rep in "${repeats_array[@]}" ; do
repeat=${dir}_${rep}
cd $repeat
repeat_dir=$(pwd)

echo "Running in $repeat_dir"
echo "Lambda is $lam"

# run the runs based on which engine
if [ $2 = "AMBER" ]; then
echo "min"
    cp min/lambda_$lam/amber.rst7 min/lambda_$lam/initial_amber.rst7
    pmemd.cuda -i min/lambda_$lam/amber.cfg -c min/lambda_$lam/amber.rst7 -ref min/lambda_$lam/amber.rst7 -p min/lambda_$lam/amber.prm7 -O -o min/lambda_$lam/amber.out -inf min/lambda_$lam/amber.info -e min/lambda_$lam/amber.en -r min/lambda_$lam/amber.rst7 -x min/lambda_$lam/amber.nc -l min/lambda_$lam/amber.log ;
echo "heat"
    pmemd.cuda -i heat/lambda_$lam/amber.cfg -c min/lambda_$lam/amber.rst7 -ref min/lambda_$lam/amber.rst7 -p heat/lambda_$lam/amber.prm7 -O -o heat/lambda_$lam/amber.out -inf heat/lambda_$lam/amber.info -e heat/lambda_$lam/amber.en -r heat/lambda_$lam/amber.rst7 -x heat/lambda_$lam/amber.nc -l heat/lambda_$lam/amber.log ;
echo "eq"
    pmemd.cuda -i eq/lambda_$lam/amber.cfg -c heat/lambda_$lam/amber.rst7 -ref heat/lambda_$lam/amber.rst7 -p eq/lambda_$lam/amber.prm7 -O -o eq/lambda_$lam/amber.out -inf eq/lambda_$lam/amber.info -e eq/lambda_$lam/amber.en -r eq/lambda_$lam/amber.rst7 -x eq/lambda_$lam/amber.nc -l eq/lambda_$lam/amber.log ;
echo "prod"
    pmemd.cuda -i lambda_$lam/amber.cfg -c eq/lambda_$lam/amber.rst7 -ref eq/lambda_$lam/amber.rst7 -p lambda_$lam/amber.prm7 -O -o lambda_$lam/amber.out -inf lambda_$lam/amber.info -e lambda_$lam/amber.en -r lambda_$lam/amber.rst7 -x lambda_$lam/amber.nc -l lambda_$lam/amber.log ;

# delete simulation data 
if [[ $keep_traj == "None" ]]; then
    echo "Removing all trajectories..."
    rm min/lambda_$lam/*.nc heat/lambda_$lam/*.nc eq/lambda_$lam/*.nc
    rm lambda_$lam/*.nc
fi

fi

if [ $2 = "GROMACS" ]; then

echo "min"
cp min/lambda_$lam/gromacs.gro min/lambda_$lam/initial_gromacs.gro
gmx grompp -f min/lambda_$lam/gromacs.mdp -c min/lambda_$lam/gromacs.gro -p min/lambda_$lam/gromacs.top -o min/lambda_$lam/gromacs.tpr
gmx mdrun -ntmpi 1 -deffnm min/lambda_$lam/gromacs ;

echo "min1"
gmx grompp -f min1/lambda_$lam/gromacs.mdp -c min/lambda_$lam/gromacs.gro -p min1/lambda_$lam/gromacs.top -o min1/lambda_$lam/gromacs.tpr
gmx mdrun -ntmpi 1 -deffnm min1/lambda_$lam/gromacs ;

echo "min2"
gmx grompp -f min2/lambda_$lam/gromacs.mdp -c min1/lambda_$lam/gromacs.gro -p min2/lambda_$lam/gromacs.top -o min2/lambda_$lam/gromacs.tpr
gmx mdrun -ntmpi 1 -deffnm min2/lambda_$lam/gromacs ;

echo "heat"
gmx grompp -f heat/lambda_$lam/gromacs.mdp -c min2/lambda_$lam/gromacs.gro -p heat/lambda_$lam/gromacs.top -o heat/lambda_$lam/gromacs.tpr
gmx mdrun -ntmpi 1 -deffnm heat/lambda_$lam/gromacs ;

echo "eq"
gmx grompp -f eq/lambda_$lam/gromacs.mdp -c heat/lambda_$lam/gromacs.gro -p eq/lambda_$lam/gromacs.top -t heat/lambda_$lam/gromacs.cpt  -o eq/lambda_$lam/gromacs.tpr
gmx mdrun -ntmpi 1 -deffnm eq/lambda_$lam/gromacs ;
echo "prod"
gmx grompp -f lambda_$lam/gromacs.mdp -c eq/lambda_$lam/gromacs.gro -p lambda_$lam/gromacs.top -t eq/lambda_$lam/gromacs.cpt -o lambda_$lam/gromacs.tpr
gmx mdrun -ntmpi 1 -deffnm lambda_$lam/gromacs ;

# delete simulation data 
if [[ $keep_traj == "None" ]]; then
    echo "Removing all trajectories..."
    rm min/lambda_$lam/*.trr heat/lambda_$lam/*.trr eq/lambda_$lam/*trr
    rm lambda_$lam/*.trr
fi

fi

if [ $2 = "SOMD" ]; then
echo "min + eq"
cd eq/lambda_$lam
somd-freenrg -c somd.rst7 -t somd.prm7 -m somd.pert -C somd.cfg -p CUDA

cd $repeat_dir
#copy the files over for restart
echo "copying restart files for production..."
cp eq/lambda_$lam/sim_restart.s3 lambda_$lam/sim_restart.s3
cp eq/lambda_$lam/sim_restart.s3.previous lambda_$lam/sim_restart.s3.previous
cp eq/lambda_$lam/SYSTEM.s3 lambda_$lam/SYSTEM.s3

echo "prod"
cd lambda_$lam
somd-freenrg -c somd.rst7 -t somd.prm7 -m somd.pert -C somd.cfg -p CUDA
cd $repeat_dir

# delete simulation data 
if [[ $keep_traj == "None" ]]; then
    echo "Removing all trajectories..."
    rm eq/lambda_$lam/*.dcd eq/lambda_$lam/*.s3* eq/lambda_$lam/latest*
    rm lambda_$lam/*.dcd lambda_$lam/*.s3* lambda_$lam/latest*
fi

fi

cd $trans_dir

done
done

echo "done."

end=`date +%s`
runtime=$((end-start))
echo "runtime was $runtime seconds, which is $((runtime/60)) minutes"


# delete trajectory files

#removing all traj
# remove amber trajectories
if [ "$2" = "AMBER" ] || [ "$2" = "ALL" ] && [ "$keep_traj" = "None" ]; then
for i in $(find $trans_dir -name 'amber.nc');
do
    rm -rf $i
done
fi

# remove gromacs trajectories
if [ "$2" = "GROMACS" ] || [ "$2" = "ALL" ] && [ "$keep_traj" = "None" ]; then
for i in $(find $trans_dir -name 'gromacs.trr');
do
    rm -rf $i
done
fi

# remove somd trajectories
if [ "$2" = "SOMD" ] || [ "$2" = "ALL" ] && [ "$keep_traj" = "None" ]; then
for i in $(find $trans_dir -name 'traj*.dcd');
do
    rm -rf $i
done
for i in $(find $trans_dir -name '*.s3*');
do
    rm -rf $i
done
fi

#removing all except lambda 0 and 1
# remove amber trajectories
if [ "$2" = "AMBER" ] || [ "$2" = "ALL" ] && [ "$keep_traj" = "0,1" ]; then
for i in $(find $trans_dir -name 'amber.nc');

do
if [[ "$i" != *"lambda_0.0000"* ]] && [[ "$i" != *"lambda_1.0000"* ]]; then
    rm -rf $i
fi
done

fi

# remove gromacs trajectories
if [ "$2" = "GROMACS" ] || [ "$2" = "ALL" ]&& [ "$keep_traj" = "0,1" ]; then
for i in $(find $trans_dir -name 'gromacs.trr');

do
if [[ "$i" != *"lambda_0.0000"* ]] && [[ "$i" != *"lambda_1.0000"* ]]; then
    rm -rf $i
fi
done

fi

# remove somd trajectories
if [ "$2" = "SOMD" ] || [ "$2" = "ALL" ]&& [ "$keep_traj" = "0,1" ]; then
for i in $(find $trans_dir -name 'traj*.dcd');

do
if [[ "$i" != *"lambda_0.0000"* ]] && [[ "$i" != *"lambda_1.0000"* ]]; then
    rm -rf $i
fi
done

for i in $(find $trans_dir -name '*.s3*');

do
if [[ "$i" != *"lambda_0.0000"* ]] && [[ "$i" != *"lambda_1.0000"* ]]; then
    rm -rf $i
fi
done

fi

#removing all except lambda 0, 0.5 and 1
# remove amber trajectories
if [ "$2" = "AMBER" ] || [ "$2" = "ALL" ] && [ "$keep_traj" = "0,0.5,1" ]; then
for i in $(find $trans_dir -name 'amber.nc');

do
if [[ "$i" != *"lambda_0.0000"* ]] && [[ "$i" != *"lambda_0.5000"* ]] && [[ "$i" != *"lambda_1.0000"* ]]; then
    rm -rf $i
fi
done

fi

# remove gromacs trajectories
if [ "$2" = "GROMACS" ] || [ "$2" = "ALL" ]&& [ "$keep_traj" = "0,0.5,1" ]; then
for i in $(find $trans_dir -name 'gromacs.trr');

do
if [[ "$i" != *"lambda_0.0000"* ]] && [[ "$i" != *"lambda_0.5000"* ]] && [[ "$i" != *"lambda_1.0000"* ]]; then
    rm -rf $i
fi
done

fi

# remove somd trajectories
if [ "$2" = "SOMD" ] || [ "$2" = "ALL" ]&& [ "$keep_traj" = "0,0.5,1" ]; then
for i in $(find $trans_dir -name 'traj*.dcd');

do
if [[ "$i" != *"lambda_0.0000"* ]] && [[ "$i" != *"lambda_0.5000"* ]] && [[ "$i" != *"lambda_1.0000"* ]]; then
    rm -rf $i
fi
done

for i in $(find $trans_dir -name '*.s3*');

do
if [[ "$i" != *"lambda_0.0000"* ]] && [[ "$i" != *"lambda_0.5000"* ]] && [[ "$i" != *"lambda_1.0000"* ]]; then
    rm -rf $i
fi
done

fi

echo "Done deleting, keeping $keep_traj trajectories."
