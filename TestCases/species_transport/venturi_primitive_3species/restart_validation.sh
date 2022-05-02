# ------------------------------------------------------------------------------------------------ #
# T. Kattmann, 28.09.2021, restart_validation.sh
#
# This script validates that the steady restarts work perfectly in SU2.
# For the primal-only restart and the primal restart for the adjoint.
# ------------------------------------------------------------------------------------------------ #

# Prepare useful variable-strings that run primal and adjoint code
code_dir="/home/kat7rng/1__SU2-code/1__Vanilla-SU2/bin"
num_cores=1
config="species3_primitiveVenturi.cfg"
adapt_config="adapt_fluid.cfg"
# Iteration number at which to compare restart with full iter. Note this value is the actual screen iteration number, so for FullIter we have to set ITER=compare_iter+1.
compare_iter=10
echo "ITER: $compare_iter"
compare_iter_plus1=$((compare_iter + 1))

# ITER for SingleZone-steady, OUTER_ITER for multizone-steady
iter_string="ITER"

run_command_primal="  mpirun -n $num_cores $code_dir/SU2_CFD    $adapt_config"
run_command_adjoint=" mpirun -n $num_cores $code_dir/SU2_CFD_AD $adapt_config"

echo $run_command_primal
echo $run_command_adjoint

# Extract MESH_FILENAME= name from the provided $config . This is necessary in order to link it into the folders
# The "^" in the grep signals the beginning of the line so that commented out lines are disregarded ( https://unix.stackexchange.com/questions/60994/how-to-grep-lines-which-does-not-begin-with-or )
# The cut command extract everyhing between the first and next "="-char ( https://stackoverflow.com/questions/15148796/get-string-after-character )
# The xargs command trims any unnecessary whitespaces from the string ( https://stackoverflow.com/questions/369758/how-to-trim-whitespace-from-a-bash-variable )
mesh=$(grep "^MESH_FILENAME=" $config | cut -d "=" -f2 | xargs)
echo $mesh

# set folder names
FullIter_dir="1__FullIter"
FullMinus1Iter_dir="2__FullMinus1Iter"
PrimalRestart_dir="3__PrimalRestart"
AdjointRestart_dir="4__AdjointRestart"

echo "Deleting old folders for a fresh test."
rm -r $FullIter_dir $FullMinus1Iter_dir $PrimalRestart_dir $AdjointRestart_dir

link_config="ln -s ../*.cfg ."
echo $link_config
link_mesh="ln -s ../$mesh"
echo $link_mesh

link_restart0="ln -s ../$FullMinus1Iter_dir/restart.csv solution.csv"
link_restart1="ln -s ../$FullMinus1Iter_dir/restart.csv solution.csv"

# ------------------------------------------------------------------------------------------------ #
# Run FullIter and FullMinus1Iter in parallel

mkdir $FullIter_dir
cd $FullIter_dir
$link_config
$link_mesh
# Change ITER number. The ^ ensures that the string start at the beginning of the line, such that other *ITER= are not changed
sed "s/^$iter_string=.*/$iter_string=$compare_iter_plus1/g" $config > $adapt_config
echo "Running full primal."
$run_command_primal > CFD.log &
cd ..

mkdir $FullMinus1Iter_dir
cd $FullMinus1Iter_dir
$link_config
$link_mesh
# Change ITER number
sed "s/^$iter_string=.*/$iter_string=$compare_iter/g" $config > $adapt_config
echo "Running full-1 primal."
$run_command_primal > CFD.log
cd ..

# ------------------------------------------------------------------------------------------------ #
# Run primal and adjoint restart in parallel

mkdir $PrimalRestart_dir
cd $PrimalRestart_dir
$link_config
$link_mesh
$link_restart0
$link_restart1
# Change ITER number to 1
sed "s/^$iter_string=.*/$iter_string=1/g" $config > $adapt_config
sed -i "s/RESTART_SOL=.*/RESTART_SOL= YES/g"  $adapt_config
#sed -i "s/RESTART_ITER=.*/RESTART_ITER= $compare_iter/g"  $adapt_config # only applies to unsteady simulations
echo "Running restarted primal."
$run_command_primal > CFD.log &
cd ..

mkdir $AdjointRestart_dir
cd $AdjointRestart_dir
$link_config
$link_mesh
$link_restart0
$link_restart1
# Change ITER number
sed "s/^$iter_string=.*/$iter_string=1/g" $config > $adapt_config
echo "Running adjoint."
$run_command_adjoint > CFD_AD.log
cd ..


# ------------------------------------------------------------------------------------------------ #
# Postproces results: Print residuals to screen. Remove delimiter using sed and the xargs cuts unnecessary spaces
delimiter=","
remove_delimiter="sed "s/$delimiter//g""
history_file="history.csv"
tail -n 1 $FullIter_dir/$history_file | $remove_delimiter | xargs > output.txt
tail -n 1 $PrimalRestart_dir/$history_file | $remove_delimiter | xargs >> output.txt
delimiter="|"
remove_delimiter="sed "s/$delimiter//g""
# retrieve line number
adj_lineNumber=$(sed -n '/rms_Flow/=' $AdjointRestart_dir/CFD_AD.log)
# $((..)) necessary to expand number
sed -n $((adj_lineNumber + 2))p $AdjointRestart_dir/CFD_AD.log | $remove_delimiter | xargs >> output.txt
cat output.txt
