for i in Ma01 Ma02 Ma03 Ma0357 Ma0364
do
cd $i
mpirun -np 4 python3 ../../usr/SU2/bin/fsi_computation.py --parallel -f fsi.cfg > log.txt &
cd ..
done
