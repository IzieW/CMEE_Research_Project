# !/bin/bash
#PBS -lwalltime=00:20:00
#PBS -lselect=1:ncpus=1:mem=1gb
# Run simulations on the cluster, move results to home

module load anaconda3/personal

cp $HOME/lenia_package.py .
cp $HOME/multiple_creatures_nutrient_B.py .
cp $HOME/orbium_parameters.csv .
cp -r $HOME/multiple_B .

echo "Script is about to run"

python3 $HOME/render_test_multiple_B.py

echo "Simulation complete, moving to Home..."

mv * $HOME

echo "Program has finished running"
