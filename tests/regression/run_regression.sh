#/bin/bash

source /usr/workspace/coldqcd/software/lassen_smpi_RR/install/env.sh
EXE=/usr/workspace/coldqcd/software/lassen_smpi_RR/install/lalibe_regression/bin/lalibe
lalibe_tests=/usr/workspace/coldqcd/software/src/lalibe_regression/tests/regression
input_decks=$lalibe_tests/input_decks
my_python=/usr/workspace/coldqcd/software/python_venv-3.7.2.lassen/bin/python

mpirun="jsrun -n 1 $EXE"

new_dir=/p/gpfs1/walkloud/c51/x_files/project_2/test_lalibe/regression
pushd $new_dir
rm *

ini_files=(
    source_prop_h5.ini.xml
    mesons_baryons_h5.ini.xml
    fh-props_h5.ini.xml
    fh-corrs_h5.ini.xml
    proton_seqprop_h5.ini.xml
    proton_formfac_h5.ini.xml
)

for ini in "${ini_files[@]}"; do
    $mpirun -i $input_decks/$ini
done

$my_python $lalibe_tests/py_scripts/perform_revision_test.py $lalibe_tests/known_results $new_dir

popd
