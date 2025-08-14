mkdir 00_inputs

for file in *.*; do
    [ -f $file ] || continue
    mv $file 00_inputs/
done

mkdir 01_convert
mkdir 02_geometry
mkdir 03_gaussian
mkdir 04_resp
mkdir 05_solute_params
mkdir 06_solvent_md



