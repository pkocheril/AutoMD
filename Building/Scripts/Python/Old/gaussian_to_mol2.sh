for file in *optfreq*.chk; do
    [ -f $file ] || continue
    formchk $file
done

for file in *optfreq*.fchk; do
    [ -f $file ] || continue
    obabel $file -O ${file%.fchk}.mol2
done