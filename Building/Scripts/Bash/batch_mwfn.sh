for file in *.chk; do
    [ -f $file ] || continue
    formchk $file
done
    
for file in *SP*.fchk; do
    [ -f $file ] || continue
    echo -e "7\n18\n1\ny\n" | Multiwfn $file
done