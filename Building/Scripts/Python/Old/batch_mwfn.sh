for file in *.fchk; do
    [ -f $file ] || continue
    echo -e "7\n18\n1\ny\n" | Multiwfn $file
done