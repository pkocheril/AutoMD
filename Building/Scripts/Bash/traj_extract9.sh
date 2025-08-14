for file in *md.trr; do
    [ -f $file ] || continue
    printf '2\n' | gmx traj -f $file -s ${file%.trr}.tpr -ox ${file%_md.trr}_xyz.xvg -of ${file%_md.trr}_f1.xvg
    printf '2\n' | gmx traj -f ${file%_md.trr}_0q.trr -s ${file%_md.trr}_0q.tpr -of ${file%_md.trr}_f0.xvg   
done
