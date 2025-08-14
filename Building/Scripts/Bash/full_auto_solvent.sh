#for file in *box*.gro; do
#    [ -f $file ] || continue
#    rm $file
#done

#for file in "#*box*#"; do
#    [ -f $file ] || continue
#    rm $file
#done

#for file in "#*solv*#"; do
#    [ -f $file ] || continue
#    rm $file
#done


for file in solv*.gro; do
    [ -f $file ] || continue
    gmx grompp -f md/em.mdp -c $file -p solv.top -o ${file%.gro}_em.tpr -maxwarn 0
done

for file in *em.tpr; do
    [ -f $file ] || continue
    gmx mdrun -v -deffnm ${file%.tpr} -nt ${nproc}
done

for file in *_em.gro; do
    [ -f $file ] || continue
    gmx grompp -f md/v.mdp -c $file -p solv.top -o ${file%_em.gro}_v.tpr -maxwarn 1
done

for file in *v.tpr; do
    [ -f $file ] || continue
    gmx mdrun -v -deffnm ${file%.tpr} -cpo solv_v.cpt -nt ${nproc}
done

for file in *v.gro; do
    [ -f $file ] || continue
    gmx grompp -f md/p.mdp -c $file -t ${file%.gro}_v.cpt -p solv.top -o ${file%_v.gro}_p.tpr -maxwarn 0
done

for file in *p.tpr; do
    [ -f $file ] || continue
    gmx mdrun -v -deffnm ${file%.tpr} -nt ${nproc}
done



for file in *v.edr; do
    [ -f $file ] || continue
    echo -e "14\n\n" | gmx energy -f $file -o temperature.xvg
done

for file in *p.edr; do
    [ -f $file ] || continue
    echo -e "22\n\n" | gmx energy -f $file -o density.xvg
done

