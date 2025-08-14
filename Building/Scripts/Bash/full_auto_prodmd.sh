for file in *MIX.gro; do
    [ -f $file ] || continue
    gmx grompp -f md/em.mdp -c $file -p topol.top -o ${file%.gro}_em.tpr -maxwarn 0
done

for file in *em.tpr; do
    [ -f $file ] || continue
    gmx mdrun -v -deffnm ${file%.tpr} -nt ${nproc}
done

for file in *_em.gro; do
    [ -f $file ] || continue
    gmx grompp -f md/v.mdp -c $file -p topol.top -o ${file%_em.gro}_v.tpr -maxwarn 0
done

for file in *v.tpr; do
    [ -f $file ] || continue
    gmx mdrun -v -deffnm ${file%.tpr} -nt ${nproc}
done

for file in *v.gro; do
    [ -f $file ] || continue
    gmx grompp -f md/p.mdp -c $file -t ${file%.gro}.cpt -p topol.top -o ${file%_v.gro}_p.tpr -maxwarn 0
done

for file in *p.tpr; do
    [ -f $file ] || continue
    gmx mdrun -v -deffnm ${file%.tpr} -nt ${nproc}
done



for file in *p.gro; do
    [ -f $file ] || continue
    gmx grompp -f md/md.mdp -c $file -t ${file%.gro}.cpt -p topol.top -o ${file%_p.gro}_md.tpr -maxwarn 0
done

for file in *p.gro; do
    [ -f $file ] || continue
    gmx grompp -f md/md.mdp -c $file -t ${file%.gro}.cpt -p topol_0q.top -o ${file%_p.gro}_0q.tpr -maxwarn 1
done

for file in *md.tpr; do
    [ -f $file ] || continue
    gmx mdrun -v -deffnm ${file%.tpr} -nt ${nproc}
done

for file in *md.trr; do
    [ -f $file ] || continue
    gmx mdrun -v -deffnm ${file%_md.trr}_0q -rerun ${file%.trr}.trr -nt ${nproc}
done
