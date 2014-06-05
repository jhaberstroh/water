#!/bin/bash

MDP_DIR="/home/jhaberstroh/Code/gromacs-scripts/mdp/water"


#mkdir em
#cd em
#grompp -c ../waterbox.gro -f $MDP_DIR/water_em.mdp -p ../topol.top -o em.tpr
#mdrun -deffnm em
#cd ..
#
#mkdir nvt
#cd nvt
#grompp -c ../em/em.gro -f $MDP_DIR/water_nvt.mdp -p ../topol.top -o nvt.tpr
#mdrun -deffnm nvt
#cd ..
#
#
#mkdir npt
#cd npt
#grompp -c ../nvt/nvt.gro -t ../nvt/nvt.cpt -f $MDP_DIR/water_npt.mdp -p ../topol.top -o npt.tpr
#mdrun -deffnm npt
#cd ..


mkdir ProductionMD
cd ProductionMD
grompp -c ../npt/npt.gro -t ../npt/npt.cpt -f $MDP_DIR/water_ProductionMD.mdp -p ../topol.top -o md.tpr
mdrun -deffnm md
cd ..





