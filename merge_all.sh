#!/bin/sh
#$ -S /bin/sh
#$ -t 1-10
#$ -tc 10
#$ -l s_vmem=4G -l mem_req=4G
#$ -cwd
#$ -o $HOME/average_physical/error_merge
#$ -e $HOME/average_physical/error_merge

#chr_num=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
#p_name=(A_philicity B_DNA_twist DNA_bending_stiffness DNA_denaturation Duplex_disrupt_energy Duplex_free_energy P_DNA_twist Propeller_twist Protein_induced_deformability_Bp Stabilizing_energy_of_Z_DNA_AS Stabilizing_energy_of_Z_DNA_SA Stacking_energy Bruckner flexibility)
p_name=(DNA_bending_stiffness DNA_denaturation Duplex_disrupt_energy Duplex_free_energy Protein_induced_deformability_Bp Stabilizing_energy_of_Z_DNA_AS Stabilizing_energy_of_Z_DNA_SA Stacking_energy Bruckner flexibility)

b=$((SGE_TASK_ID%10))

i_p_name=${p_name[$b]}

cat Result/${i_p_name}/${i_p_name}_chr*.txt > all_data/${i_p_name}_chr_all.txt
