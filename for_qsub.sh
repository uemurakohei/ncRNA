#!/bin/sh
#$ -S /bin/sh
#$ -t 1-240
#$ -tc 240
#$ -l s_vmem=32G -l mem_req=32G
#$ -cwd
#$ -o $HOME/average_physical/error_output
#$ -e $HOME/average_physical/error_output

chr_num=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
#p_name=(A_philicity B_DNA_twist DNA_bending_stiffness DNA_denaturation Duplex_disrupt_energy Duplex_free_energy P_DNA_twist Propeller_twist Protein_induced_deformability_Bp Stabilizing_energy_of_Z_DNA_AS Stabilizing_energy_of_Z_DNA_SA Stacking_energy Bruckner flexibility)
p_name=(DNA_bending_stiffness DNA_denaturation Duplex_disrupt_energy Duplex_free_energy Protein_induced_deformability_Bp Stabilizing_energy_of_Z_DNA_AS Stabilizing_energy_of_Z_DNA_SA Stacking_energy Bruckner flexibility)

a=$((SGE_TASK_ID%24))
b=$((SGE_TASK_ID/24%10))

i_chr=${chr_num[$a]}
i_p_name=${p_name[$b]}

echo $i_chr
echo $i_p_name

$HOME/local/bin/R --vanilla --args $i_p_name "chr${i_chr}" ${HOME}/FANTOMCATlv3robust/lncRNA_all/chr${i_chr}.bed ${HOME}/promoter2/${i_p_name}/chr${i_chr}.txt Result <make_average_data.R> error.log
