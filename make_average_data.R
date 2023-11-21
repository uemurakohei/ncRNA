PPname <- as.character(commandArgs(trailingOnly=TRUE)[1]);#
chrxx <- as.character(commandArgs(trailingOnly=TRUE)[2]);#chrxx<-"chrX"
path.FANTOM_CAGE <- as.character(commandArgs(trailingOnly=TRUE)[3]);#
path.Physical_properties <- as.character(commandArgs(trailingOnly=TRUE)[4]);#
path.out <- as.character(commandArgs(trailingOnly=TRUE)[5]);#

print(c(PPname, chrxx, path.FANTOM_CAGE, path.Physical_properties, path.out))

dir.create(path.out)
setwd(path.out);
if(!sum(list.files(path.out) == PPname)){
  dir.create(PPname);
}
setwd(PPname);

distance_fromTSS.U <- 150
distance_fromTSS.D <- 150

antisense_adjust_num <- 1
if(length(grep("Bruckner", path.Physical_properties))){
  antisense_adjust_num <- 2
}else if(length(grep("flexibility", path.Physical_properties))){
  antisense_adjust_num <- 3
}

FANTOM_CAGE <- read.table(path.FANTOM_CAGE , stringsAsFactors=F);
Physical_properties <- scan(path.Physical_properties, what=1)

FANTOM_CAGE_ex <- FANTOM_CAGE[,c(1,4,5,6,7,8)]
colnames(FANTOM_CAGE_ex) <- c("chr", "gene_id", "TIE_score", "strand", "start", "end");

dat_TSS_sense <- FANTOM_CAGE_ex[which(FANTOM_CAGE_ex[,"strand"]=="+"),"start"] +1;
dat_TSS_anti <- FANTOM_CAGE_ex[which(FANTOM_CAGE_ex[,"strand"]=="-"),"end"];

dat_gene_id_sense <- FANTOM_CAGE_ex[which(FANTOM_CAGE_ex[,"strand"]=="+"),"gene_id"]
dat_gene_id_anti <- FANTOM_CAGE_ex[which(FANTOM_CAGE_ex[,"strand"]=="-"),"gene_id"]

if(!(length(dat_TSS_sense) == length(dat_gene_id_sense))){
  print("the number of TSS and id incorrect")
  q("no")
}

if(!(length(dat_TSS_anti) == length(dat_gene_id_anti))){
  print("the number of TSS and id incorrect")
  q("no")
}


PP_all <- c(); gene_id.all <- c();
for(i_tss_sense in 1:length(dat_TSS_sense)){
  i_tss_Pos <- as.numeric(dat_TSS_sense[i_tss_sense])
  if(((i_tss_Pos-distance_fromTSS.U)<0)|((i_tss_Pos+distance_fromTSS.D)>length(Physical_properties))){next}
  i_PP <- Physical_properties[(i_tss_Pos-distance_fromTSS.U):(i_tss_Pos+distance_fromTSS.D)];
  if(sum(i_PP == 55555)){next}
  gene_id.all <- c(gene_id.all, dat_gene_id_sense[i_tss_sense])
  PP_all <- rbind(PP_all, i_PP)
}
for(i_tss_anti in 1:length(dat_TSS_anti)){
  i_tss_Pos <- as.numeric(dat_TSS_anti[i_tss_anti])
  if(((i_tss_Pos-distance_fromTSS.D)<0)|((i_tss_Pos+distance_fromTSS.U)>length(Physical_properties))){next}
  i_PP <- Physical_properties[(i_tss_Pos+distance_fromTSS.U-antisense_adjust_num):(i_tss_Pos-distance_fromTSS.D-antisense_adjust_num)];
  if(sum(i_PP == 55555)){next}
  gene_id.all <- c(gene_id.all, dat_gene_id_anti[i_tss_anti])
  PP_all <- rbind(PP_all, i_PP)
}

rownames(PP_all) <- gene_id.all

pdf(sprintf("%s_%s.pdf",PPname,chrxx))
plot(apply(PP_all, 2, mean), type="l")
dev.off()

#PP_all_norm <- (PP_all-min(PP_all))*(1-(-1))/(max(PP_all)-min(PP_all))+(-1);
#plot(apply(PP_all_norm, 2, mean), typel="l")

write.table(PP_all, col.names=F, row.names=T, file=sprintf("%s_%s.txt",PPname,chrxx))
