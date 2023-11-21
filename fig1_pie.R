library(openxlsx)

data<-read.xlsx("/Users/uemurakohei/Desktop/Supplementary data (paper2)/Suppl. 1 lncRNA element.xlsx", sheet=1)

cols <- c(rainbow(10, alpha=0.8), rainbow((nrow(data)-1), alpha=0.8)[-c(1:10)])
names(cols)<-data[,1][1:(length(data[,1])-1)]

for(i_plot in 2:ncol(data)){

  pie_e<-data[1:(nrow(data)-1),i_plot]
  names(pie_e)<-names(cols)<-data[,1][1:(length(data[,1])-1)]

  sort_pie_e <- sort(pie_e, decreasing=T)
  sort_cols <- cols[names(sort_pie_e)]

  pdf(sprintf("%s.pdf", colnames(data)[i_plot])
  , height = 1, width = 1
  )
  par(mar = c(0, 0, 0, 0))
  par(oma = c(0, 0, 0, 0))
  par(family = "Helvetica",cex = 0.8)

  pie(sort_pie_e,radius=1,
    labels="",
    #main=sprintf("%s-antisense", c_DHS),
    col=sort_cols,clockwise=T,border="#ffffff",cex.main=2)
  dev.off()

}
