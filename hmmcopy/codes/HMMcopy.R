library(HMMcopy)
library(grid)
library(ggplot2)
#library(lattice)
library(gridExtra)


args <- commandArgs(trailingOnly = TRUE)
rfile <- args[1]
gfile <- args[2]
mfile <- args[3]


# setwd("E:/liangcai")
# rfile <- "readcounts.seg.wig"
# gfile <- "genome.gc.wig"
# mfile <- "genome.fa.map.w5000.wig"
# Correcting  and   visualizing   tumour   copy   number  proﬁles

tumor_reads <- wigsToRangedData(rfile,gfile,mfile)
tumor_reads[1000:1010,]
#tumor_reads$chr <- lapply(tumor_reads$chr,factor)
tumour_copy<- correctReadcount(tumor_reads,samplesize = 5000)
tumour_copy[1000:1010,]
pdf("plotCorrection.pdf",width = 12,height =4)
par(mar=c(4,4,2,0))
plotCorrection(tumour_copy, pch=".")
dev.off()

for (i in c(1:22)){
  i1=as.character(i)
  file_name = paste0(paste0("chr",i1,sep=""),".pdf",sep="")
  print(file_name)
  tumour_copy1 = tumour_copy[chr == i1,]
  tumour_segments <- HMMsegment(tumour_copy1)
  pdf(file_name)
  par(mfrow =c(1,1))
  par(cex.main=0.5,cex.lab =0.5,cex.axis=0.5,mar =c(2,1.5,0,0),mgp=c(1,0.5, 0))
  plotSegments(tumour_copy1,tumour_segments,pch = ".",
               ylab= "Tumour Copy Number",xlab = "Chromosome Position")
  cols<-stateCols()
  legend("topleft", c("HOMD", "HETD", "NEUT", "GAIN","AMPL","HLAMP"),
         fill = cols, horiz = TRUE,bty = "n",cex = 0.5)
  dev.off()
}



plotSegments <- function(correctOutput, segmentOutput,
                         chr = correctOutput$chr[1], ...){
  if (is.null(segmentOutput$segs)) {
    warning("Processed segments now found, automatically processing")
    segmentOutput$segs <- processSegments(segments$segs,
                                          correctOutput$chr, correctOutput$start, correctOutput$end,
                                          correctOutput$copy)
  }
  
  segs <- segmentOutput$segs
  correctOutput$state <- segmentOutput$state
  cols <- stateCols()
  range <- quantile(correctOutput$copy, na.rm = TRUE, prob = c(0.01, 0.99))
  
  a <- subset(correctOutput, chr == chr)
  b <- subset(segs, chr == chr)
  #plot(a$start, a$copy,col = cols[as.numeric(as.character(a$state))], ylim = range, ...)
  if (chr==1){
    p1 <- ggplot(data = a, aes(start,copy))+
      geom_point(colour = cols[as.numeric(as.character(a$state))],size = 0.25)+ylim(-0.8,1)+theme_classic()+ 
      scale_x_discrete(breaks=NULL)+labs(x = paste0("chr",chr,sep=""), y = "Tumor Copy Number")+
      theme(plot.margin = unit(c(0,0,0,0),units = "lines"))
  }else{
    p1 <- ggplot(data = a, aes(start,copy))+
      geom_point(colour = cols[as.numeric(as.character(a$state))],
                                                       size = 0.25)+
      ylim(-0.8,1)+theme_classic()+ scale_x_discrete(breaks=NULL)+
      labs(x = paste0("chr",chr,sep=""), y = "")+theme(plot.margin = unit(c(0,0,0,0),units = "lines"),
                                                       axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank())
  }
  

  for (k in 1:nrow(b)){
    p1 <- p1 + geom_segment(x=b$start[k], y=b$median[k],xend=b$end[k],yend=b$median[k],colour = "green",linetype = 1,size =1)
  }

  #for (k in 1:nrow(b)){lines(c(b$start[k], b$end[k]), rep(b$median[k], 2), lwd = 3,col = "green")}
  return (p1)
}

plotchromosome <- function(i, tumour_copy = tumour_copy){
  i1=as.character(i)
  tumour_copy1 = tumour_copy[chr == i1,]
  tumour_segments <- HMMsegment(tumour_copy1)
  p1 = plotSegments(tumour_copy1,tumour_segments,pch = ".",
                    ylab= "",xlab = "Chromosome Position")
  return (p1)
}


p1 <- plotchromosome(i=1,tumour_copy = tumour_copy)
p2 <- plotchromosome(i=2,tumour_copy = tumour_copy)
p3 <- plotchromosome(i=3,tumour_copy = tumour_copy)
p4 <- plotchromosome(i=4,tumour_copy = tumour_copy)
p5 <- plotchromosome(i=5,tumour_copy = tumour_copy)
p6 <- plotchromosome(i=6,tumour_copy = tumour_copy)
p7 <- plotchromosome(i=7,tumour_copy = tumour_copy)
p8 <- plotchromosome(i=8,tumour_copy = tumour_copy)
p9 <- plotchromosome(i=9,tumour_copy = tumour_copy)
p10 <- plotchromosome(i=10,tumour_copy = tumour_copy)
p11 <- plotchromosome(i=11,tumour_copy = tumour_copy)
p12 <- plotchromosome(i=12,tumour_copy = tumour_copy)
p13 <- plotchromosome(i=13,tumour_copy = tumour_copy)
p14 <- plotchromosome(i=14,tumour_copy = tumour_copy)
p15 <- plotchromosome(i=15,tumour_copy = tumour_copy)
p16 <- plotchromosome(i=16,tumour_copy = tumour_copy)
p17 <- plotchromosome(i=17,tumour_copy = tumour_copy)
p18 <- plotchromosome(i=18,tumour_copy = tumour_copy)
p19 <- plotchromosome(i=19,tumour_copy = tumour_copy)
p20 <- plotchromosome(i=20,tumour_copy = tumour_copy)
p21 <- plotchromosome(i=21,tumour_copy = tumour_copy)
p22 <- plotchromosome(i=22,tumour_copy = tumour_copy)
px <- plotchromosome(i="X",tumour_copy = tumour_copy)


pdf("Hmmcopy.pdf",width=16,height=4)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,px,nrow=1)
dev.off()



# list1=list()
# for (i in c(1:4)){
#   i1=as.character(i)
#   tumour_copy1 = tumour_copy[chr == i1,]
#   tumour_segments <- HMMsegment(tumour_copy1)
#   p1 = plotSegments(tumour_copy1,tumour_segments,pch = ".",
#                      ylab= "",xlab = "Chromosome Position")
#   list1 <-  list(list1,p1)
# }
# grid.newpage()
# pushViewport(viewport(layout = grid.layout(1,length(list1))))
# for (i in (1:length(list1))){
#   print(list1[i],vp = viewport(layout.pos.row = 1,layout.pos.col = i))
# }
# 
# grid.arrange(list1,nrow=1)


# 
# grid.newpage()
# pushViewport(viewport(layout = grid.layout(3,3)))
# for (i in c(1:3)){
#   i1=as.character(i)
#   print(i1)
#   tumour_copy1 = tumour_copy[chr == i1,]
#   tumour_segments <- HMMsegment(tumour_copy1)
#   vp=viewport(plotSegments(tumour_copy1,tumour_segments,pch = ".",
#                            ylab= "Tumour Copy Number",xlab = "Chromosome Position"),layout.pos.row = 1,layout.pos.col = i)
# }

# tumour_copy1 = tumour_copy[chr == 1,]
# tumour_segments <- HMMsegment(tumour_copy1)
# p1 <- plotSegments(tumour_copy1,tumour_segments,pch = ".",
#                    ylab= "Tumour Copy Number",xlab = "Chromosome Position")
# p1
# 
# tumour_copy1 = tumour_copy[chr == 2,]
# tumour_segments <- HMMsegment(tumour_copy1)
# p2 <- plotSegments(tumour_copy1,tumour_segments,pch = ".",
#                    ylab= "Tumour Copy Number",xlab = "Chromosome Position")
# grid.newpage()
# pushViewport(viewport(layout = grid.layout(1,2)))
# print(p1,vp = viewport(layout.pos.row = 1,layout.pos.col = 1))
# print(p2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
# 
# multiplot(p1,p2,cols = 2)



