arrArgs<-commandArgs(trailingOnly = TRUE);

infile<-as.character(arrArgs[1])
data<-read.csv(file=infile,sep="\t",head=FALSE)

binfile<-as.character(arrArgs[2])
bins<-read.csv(file=binfile,sep="\t",head=FALSE)
options(scipen=999)

windows<-{}

for (x in 1:length(bins[,1])){

focal<-subset(data,data$V1>=bins[,2][x] & data$V1<=bins[,3][x])

p5_p4_af<-mean(abs(focal$V6-focal$V5))
p4_p3_af<-mean(abs(focal$V5-focal$V4))

#focal$V2 - pop1
#focal$V3 - pop2
#focal$V4 - pop3
#focal$V5 - pop4
#focal$V6 - pop5

#p5 (V6) versus p4 (V5) is admixed pop versus sister
#p3 (V4) versus p4 (V5) is unadmixed pop versus sister

p5_p4_fixed<-length(subset(focal$V6,focal$V5>0.9 & focal$V6<0.1 | focal$V5<0.1 & focal$V6>0.9))
p4_p3_fixed<-length(subset(focal$V5,focal$V4>0.9 & focal$V5<0.1 | focal$V4<0.1 & focal$V5>0.9))
p5_p4_fixed_derived<-length(subset(focal$V6, focal$V5<0.1 & focal$V6>0.9 ))
p4_p3_fixed_derived<-length(subset(focal$V5, focal$V5>0.9 & focal$V4<0.1 ))

if(length(focal[,1])>0){
windows<-rbind(windows,cbind(bins[x,],p5_p4_af,p4_p3_af,p5_p4_fixed,p5_p4_fixed_derived,p4_p3_fixed,p4_p3_fixed_derived))
}#if data

}#go through all bins

write.table(windows,file=paste(infile,"_windows.txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
