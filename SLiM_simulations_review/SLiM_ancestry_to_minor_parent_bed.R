arrArgs<-commandArgs(trailingOnly = TRUE);

options(scipen=999)

infile<-as.character(arrArgs[1])
data<-read.table(file=infile,sep="\t",head=TRUE,as.is=T)

chrlength<-31064592

warning("chromosome length hard coded to match SLiM, update if needed\n")

data<-subset(data,data$end<=chrlength)

trim_data<-data[-c(1,2,3)]

#head(trim_data)

numcol<-length(trim_data[1,])
keep<-{}
for(y in 2:length(trim_data[,1])){

focal<-trim_data[y,]
#write.table(focal,col.names=FALSE)

if(mean(t(focal)/2) < 0.1){
keep<-rbind(keep,data[y,])
}

}

write.table(keep,file=paste(infile,"_fixed_minor.bed",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
