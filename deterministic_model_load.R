# calculate change in frequency of a single locus under selection or drift
# define the selection coefficient (0-1)
s<- -0.001
#starting allele frequency for the simulation
q<-0.25
#dominance (0-1)
h<-0.5
#set number of generations to simulate 
num_gen<-5000

#initialize frequency array
freq<-{}

#run selection for num_gen
for (x in 1:num_gen){
q = (q*(1-q)*(1+h*s) + q^2*(1+s))/((1-q)^2 + 2*q*(1-q)*(1+h*s) + q^2*(1+s))
freq<-c(freq,q)
} 

plot(1:5000,freq,type="l",lwd=2,col="darkblue",ylab="Minor parent ancestry",xlab="Generations since admixture",ylim=c(0,1),xlim=c(0,5000),cex.lab=1.4,cex.axis=1.2)
