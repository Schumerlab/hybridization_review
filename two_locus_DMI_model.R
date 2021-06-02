
#####DMI

#Simulates change in allele frequency as a result of two locus selection based on Karlin 1975, Lewontin & Kojima

# % parent 1 ancestry
M=0.25
B=1-M

#set number of generations to simulate
numreps=2000

# initial haplotype frequencies
x<-c((1-B)*M, (1-B)*(1-M), B*M, B*(1-M))
ancestry_track1<-{}
ancestry_track2<-{}

# Set fitness and dominance
s1<-0.1
s2<-0
hbc=0
ha=0.5

#fitness function for the diagonal is multiplicative here for coevolving incompatibilities but can be set differently
diagonal=(1-((1-ha)*hbc*s1))*(1-(ha*(1-hbc)*s2))

# fitness table
# w11 and w44 are conspecific homozygotes; w22 and w33 are heterospecific homozygotes, and so on

w<-function(s){matrix(c(
    1,       1-(1-hbc)*s2,         1-(1-ha)*s1,     diagonal,
    1-(1-hbc)*s2,       1-s2,      diagonal,     1-ha*s2,
    1-(1-ha)*s1,  diagonal,    1-s1,     1-hbc*s1,
    diagonal,  1-ha*s2,    1-hbc*s1,     1
    ), nrow=4,ncol=4)}
    
    
# intitialise and store haplotype frequency over generations
X<-x
d <- x[1]*x[4] - x[2]*x[3]
w_mean<-0
for (j in 1:4){
	for (k in 1:4){
		w_mean = w_mean + as.numeric(x[j]*x[k]*w(s)[j,k])
	}
}
w_bar<-w_mean


# start simulation here
for (z in 1:numreps){ # number of generations
    # disequilibrium
	D = x[1]*x[4] - x[2]*x[3]; D
	d<-c(d,as.numeric(D)) # tracks D over time
    
    #from allele frequencies assuming random mating
    f_M1=x[1]+x[2]
    f_M2=x[1]+x[3]
    
    #store allele frequencies
    ancestry_track1<-c(ancestry_track1,f_M1)
    ancestry_track2<-c(ancestry_track2,f_M2)
    
    # marginal fitnesses
	wstar<-{}
	for (i in 1:4){
        wstar<-c(wstar, x[1]*w(s)[i,1] + x[2]*w(s)[i,2] + x[3]*w(s)[i,3] + x[4]*w(s)[i,4])
	}
    # mean fitness
	w_mean<-0
	for (j in 1:4){
		for (k in 1:4){
			w_mean = w_mean + as.numeric(x[j]*x[k]*w(s)[j,k])
		}
	}
	w_bar<-c(w_bar, w_mean)
    # new haplotype frequencies
	y<-{}
	y<-c(y,as.numeric(round(x[1]*(wstar[1]/w_mean)-(0.5*D*w(s)[1,4]/w_mean),10)))
	y<-c(y,as.numeric(round(x[2]*(wstar[2]/w_mean)+(0.5*D*w(s)[1,4]/w_mean),10)))
	y<-c(y,as.numeric(round(x[3]*(wstar[3]/w_mean)+(0.5*D*w(s)[1,4]/w_mean),10)))
	y<-c(y,as.numeric(round(x[4]*(wstar[4]/w_mean)-(0.5*D*w(s)[1,4]/w_mean),10)))
	x<-y
	X<-rbind(X,as.numeric(y)) # store information
} # for numreps


plot(ancestry_track1,ylim=c(0,0.6),col="darkblue",pch=20,cex.lab=1.4,cex.axis=1.2,ylab="Minor parent ancestry",xlab="Generations since admixture",type="l",lwd=3,xlim=c(0,2000))

points(ancestry_track2,col="darkblue",pch=20,cex.lab=1.3,cex.axis=1.2,ylab="Proportion Parent1",xlab="Generations",type="l",lwd=3,lty=2)


M=0.4
B=1-M

#set number of generations to simulate
numreps=500

# initial haplotype frequencies
x<-c((1-B)*M, (1-B)*(1-M), B*M, B*(1-M))
ancestry_track1<-{}
ancestry_track2<-{}

# Set fitness and dominance
s1<-0.1
s2<-0.1
hbc=0.5
ha=0.5

#fitness function for the diagonal is multiplicative here for coevolving incompatibilities but can be set differently
diagonal=(1-((1-ha)*hbc*s1))*(1-(ha*(1-hbc)*s2))

# fitness table
# w11 and w44 are conspecific homozygotes; w22 and w33 are heterospecific homozygotes, and so on

w<-function(s){matrix(c(
  1,       1-(1-hbc)*s2,         1-(1-ha)*s1,     diagonal,
  1-(1-hbc)*s2,       1-s2,      diagonal,     1-ha*s2,
  1-(1-ha)*s1,  diagonal,    1-s1,     1-hbc*s1,
  diagonal,  1-ha*s2,    1-hbc*s1,     1
), nrow=4,ncol=4)}


# intitialise and store haplotype frequency over generations
X<-x
d <- x[1]*x[4] - x[2]*x[3]
w_mean<-0
for (j in 1:4){
  for (k in 1:4){
    w_mean = w_mean + as.numeric(x[j]*x[k]*w(s)[j,k])
  }
}
w_bar<-w_mean


# start simulation here
for (z in 1:numreps){ # number of generations
  # disequilibrium
  D = x[1]*x[4] - x[2]*x[3]; D
  d<-c(d,as.numeric(D)) # tracks D over time
  
  #from allele frequencies assuming random mating
  f_M1=x[1]+x[2]
  f_M2=x[1]+x[3]
  
  #store allele frequencies
  ancestry_track1<-c(ancestry_track1,f_M1)
  ancestry_track2<-c(ancestry_track2,f_M2)
  
  # marginal fitnesses
  wstar<-{}
  for (i in 1:4){
    wstar<-c(wstar, x[1]*w(s)[i,1] + x[2]*w(s)[i,2] + x[3]*w(s)[i,3] + x[4]*w(s)[i,4])
  }
  # mean fitness
  w_mean<-0
  for (j in 1:4){
    for (k in 1:4){
      w_mean = w_mean + as.numeric(x[j]*x[k]*w(s)[j,k])
    }
  }
  w_bar<-c(w_bar, w_mean)
  # new haplotype frequencies
  y<-{}
  y<-c(y,as.numeric(round(x[1]*(wstar[1]/w_mean)-(0.5*D*w(s)[1,4]/w_mean),10)))
  y<-c(y,as.numeric(round(x[2]*(wstar[2]/w_mean)+(0.5*D*w(s)[1,4]/w_mean),10)))
  y<-c(y,as.numeric(round(x[3]*(wstar[3]/w_mean)+(0.5*D*w(s)[1,4]/w_mean),10)))
  y<-c(y,as.numeric(round(x[4]*(wstar[4]/w_mean)-(0.5*D*w(s)[1,4]/w_mean),10)))
  x<-y
  X<-rbind(X,as.numeric(y)) # store information
} # for numreps

plot(ancestry_track1,ylim=c(0,0.6),col="darkblue",pch=20,cex.lab=1.4,cex.axis=1.2,ylab="Minor parent ancestry",xlab="Generations since admixture",type="l",lwd=3,xlim=c(0,500))

points(ancestry_track2,col="darkblue",pch=20,cex.lab=1.3,cex.axis=1.2,ylab="Proportion Parent1",xlab="Generations",type="l",lwd=3,lty=2)


