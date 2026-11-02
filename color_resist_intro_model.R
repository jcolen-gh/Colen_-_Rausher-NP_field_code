
#this first loop sets the matrix of all possible gametes
zn<-5 #number of size loci
nc<-2^(zn+1)
#holds each gamete
cmb<-matrix(data=NA, nrow=nc, ncol=zn+1)
ilist<-rep(NA,nc)
index<-1
for (i in 1:nc){ #generalized equations for generating combinations of zn size and 1 color locus gametes
  for (j in 1:(zn+1)){
    if (i%%(2^(j))>=2^(j-1)){
      cmb[i,j]<-1
    } else {
      cmb[i,j]<-0
    } 
    index<-index+cmb[i,j]*2^(j-1)
  }
  ilist[i]<-index #gives unique identifier to each gamete, corresponding with cmb
  index<-1
}

t1<-split(cmb,row(cmb)) #used for faster indexing later on




#need to make pairwise matrix, with row X column as each row
#then separately make gamete distributions, fitness assessments, and bases for recursion
gmtrx<-matrix(data=NA, nrow=nc^2, ncol=nc+3)

#Scenario 1:for seed set by weight -flower size
Ww=20.4 #fitness of white-limbed lines on seed set (intercept)
#B=1.8782 #represents raw output of flower size on seed set ###(0.032,0.034) is output
B=1.8782*10/9.4755 #normalized for range of 10; i.e. B*10/(x-range of s.c on flower size)
#Scenario 2: by weight - sugar concentration
Ww=15.6 #white fitness for sugar concentration
#B=1.177 ### (0.024,0.026) is output
B=1.177*10/15.52667#normalized for range of 10; i.e. B*10/(x-range of s.c on flower size)

#3) for seed set by est number -flower size
Ww=957.4 #b/c white lines are slightly negatively sloped, this is the exact value where the lines cross
#B=88.79 ### (0.32, 0.034) is output
B=88.79*10/9.4755 #normalized as above
#4) by est.number - s.c.
Ww=771.1 
#B=64.29 ### (0.028, 0.03)
B=64.29*10/15.52667
#represents normalized for differences between max step-sizes (i.e. largest range of z -5 to 5) and max differences in seed set

phold=NA #for starting hold
gc() #next loop creates the matrix structure for the fitness and likelihood of gamete formation for each genotype
for (i in 1:nrow(gmtrx)){
  gmtrx[i,1]<-ilist[ceiling(i/2^(zn+1))] #creates gamete crosses
  if (i%%(2^(zn+1))==0){
    gmtrx[i,2]<-ilist[2^(zn+1)]
  } else{
    gmtrx[i,2]<-ilist[i%%2^(zn+1)]
  }
  if (t1[[gmtrx[i,1]]][1]==0&&t1[[gmtrx[i,2]]][1]==0){ #line is white
    W=Ww
  }  else {
    zs1<-t1[[gmtrx[i,1]]][2:(zn+1)]
    zs2<-t1[[gmtrx[i,1]]][2:(zn+1)]
    zscore<-sum(zs1)+sum(zs2)-zn
    W=Ww+B*zscore
  }
  gmtrx[i,3]<-W
  #first find differences between characters
  dspots<-abs(t1[[gmtrx[i,1]]]-t1[[gmtrx[i,2]]])
  ldiff<-sum(dspots)
  imx<-matrix(data=NA, nrow=2^ldiff, ncol=length(dspots)) #next, produce possible gametes created (where to place distributions)
  counter=0
  for (j in 1:length(dspots)){
    nrowvec<-seq(1,2^ldiff)%%2^(counter+1)>=2^(counter) #i%%(2^(j))>=2^(j-1)
    if (dspots[j]==0) {
      imx[,j]=t1[[gmtrx[i,1]]][j]
    } else {
      imx[nrowvec==TRUE,j]=t1[[gmtrx[i,1]]][j]
      imx[nrowvec==FALSE,j]=t1[[gmtrx[i,2]]][j]
      counter=counter+1
    }
  }
  for (j in 1:2^ldiff){ #now to match gametes with right columns
    for (k in 1:nrow(cmb)){
      if (all(imx[j,]==t1[[k]])){ #try to match vector
        phold=k
      }
    }
    gmtrx[i,3+phold]<-(0.5)^ldiff
    phold=0
  }
  #next need to index numbers to frequencies of the haplotypes, and run the equation
}

################################ NON-ITERATED CODE
#need to set initial frequencies
pstart=0.6 #arbitrary, will try across all combinations, possibly through for loops
z1s=0.2
z2s=0.4
z3s=0.6
z4s=0.8
z5s=1.0
#z6s=0.10
#z7s=0.01
startvec=c(pstart,z1s,z2s,z3s,z4s,z5s) #modify depending on number of starting loci
ave(startvec)
hold=rep(NA,length(startvec))
#assuming HWE, frequencies of haplotypes should be evenly distributed given frequencies
freq<-rep(NA,nc)
for (i in 1:nrow(cmb)){
  for (j in 1:length(startvec)) {
    hold[j]<-cmb[i,j]*startvec[j]+(1-cmb[i,j])*(1-startvec[j])
  }
  print(hold)
  freq[i]<- prod(hold)
}

#For migration model, do this instead:
m<-0.1 #migration paramater
freq<-rep(0,nc)
freq[nc]<-m
whapl<-freq
freq[nc-1]<-(1-m)

#for migration determined on color value
m0<-0.01 #low migration favored in all conditions - when flowers are pink
m1<-0.19 #high migration when lines are white

m<-m0
freq<-rep(0,nc)
freq[nc]<-m
freq[nc-1]<-(1-m)


#no continual migration
gen=500 #number of generations simulated
freqmatrix<-matrix(data=NA, nrow=gen,ncol=length(freq))
dfreqmatrix<-matrix(data=NA, nrow=gen,ncol=length(freq))
wbarvec<-rep(NA, length(gen))
gc()
for (i in 1:gen){
    tvar<-rep(NA,nrow(gmtrx))
    nfreq<-rep(NA,ncol(gmtrx)-3)
    for (j in 1:(ncol(gmtrx)-3)){
      for (k in 1:nrow(gmtrx)){
        if (is.na(gmtrx[k,j+3])){} else{
          tvar[k]<-gmtrx[k,j+3]*(gmtrx[k,3]/mean(gmtrx[,3]))*freq[gmtrx[k,1]]*freq[gmtrx[k,2]]
          print(c(i,j,k))
        }
      }
      nfreq[j]=sum(tvar, na.rm=T)
      print(sum(tvar, na.rm=T))
      tvar<-rep(NA,nrow(gmtrx))
    }
    #bug here, nfreq does not sum to 1
    #for now, normalize haplotype frequencies so that sums to 1
    newfreq<-nfreq/sum(nfreq)
    wbarvec[i]<-sum(nfreq)
    freqmatrix[i,]<-freq
    freq<-newfreq
}

#with continual migration
#For migration model, do this instead:
m<-0.04 #migration paramater
freq<-rep(0,nc)
freq[nc]<-m
whapl<-freq
freq[nc-1]<-(1-m)
gen=500 #number of generations simulated
freqmatrix<-matrix(data=NA, nrow=gen,ncol=length(freq))
dfreqmatrix<-matrix(data=NA, nrow=gen,ncol=length(freq))
pfreq2<-wbarvec<-rep(NA, length(gen))
gc()
for (i in 1:gen){
  tvar<-rep(NA,nrow(gmtrx))
  nfreq<-rep(NA,ncol(gmtrx)-3)
  for (j in 1:(ncol(gmtrx)-3)){
    for (k in 1:nrow(gmtrx)){
      if (is.na(gmtrx[k,j+3])){
        #print(NA)
      } else{
        tvar[k]<-(gmtrx[k,j+3]*gmtrx[k,3]*freq[gmtrx[k,1]]*freq[gmtrx[k,2]])
        #print(c(i,j,k))
      }
    }
    nfreq[j]=sum(tvar, na.rm=T)
    tvar<-rep(NA,nrow(gmtrx))
  }
  print(i)
  #bug here, nfreq does not sum to 1
  #for now, normalize haplotype frequencies so that sums to 1
  newfreq<-nfreq/sum(nfreq)
  pfreq2[i]<-sum(newfreq[c(TRUE,FALSE)])
  wbarvec[i]<-sum(nfreq)
  dhap<-newfreq-freq
  freqmatrix[i,]<-freq
  dfreqmatrix[i,]<-dhap
  freq<-newfreq*(1-m)+whapl
}

#no LD

#all odds are freq(white)/32
#all evens are freq(pink)/32
m<-0.015 #migration paramater
freq<-whapl<-rep(0,nc)
for (i in 1:length(freq)){
  if(i%%2==0) {
    freq[i]<-whapl[i]<-2*m/length(freq)
  } else if(i%%2==1) {
    freq[i]<-2*(1-m)/length(freq)
  }
}

gen=500 #number of generations simulated
freqmatrix<-matrix(data=NA, nrow=gen,ncol=length(freq))
dfreqmatrix<-matrix(data=NA, nrow=gen,ncol=length(freq))
wbarvec<-rep(NA, length(gen))
gc()
#run program from above with different LDs, modifying migration

#cycled migration 
nm<-5 #number of cycles that migration happens in
gen=500 #number of generations simulated
freqmatrix<-matrix(data=NA, nrow=gen,ncol=length(freq))
dfreqmatrix<-matrix(data=NA, nrow=gen,ncol=length(freq))
wbarvec<-rep(NA, length(gen))
gc()
for (i in 1:gen){
  tvar<-rep(NA,nrow(gmtrx))
  nfreq<-rep(NA,ncol(gmtrx)-3)
  for (j in 1:(ncol(gmtrx)-3)){
    for (k in 1:nrow(gmtrx)){
      if (is.na(gmtrx[k,j+3])){} else{
        tvar[k]<-gmtrx[k,j+3]*gmtrx[k,3]*freq[gmtrx[k,1]]*freq[gmtrx[k,2]]
        #print(c(i,j,k))
      }
    }
    nfreq[j]=sum(tvar, na.rm=T)
    tvar<-rep(NA,nrow(gmtrx))
  }
  #bug here, nfreq does not sum to 1
  #for now, normalize haplotype frequencies so that sums to 1
  newfreq<-nfreq/sum(nfreq)
  wbarvec[i]<-sum(nfreq)
  dhap<-newfreq-freq
  freqmatrix[i,]<-freq
  dfreqmatrix[i,]<-dhap
  if (i%%nm==0){
    freq<-newfreq*(1-m)+whapl
  } else{
    freq<-newfreq
  }
  
}

#for dependent migration

#with continual migration
gen=500 #number of generations simulated
freqmatrix<-matrix(data=NA, nrow=gen,ncol=length(freq))
dfreqmatrix<-matrix(data=NA, nrow=gen,ncol=length(freq))
wbarvec<-migr<-rep(NA, length(gen))
gc()
for (i in 1:gen){
  tvar<-rep(NA,nrow(gmtrx))
  nfreq<-rep(NA,ncol(gmtrx)-3)
  for (j in 1:(ncol(gmtrx)-3)){
    for (k in 1:nrow(gmtrx)){
      if (is.na(gmtrx[k,j+3])){} else{
        tvar[k]<-gmtrx[k,j+3]*gmtrx[k,3]*freq[gmtrx[k,1]]*freq[gmtrx[k,2]]
        #print(c(i,j,k))
      }
    }
    nfreq[j]=sum(tvar, na.rm=T)
    tvar<-rep(NA,nrow(gmtrx))
  }
  #bug here, nfreq does not sum to 1
  #for now, normalize haplotype frequencies so that sums to 1
  newfreq<-nfreq/sum(nfreq)
  wbarvec[i]<-sum(nfreq)
  dhap<-newfreq-freq
  freqmatrix[i,]<-freq
  dfreqmatrix[i,]<-dhap
  pfq<-sum(newfreq[(c(TRUE,FALSE))])
  wfq<-sum(newfreq[(c(FALSE, TRUE))])
  m<-m0*pfq+m1*wfq
  migr[i]<-m
  whapf<-rep(0, nc)
  whapf[nc]<-m
  freq<-newfreq*(1-m)+whapf
}

nm<-10 #number of cycles that migration happens in
gen=500 #number of generations simulated
freqmatrix<-matrix(data=NA, nrow=gen,ncol=length(freq))
dfreqmatrix<-matrix(data=NA, nrow=gen,ncol=length(freq))
wbarvec<-migr<-rep(NA, length(gen))
gc()
for (i in 1:gen){
  tvar<-rep(NA,nrow(gmtrx))
  nfreq<-rep(NA,ncol(gmtrx)-3)
  for (j in 1:(ncol(gmtrx)-3)){
    for (k in 1:nrow(gmtrx)){
      if (is.na(gmtrx[k,j+3])){} else{
        tvar[k]<-gmtrx[k,j+3]*gmtrx[k,3]*freq[gmtrx[k,1]]*freq[gmtrx[k,2]]
        #print(c(i,j,k))
      }
    }
    nfreq[j]=sum(tvar, na.rm=T)
    tvar<-rep(NA,nrow(gmtrx))
  }
  #bug here, nfreq does not sum to 1
  #for now, normalize haplotype frequencies so that sums to 1
  newfreq<-nfreq/sum(nfreq)
  wbarvec[i]<-sum(nfreq)
  dhap<-newfreq-freq
  freqmatrix[i,]<-freq
  dfreqmatrix[i,]<-dhap
  pfq<-sum(newfreq[(c(TRUE,FALSE))])
  wfq<-sum(newfreq[(c(FALSE, TRUE))])
  m<-m0*pfq+m1*wfq
  migr[i]<-m
  whapf<-rep(0, nc)
  whapf[nc]<-m
  if (i%%nm==0){
    freq<-newfreq*(1-m)+whapf
  } else{
    freq<-newfreq
  }
}


#now based on size
m0<-0.01 #low migration when flowers are largest
m1<-0.21 #high migration when flowers are smallest

m<-m0
freq<-rep(0,nc)
freq[nc]<-m
whapl<-freq
freq[nc-1]<-(1-m)
nsize<-rowSums(cmb[,2:(zn+1)])
#when nsize=0,eff m/haplotype is m0, and when nsize=5, effective migration is m1
em<-m0+m1*(zn-nsize)/zn
wempt<-rep(0,nc)
wempt[nc]<-1


gen=500 #number of generations simulated
freqmatrix<-matrix(data=NA, nrow=gen,ncol=length(freq))
dfreqmatrix<-matrix(data=NA, nrow=gen,ncol=length(freq))
wbarvec<-migr<-rep(NA, length(gen))
gc()
for (i in 1:gen){
  tvar<-rep(NA,nrow(gmtrx))
  nfreq<-rep(NA,ncol(gmtrx)-3)
  for (j in 1:(ncol(gmtrx)-3)){
    for (k in 1:nrow(gmtrx)){
      if (is.na(gmtrx[k,j+3])){} else{
        tvar[k]<-gmtrx[k,j+3]*gmtrx[k,3]*freq[gmtrx[k,1]]*freq[gmtrx[k,2]]
        #print(c(i,j,k))
      }
    }
    nfreq[j]=sum(tvar, na.rm=T)
    tvar<-rep(NA,nrow(gmtrx))
  }
  #bug here, nfreq does not sum to 1
  #for now, normalize haplotype frequencies so that sums to 1
  newfreq<-nfreq/sum(nfreq)
  wbarvec[i]<-sum(nfreq)
  dhap<-newfreq-freq
  freqmatrix[i,]<-freq
  dfreqmatrix[i,]<-dhap
  migr[i]<-m
  m<-weighted.mean(em,newfreq)
  freq<-newfreq*(1-m)+m*wempt
}

#cycled
nm<-10 #number of cycles that migration happens in
gen=500 #number of generations simulated
freqmatrix<-matrix(data=NA, nrow=gen,ncol=length(freq))
dfreqmatrix<-matrix(data=NA, nrow=gen,ncol=length(freq))
wbarvec<-migr<-rep(NA, length(gen))
gc()
for (i in 1:gen){
  tvar<-rep(NA,nrow(gmtrx))
  nfreq<-rep(NA,ncol(gmtrx)-3)
  for (j in 1:(ncol(gmtrx)-3)){
    for (k in 1:nrow(gmtrx)){
      if (is.na(gmtrx[k,j+3])){} else{
        tvar[k]<-gmtrx[k,j+3]*gmtrx[k,3]*freq[gmtrx[k,1]]*freq[gmtrx[k,2]]
        #print(c(i,j,k))
      }
    }
    nfreq[j]=sum(tvar, na.rm=T)
    tvar<-rep(NA,nrow(gmtrx))
  }
  #bug here, nfreq does not sum to 1
  #for now, normalize haplotype frequencies so that sums to 1
  newfreq<-nfreq/sum(nfreq)
  wbarvec[i]<-sum(nfreq)
  dhap<-newfreq-freq
  freqmatrix[i,]<-freq
  dfreqmatrix[i,]<-dhap
  pfq<-sum(newfreq[(c(TRUE,FALSE))])
  wfq<-sum(newfreq[(c(FALSE, TRUE))])
  m<-m0*pfq+m1*wfq
  migr[i]<-m
  whapf<-rep(0, nc)
  whapf[nc]<-m
  if (i%%nm==0){
    freq<-newfreq*(1-m)+whapf
  } else{
    freq<-newfreq
  }
  
}

#first plot - haplotype freq of each allele state by time



plot(freqmatrix[,2^(zn+1)-1], pch=20, xlab="Number of Genrations", ylab="Frequency of Haplotype with Color and all Size Loci") #can modify for change

plot(x=-1, pch=20, xlab="Number of Genrations", ylab="Haplotype Frequency", xlim=c(0,500), ylim=c(0,1)) #can modify for change
for (i in 1:ncol(freqmatrix)){
  if (t1[[i]][1]==1){
    color<-2
  } else {
    color<-1
  }
  shape<-20
  points(freqmatrix[,i], col=color, pch=shape, cex=0.6)
}

#next color allele frequency +size frequency
pfreq<-freqmatrix[,c(TRUE,FALSE)]
wfreq<-freqmatrix[,c(FALSE, TRUE)]
pft<-rowSums(pfreq)
wft<-rowSums(wfreq)

cmigr<-rep(NA,500)
for (i in 1:length(cmigr)){
  cmigr[i]<-(1-m)^i
}
cmigr<-seq()
plot(pft, pch=20, xlab="Number of Generations", ylab="Frequency of Color Allele", col="magenta", ylim=c(0,1), type="l", lwd=3)
points(pfreq2, pch=20, xlab="Number of Generations", ylab="Frequency of Color Allele", col="red4", ylim=c(0,1), type="l", lwd=3)
points(cmigr, pch=20, col="grey", lwd=0.5)

lines(wft,type="l", pch=20, lwd=3)
h<-recordPlot()
h
#points(pft, pch=20, col="violet")
#points(wft, pch=20, col="lightgrey")

plot(migr, pch=20)

dpfreq<-dfreqmatrix[,c(TRUE,FALSE)]
dwfreq<-dfreqmatrix[,c(FALSE, TRUE)]
dpft<-rowSums(dpfreq)
dwft<-rowSums(dwfreq)

plot(dpft, pch=20, xlab="Number of Generations", ylab="Frequency of Color Allele", col="magenta", ylim=c(min(dwft),max(dpft)))
points(dwft, pch=20)

#size loci - done manually
olds1<-olds2<-olds3<-olds4<-olds5<-olds6<-olds7<-NA
s1s<-s2s<-s3s<-s4s<-s5s<-s6s<-s7s<-NA
for (i in 1:nrow(cmb)){
  if (cmb[i,2]==1){
    s1s<-i
    s1v<-c(olds1,s1s)
    olds1<-s1v
  }  else{}
  if (cmb[i,3]==1){
    s2s<-i
    s2v<-c(olds2,s2s)
    olds2<-s2v
  }  else{}
  if (cmb[i,4]==1){
    s3s<-i
    s3v<-c(olds3,s3s)
    olds3<-s3v
  }  else{}
  if (cmb[i,5]==1){
    s4s<-i
    s4v<-c(olds4,s4s)
    olds4<-s4v
  }  else{}
  if (cmb[i,6]==1){
    s5s<-i
    s5v<-c(olds5,s5s)
    olds5<-s5v
  }  else{}
  if (0==1){
    s6s<-i
    s6v<-c(olds6,s6s)
    olds6<-s6v
  }  else{}
  if (0==1){
    s7s<-i
    s7v<-c(olds7,s7s)
    olds7<-s7v
  }  else{}
}
s1<-s1v[2:(nrow(cmb)/2+1)]
s2<-s2v[2:(nrow(cmb)/2+1)]
s3<-s3v[2:(nrow(cmb)/2+1)]
s4<-s4v[2:(nrow(cmb)/2+1)]
s5<-s5v[2:(nrow(cmb)/2+1)]
#s6<-s6v[2:(nrow(cmb)/2+1)]
#s7<-s7v[2:(nrow(cmb)/2+1)]

s1freq<-freqmatrix[,s1]
s1ft<-rowSums(s1freq)
#g
#points(s1ft, pch=20,col="red")
#points(s1ft, pch=20,col="coral")
plot(s1ft, pch=20, xlab="Number of Generations", ylab="Frequency of Size Allele", col="black", ylim=c(0,1))
s2freq<-freqmatrix[,s2]
s2ft<-rowSums(s2freq)
points(s2ft, pch=20, col="darkorange")
points(s2ft, pch=20, col="orange")
s3freq<-freqmatrix[,s3]
s3ft<-rowSums(s3freq)
points(s3ft, pch=20, col="darkgreen")
points(s3ft, pch=20, col="green")
s4freq<-freqmatrix[,s4]
s4ft<-rowSums(s4freq)
points(s4ft, pch=20, col="blue")
points(s4ft, pch=20, col="lightblue")
s5freq<-freqmatrix[,s5]
s5ft<-rowSums(s5freq)
points(s5ft, pch=20, col="purple")
points(s5ft, pch=20, col="violet")
points(pft, pch=20, col="magenta")

saviorg<-c(pft[500],s5ft[500])

g<-recordPlot()
s6freq<-freqmatrix[,s6]
s6ft<-rowSums(s6freq)
points(s6ft, pch=20, col="goldenrod")
s7freq<-freqmatrix[,s7]
s7ft<-rowSums(s7freq)
points(s7ft, pch=20, col="burlywood")
g<-recordPlot()

f
plot(s1ft, pft, pch=20, main="for 1 initial condition, change of size allele freq with color", xlab="Frequency of Major Size Allele" , ylab="Pink Frequency", type="l", lwd=2)
points(s1ft, pft, pch=20, col="red", type="l", lwd=2)
points(s1ft, pft, pch=20, col="coral")
points(s1ft[1],pft[1], pch=1, cex=2)
points(s2ft, pft, pch=20, col="darkorange")
points(s2ft, pft, pch=20, col="orange")
points(s2ft[1],pft[1], pch=16)
points(s3ft, pft, pch=20, col="darkgreen")
points(s3ft, pft, pch=20, col="green")
points(s3ft[1],pft[1], pch=16)
points(s4ft, pft, pch=20, col="blue")
points(s4ft, pft, pch=20, col="lightblue")
points(s4ft[1],pft[1], pch=16)
points(s5ft, pft, pch=20, col="purple")
points(s5ft, pft, pch=20, col="violet")
points(s5ft[1],pft[1], pch=16)
points(s6ft, pft, pch=20, col="goldenrod")
points(s7ft, pft, pch=20, col="burlywood")

f<-recordPlot()


#contour map?
library(plot3D)
colfunc <- colorRampPalette(c(l, l+1))
ourcolpal<-cbind(sort(wbarvec),colfunc(length(wbarvec)))

snav<-(s1ft+s2ft+s3ft+s4ft+s5ft)/5
plot(snav, pft, pch=20, main="Change of average of size allele freq with color", xlab="Average Frequency of Size Allele" , ylab="Pink Frequency", xlim=c(0,1), ylim=c(0,max(pft)), col=ourcolpal[match(wbarvec, ourcolpal),2])
points(snav[1],pft[1], pch=16,cex=1)
e<-recordPlot()

###############################FOR ITERATED, CONTINUAL MIGRATION

#with continual migration
gen=500 #number of generations simulated
colfuncl <- colorRampPalette(c("black","lightgrey"))
fixcolol<-colfuncl(20)
plot(x=-1, pch=20, main="Trajectory of size, color allele", xlab="Frequency of Major Size Allele" , ylab="Pink Frequency", xlim=c(0,1), ylim=c(0,1), type="l", lwd=2)
ff<-recordPlot()
grads<-20
holdends<-holdstart<-data.frame(NA, nrow=grads, ncol=5)
pvdata<-matrix(NA, nrow=grads*gen,ncol=5)
oldtdf<-NA
for (z in 1:grads){
  m<-0.002*(grads+1)-0.002*z
  freq<-rep(0,nc)
  freq[nc]<-m
  whapl<-freq
  freq[nc-1]<-(1-m)
  freqmatrix<-matrix(data=NA, nrow=gen,ncol=length(freq))
  wbarvec<-rep(NA, length(gen))
  print(z)
  for (i in 1:gen){
    tvar<-rep(NA,nrow(gmtrx))
    nfreq<-rep(NA,ncol(gmtrx)-3)
    for (j in 1:(ncol(gmtrx)-3)){
      for (k in 1:nrow(gmtrx)){
        if (is.na(gmtrx[k,j+3])){} else{
          tvar[k]<-gmtrx[k,j+3]*gmtrx[k,3]*freq[gmtrx[k,1]]*freq[gmtrx[k,2]]
          #print(c(i,j,k))
        }
      }
      nfreq[j]=sum(tvar, na.rm=T)
      tvar<-rep(NA,nrow(gmtrx))
    }
    #normalize haplotype frequencies so that sums to 1
    newfreq<-nfreq/sum(nfreq)
    wbarvec[i]<-sum(nfreq)
    freqmatrix[i,]<-newfreq
    freq<-newfreq*(1-m)+whapl
  }
  pfreq<-freqmatrix[,c(TRUE,FALSE)]
  pft<-rowSums(pfreq)
  #h
  #points(pft, pch=20, xlab="Number of Generations", ylab="Frequency of Color Allele", col=coldet, ylim=c(0,1), type="l", lwd=3)
  #h<-recordPlot()
  #print(h)
  print(pft[500])
  olds1<-NA
  s1s<-NA
  for (i in 1:nrow(cmb)){
    if (cmb[i,2]==1){
      s1s<-i
      s1v<-c(olds1,s1s)
      olds1<-s1v
    }  else{}
  }
  s1<-s1v[2:(nrow(cmb)/2+1)]
  s1freq<-freqmatrix[,s1]
  s1ft<-rowSums(s1freq)
  #g
  #points(s1ft, pch=20, col=coldet, type="l")
  #g<-recordPlot()
  print(s1ft[500])
  #print(g)
  ff
  colfuncl3 <- colorRampPalette(c("#C3FD8B","#A4863A"))
  fixcolol3<-colfuncl3(100)
  colfuncl4 <- colorRampPalette(c("#A4863A","deeppink3"))
  fixcolol4<-colfuncl4(400)
  fixcolol2<-c(fixcolol3,fixcolol4) #to make gradient more striking
  tdf<-as.data.frame(cbind(s1ft,pft,fixcolol2, seq(1,gen), rep(z,gen)))
  points(as.numeric(tdf$pft)~as.numeric(tdf$s1ft), pch=20, col=tdf$fixcolol2)
  holdends[z,]<-tdf[500,]
  holdstart[z,]<-tdf[1,]
  holdstart[z,3]<-fixcolol[z]
  #points(s1ft[1],pft[1], pch=20,col=fixcolol[z], cex=1.5)
  #points(as.numeric(tdf$pft[75])~as.numeric(tdf$s1ft[75]), pch=20,col=tdf$fixcolol2[75], cex=2.5) #test point for vis
  pvdata<-rbind(oldtdf,tdf)
  oldtdf<-pvdata
  ff<-recordPlot()
}
colnames(holdends)<-colnames(holdstart)<-colnames(tdf)[1:3]
points(as.numeric(holdends$pft)~as.numeric(holdends$s1ft), pch=20, col=holdends$fixcolol2, cex=1.5)
points(as.numeric(holdstart$pft)~as.numeric(holdstart$s1ft), pch=20, col=holdstart$fixcolol2, cex=1.5)

plot(as.numeric(pvdata$pft)~as.numeric(pvdata$s1ft), pch=20, col=pvdata$fixcolol2, main="Trajectory of size, color allele", xlab="Frequency of Major Size Allele" , ylab="Frequency of Color Allele")
plot(as.numeric(pvdata$pft)~as.numeric(pvdata$V4), pch=20, col=pvdata$fixcolol2, main="Trajectory of size, color allele", ylab="Frequency of Color Allele" , xlab="Generation",cex=0.8)
plot(as.numeric(pvdata$s1ft)~as.numeric(pvdata$V4), pch=20, col=pvdata$fixcolol2, main="Trajectory of size, color allele", ylab="Frequency of Major Size Allele" , xlab="Generation",cex=0.8)
#need to add starting and ending points to 3 graphs by doing z+500 trick

