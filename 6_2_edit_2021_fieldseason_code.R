getwd()
setwd("C:/Users/jonat/OneDrive/Documents/Grad School Research Files")
getwd()
masterdata <- read.csv("C:/Users/jonat/OneDrive/Documents/Grad School Research Files/fielddatamaster3.csv")

#reminder, make sure to look at observed color as well as line color!

install.packages("lubridate")
install.packages("lme4")
install.packages("lmerTest")
install.packages("car")
install.packages("bbmle")
install.packages("plot3D")
install.packages("MBA")

library(lme4)
library(lmerTest)
library(car)

#checking which field observations don't match up with which line colors
n <-nrow(masterdata)
checksheet <-dvec <- rep(NA,n)
for(i in 1:n) {
  if (as.character(masterdata[i,3]) %in% as.character(masterdata[i,10])) {
    next
  } else {
    nrow <- masterdata[i,]
    if (is.na(nrow[10])){
      next
    } else {
      checksheet <-rbind(checksheet,nrow)
    }
  }
}
checklist <- checksheet[,1:3]
checklist2 <- checksheet[,10]
checklist <-cbind(checklist, checklist2)
checklist[,2] # 33 total samples, 17 are either w/L or L/W misIDs, so 16 represent P/wor W/pflips
# 10 of those are from line 610
#other 6 are A4H, E4G, E4H, F3L, G1B, M2A
#E4G and E4H are likely flipped
pwchecker <- rep(NA,nrow (checklist))
for(i in 1:nrow(checklist)) {
  if (isTRUE(checklist[i,3] == "L"))  {
    next
  } else if (isTRUE(checklist[i,4] == "L")) {
    next
  } else {
    nrow <- checklist[i,]
    pwchecker <-rbind(pwchecker,nrow)
  }
}
nrow(pwchecker) #17 samples represent incorret color plantings, only 2 are w -> p; 10 are line 610


#checking RIL line 610
line610 <- subset(masterdata, RILID==610)

plot(masterdata$Actual..~masterdata$Projected.., xlab="Projected Seed Number",ylab="Actual Seed Number",pch=floor(1.72-as.numeric(masterdata$Actual)/as.numeric(masterdata$Projected..))%%2*15+1)
abline(lm(masterdata$Actual..~masterdata$Projected..))
abline(0,1,col="red")
plot(log(masterdata$Actual..)~log(masterdata$Projected..), xlab="Log Projected Seed Number",ylab="Log Actual Seed Number",pch=floor(1.72-as.numeric(masterdata$Actual)/as.numeric(masterdata$Projected..))%%2*15+1)
abline(lm(log(masterdata$Actual..)~log(masterdata$Projected..)))
abline(0,1,col="red")

#calculating differnces in early growth based off of flower color
summary(lm(masterdata$Early.Growth..cm.~masterdata$RILcolor))
boxplot(masterdata$Early.Growth..cm.~masterdata$RILcolor)
nolac <- subset(masterdata, RILcolor!="L")
summary(lm(nolac$Early.Growth..cm.~nolac$RILcolor))
Anova(lmer(Early.Growth..cm.~RILcolor+(1|RILID), data=nolac))
summary(lm(nolac$Shoot..~nolac$RILcolor))

summary(lm(nolac$Herbiv~nolac$RILcolor))
#following up on herbiv with a nested ANOVA wtith treatment as block
#defining block as variable

#from ft code
blockID <- gsub('\\b(\\pL)|.','\\U\\1',nolac$X,perl = TRUE)#adapted from stackoverflow
nolac<-cbind(nolac,blockID)
nolac<-nolac[1:798,]
nolac$Total.seed.weight<-as.numeric(nolac$Total.seed.weight)+0.000001
nolac$Projected..<-round(as.numeric(nolac$Projected..))



#instead, let's think of a fitness function for lac growth
library(plot3D)
lacdata =subset(masterdata, RILcolor=="L")
xcord<-ycord<-rep(NA, nrow(lacdata))

lettodist <- function(x){
  if(nchar(x)==3){
    latnumb=NA
    longnumb=NA
    chars <- unlist(strsplit(x, ""))
    if (chars[1]=="A"||chars[1]=="E"||chars[1]=="I"||chars[1]=="M"){
      latnumb=0
    } else if (chars[1]=="B"||chars[1]=="F"||chars[1]=="J"||chars[1]=="N"){
      latnumb=48
    }else if (chars[1]=="C"||chars[1]=="G"||chars[1]=="K"||chars[1]=="O"){
      latnumb=96
    }else if (chars[1]=="D"||chars[1]=="H"||chars[1]=="L"||chars[1]=="P"){
      latnumb=144
    }
    if (chars[1]=="A"||chars[1]=="B"||chars[1]=="C"||chars[1]=="D"){
      longnumb=0
    } else if (chars[1]=="E"||chars[1]=="F"||chars[1]=="G"||chars[1]=="H"){
      longnumb=25
    } else if (chars[1]=="I"||chars[1]=="J"||chars[1]=="K"||chars[1]=="L"){
      longnumb=50
    } else if (chars[1]=="M"||chars[1]=="N"||chars[1]=="O"||chars[1]=="P"){
      longnumb=75
    } else {
      latnumb=longnumb=NA
    }
    if (is.na(as.numeric(chars[2])=="FALSE")){
      latnumb=longnumb=NA
    } else {
      longnumb=longnumb+5*(as.numeric(chars[2])-1)
    }
    
    if (chars[3]=="B"){
      latnumb=latnumb+3
    } else if (chars[3]=="C"){
      latnumb=latnumb+6
    } else if (chars[3]=="D"){
      latnumb=latnumb+9
    } else if (chars[3]=="E"){
      latnumb=latnumb+12
    } else if (chars[3]=="F"){
      latnumb=latnumb+15
    } else if (chars[3]=="G"){
      latnumb=latnumb+18
    } else if (chars[3]=="H"){
      latnumb=latnumb+21
    } else if (chars[3]=="I"){
      latnumb=latnumb+24
    } else if (chars[3]=="J"){
      latnumb=latnumb+27
    } else if (chars[3]=="K"){
      latnumb=latnumb+30
    } else if (chars[3]=="L"){
      latnumb=latnumb+33
    } else if (chars[3]=="M"){
      latnumb=latnumb+36
    } else if (chars[3]=="N"){
      latnumb=latnumb+39
    } else if (chars[3]=="O"){
      latnumb=latnumb+42
    }
  } else{
    latnumb =longnumb<-"NA"
  }
  return(c(latnumb, longnumb))
}


for (i in 1:nrow(lacdata)){
  xci <- lettodist(lacdata$X[i])[1]
  yci <- lettodist(lacdata$X[i])[2]
  xcord[i]<-as.numeric(xci)
  ycord[i]<-as.numeric(yci)
}
lacdata<-cbind(lacdata,xcord, ycord)


growthmap=lm(lacdata$Early.Growth..cm.~lacdata$xcord*lacdata$ycord+I(lacdata$xcord^2)+I(lacdata$ycord^2))
growthmap3=lm(lacdata$Early.Growth..cm.~lacdata$xcord*lacdata$ycord+I(lacdata$xcord^2)+I(lacdata$ycord^2)+I(lacdata$xcord^3)+I(lacdata$ycord^3)+lacdata$ycord*I(lacdata$xcord^2)+lacdata$xcord*I(lacdata$ycord^2))
growthmap4=lm(lacdata$Early.Growth..cm.~lacdata$xcord*lacdata$ycord+I(lacdata$xcord^2)+I(lacdata$ycord^2)
              +I(lacdata$xcord^3)+I(lacdata$ycord^3)+lacdata$ycord*I(lacdata$xcord^2)+lacdata$xcord*I(lacdata$ycord^2)
              +I(lacdata$xcord^4)+I(lacdata$ycord^4)+lacdata$ycord*I(lacdata$xcord^3)+lacdata$xcord*I(lacdata$ycord^3)
              +I(lacdata$ycord^2)*I(lacdata$xcord^2))
gm0<-lm(lacdata$Early.Growth..cm.~1)
gstep<-step(growthmap4, direction=c("both"), scope=formula(gm0))
growthmapbf<- lm(lacdata$Early.Growth..cm.~lacdata$xcord + lacdata$ycord + 
                   I(lacdata$xcord^2) + I(lacdata$ycord^2) + I(lacdata$ycord^3) + 
                   I(lacdata$xcord^4) + I(lacdata$ycord^4) + lacdata$xcord:I(lacdata$ycord^2) + 
                   I(lacdata$xcord^2):I(lacdata$ycord^2))#from output in step func


growthre= function(x,y)
{coef(growthmap)[1] + coef(growthmap)[2]*x + coef(growthmap)[3]*y + coef(growthmap)[4]*x^2 + coef(growthmap)[5]*y^2 + coef(growthmap)[6]*x*y}
growthre3= function(x,y)
{coef(growthmap3)[1] + coef(growthmap3)[2]*x + coef(growthmap3)[3]*y + coef(growthmap3)[4]*x^2 + coef(growthmap3)[5]*y^2 + coef(growthmap3)[6]*x^3+coef(growthmap3)[7]*y^3+coef(growthmap3)[8]*x*y+coef(growthmap3)[9]*y*x^2+coef(growthmap3)[10]*x*y^2}
growthre4=function(x,y)
{coef(growthmap4)[1] + coef(growthmap4)[2]*x + coef(growthmap4)[3]*y + coef(growthmap4)[4]*x^2 + coef(growthmap4)[5]*y^2 + coef(growthmap4)[6]*x^3+coef(growthmap4)[7]*y^3+coef(growthmap4)[8]*x^4+coef(growthmap4)[9]*y^4+coef(growthmap4)[10]*x*y+coef(growthmap4)[11]*y*x^2+coef(growthmap4)[12]*x*y^2+coef(growthmap4)[13]*y*x^3+coef(growthmap4)[14]*x*y^3+coef(growthmap4)[15]*(x^2)*y^2}
growthrebf=function(x,y)
{coef(growthmapbf)[1] + coef(growthmapbf)[2]*x + coef(growthmapbf)[3]*y + coef(growthmapbf)[4]*x^2 + coef(growthmapbf)[5]*y^2 + coef(growthmapbf)[6]*y^3+coef(growthmapbf)[7]*x^4+coef(growthmapbf)[8]*y^4+coef(growthmapbf)[9]*x*y^2+coef(growthmapbf)[10]*(x^2)*y^2}

xaxis=seq(0,186, 5)
yaxis=seq(0,95, 3)

z_axis = outer(X = xaxis,Y = yaxis, growthre)
z_axis3 = outer(X = xaxis,Y = yaxis, growthre3)
z_axis4 = outer(X = xaxis,Y = yaxis, growthre4)
z_axisbf= outer(X = xaxis,Y = yaxis, growthrebf)


contour2D(x = xaxis,y = yaxis, z = z_axis, xlab = "x-coord", ylab = "y-coord", clab = "phyt func")
colfunc <- colorRampPalette(c("cyan", "red"))
ourcolpal<-cbind(sort(lacdata$Early.Growth..cm),colfunc(length(lacdata$Early.Growth..cm)))
points(x=lacdata$xcord, y=lacdata$ycord, col=ourcolpal[match(lacdata$Early.Growth..cm, ourcolpal),2], cex=0.5, pch=19)

contour2D(x = xaxis,y = yaxis, z = z_axis3, xlab = "x-coord", ylab = "y-coord", clab = "phyt func")
contour2D(x = xaxis,y = yaxis, z = z_axis4, xlab = "x-coord", ylab = "y-coord", clab = "phyt func")

#bf
contour2D(x = xaxis,y = yaxis, z = z_axisbf, xlab = "x-coord", ylab = "y-coord", clab = "phyt func")



persp3D(x = xaxis, y = yaxis, z = z_axis, phi = 45, theta = 55, xlab = "x-coord", ylab = "y-coord", zlab = "early growth (cm)", clab = "Phytometer", contour = TRUE, cex.axis = 0.75, cex.lab = 0.75,  zlim=(c(0,30))) #phi and theta change angle
scatter3D(x=lacdata$xcord, y=lacdata$ycord, z=lacdata$Early.Growth..cm., add = TRUE, type = "h", pch = 19, lwd = 0.2, colkey = FALSE)

persp3D(x = xaxis, y = yaxis, z = z_axis3, phi = 45, theta = 55, xlab = "x-coord", ylab = "y-coord", zlab = "early growth (cm)", clab = "Phytometer", contour = TRUE, cex.axis = 0.75, cex.lab = 0.75,  zlim=(c(0,30))) #phi and theta change angle
scatter3D(x=lacdata$xcord, y=lacdata$ycord, z=lacdata$Early.Growth..cm., add = TRUE, type = "h", pch = 19, lwd = 0.2, colkey = FALSE)

persp3D(x = xaxis, y = yaxis, z = z_axis4, phi = 45, theta = 55, xlab = "x-coord", ylab = "y-coord", zlab = "early growth (cm)", clab = "Phytometer", contour = TRUE, cex.axis = 0.75, cex.lab = 0.75,  zlim=(c(0,30))) #phi and theta change angle
scatter3D(x=lacdata$xcord, y=lacdata$ycord, z=lacdata$Early.Growth..cm., add = TRUE, type = "h", pch = 19, lwd = 0.2, colkey = FALSE)

persp3D(x = xaxis, y = yaxis, z = z_axisbf, phi = 45, theta = 55, xlab = "x-coord", ylab = "y-coord", zlab = "early growth (cm)", clab = "Phytometer", contour = TRUE, cex.axis = 0.75, cex.lab = 0.75,  zlim=(c(0,30))) #phi and theta change angle
scatter3D(x=lacdata$xcord, y=lacdata$ycord, z=lacdata$Early.Growth..cm., add = TRUE, type = "h", pch = 19, lwd = 0.2, colkey = FALSE)

scatter3D(x=lacdata$xcord, y=lacdata$ycord, z=lacdata$Early.Growth..cm., cex=1, pch=20, phi = 45, theta = 55, type = "h", xlab = "x-coord", ylab = "y-coord", zlab = "early growth (cm)", clab = "Phytometer")
#need to compare models
gm0<-lm(lacdata$Early.Growth..cm.~1)
gstep<-step(growthmap4, direction=c("both"), scope=formula(gm0))
growthmapbf<- lm(lacdata$Early.Growth..cm.~lacdata$xcord + lacdata$ycord + 
                   I(lacdata$xcord^2) + I(lacdata$ycord^2) + I(lacdata$ycord^3) + 
                   I(lacdata$xcord^4) + I(lacdata$ycord^4) + lacdata$xcord:I(lacdata$ycord^2) + 
                   I(lacdata$xcord^2):I(lacdata$ycord^2))

#now let's reduce lacdata to only sampled plots:
lacdatasm<-subset(lacdata, xcord<91)
lacdatasm<-subset(lacdatasm, ycord<46)
growthmapsm=lm(lacdatasm$Early.Growth..cm.~lacdatasm$xcord*lacdatasm$ycord+I(lacdatasm$xcord^2)+I(lacdatasm$ycord^2))
growthresm= function(x,y)
{coef(growthmapsm)[1] + coef(growthmapsm)[2]*x + coef(growthmapsm)[3]*y + coef(growthmapsm)[4]*x^2 + coef(growthmapsm)[5]*y^2 + coef(growthmapsm)[6]*x*y}
xaxis=seq(0,91, 1)
yaxis=seq(0,46, 1)

z_axis = outer(X = xaxis,Y = yaxis, growthresm)
contour2D(x = xaxis,y = yaxis, z = z_axis, xlab = "x-coord", ylab = "y-coord", clab = "phyt func")
colfunc <- colorRampPalette(c("cyan", "red"))
ourcolpal<-cbind(sort(lacdatasm$Early.Growth..cm),colfunc(length(lacdatasm$Early.Growth..cm)))
points(x=lacdatasm$xcord, y=lacdatasm$ycord, col=ourcolpal[match(lacdatasm$Early.Growth..cm, ourcolpal),2], cex=0.5, pch=19)
persp3D(x = xaxis, y = yaxis, z = z_axis, phi = 45, theta = 55, xlab = "x-coord", ylab = "y-coord", zlab = "early growth (cm)", clab = "Phytometer", contour = TRUE, cex.axis = 0.75, cex.lab = 0.75,  zlim=(c(0,30))) #phi and theta change angle
scatter3D(x=lacdatasm$xcord, y=lacdatasm$ycord, z=lacdatasm$Early.Growth..cm., add = TRUE, type = "h", pch = 19, lwd = 0.2, colkey = FALSE)

#now with more complex models
growthmapsm3=lm(lacdatasm$Early.Growth..cm.~lacdatasm$xcord*lacdatasm$ycord+I(lacdatasm$xcord^2)+I(lacdatasm$ycord^2)+I(lacdatasm$xcord^3)+I(lacdatasm$ycord^3)+lacdatasm$ycord*I(lacdatasm$xcord^2)+lacdatasm$xcord*I(lacdatasm$ycord^2))
growthmapsm4=lm(lacdatasm$Early.Growth..cm.~lacdatasm$xcord*lacdatasm$ycord+I(lacdatasm$xcord^2)+I(lacdatasm$ycord^2)
                +I(lacdatasm$xcord^3)+I(lacdatasm$ycord^3)+lacdatasm$ycord*I(lacdatasm$xcord^2)+lacdatasm$xcord*I(lacdatasm$ycord^2)
                +I(lacdatasm$xcord^4)+I(lacdatasm$ycord^4)+lacdatasm$ycord*I(lacdatasm$xcord^3)+lacdatasm$xcord*I(lacdatasm$ycord^3)
                +I(lacdatasm$ycord^2)*I(lacdatasm$xcord^2))
gm0<-lm(lacdatasm$Early.Growth..cm.~1)
gstep<-step(growthmapsm4, direction=c("both"), scope=formula(gm0))
growthmapsmbf<- lm(lacdatasm$Early.Growth..cm.~lacdatasm$xcord + lacdatasm$ycord + 
                     I(lacdatasm$xcord^2) + I(lacdatasm$ycord^2) + I(lacdatasm$xcord^3) +I(lacdatasm$ycord^3) + 
                     I(lacdatasm$ycord^4) +lacdatasm$ycord*lacdatasm$xcord+ lacdatasm$ycord*I(lacdatasm$xcord^2) + 
                     lacdatasm$ycord*I(lacdatasm$xcord^3))#from output in step func

growthresm3= function(x,y)
{coef(growthmapsm3)[1] + coef(growthmapsm3)[2]*x + coef(growthmapsm3)[3]*y + coef(growthmapsm3)[4]*x^2 + coef(growthmapsm3)[5]*y^2 + coef(growthmapsm3)[6]*x^3+coef(growthmapsm3)[7]*y^3+coef(growthmapsm3)[8]*x*y+coef(growthmapsm3)[9]*y*x^2+coef(growthmapsm3)[10]*x*y^2}
growthresm4=function(x,y)
{coef(growthmapsm4)[1] + coef(growthmapsm4)[2]*x + coef(growthmapsm4)[3]*y + coef(growthmapsm4)[4]*x^2 + coef(growthmapsm4)[5]*y^2 + coef(growthmapsm4)[6]*x^3+coef(growthmapsm4)[7]*y^3+coef(growthmapsm4)[8]*x^4+coef(growthmapsm4)[9]*y^4+coef(growthmapsm4)[10]*x*y+coef(growthmapsm4)[11]*y*x^2+coef(growthmapsm4)[12]*x*y^2+coef(growthmapsm4)[13]*y*x^3+coef(growthmapsm4)[14]*x*y^3+coef(growthmapsm4)[15]*(x^2)*y^2}
growthrebfsm=function(x,y)
{coef(growthmapsmbf)[1] + coef(growthmapsmbf)[2]*x + coef(growthmapsmbf)[3]*y + coef(growthmapsmbf)[4]*x^2 + coef(growthmapsmbf)[5]*y^2+ coef(growthmapsmbf)[6]*x^3 + coef(growthmapsmbf)[7]*y^3+
    +coef(growthmapsmbf)[8]*y^4+coef(growthmapsmbf)[9]*x*y+coef(growthmapsmbf)[10]*(x^2)*y+coef(growthmapsmbf)[11]*(x^3)*y}


z_axis3 = outer(X = xaxis,Y = yaxis, growthresm3)
z_axis4 = outer(X = xaxis,Y = yaxis, growthresm4)
z_axisbf = outer(X = xaxis,Y = yaxis, growthrebfsm)
contour2D(x = xaxis,y = yaxis, z = z_axis3, xlab = "x-coord", ylab = "y-coord", clab = "phyt func")
contour2D(x = xaxis,y = yaxis, z = z_axis4, xlab = "x-coord", ylab = "y-coord", clab = "phyt func")
contour2D(x = xaxis,y = yaxis, z = z_axisbf, xlab = "x-coord", ylab = "y-coord", clab = "phyt func")


persp3D(x = xaxis, y = yaxis, z = z_axis3, phi = 45, theta = 55, xlab = "x-coord", ylab = "y-coord", zlab = "early growth (cm)", clab = "Phytometer", contour = TRUE, cex.axis = 0.75, cex.lab = 0.75,  zlim=(c(0,30))) #phi and theta change angle
scatter3D(x=lacdatasm$xcord, y=lacdatasm$ycord, z=lacdatasm$Early.Growth..cm., add = TRUE, type = "h", pch = 19, lwd = 0.2, colkey = FALSE)

persp3D(x = xaxis, y = yaxis, z = z_axis4, phi = 45, theta = 55, xlab = "x-coord", ylab = "y-coord", zlab = "early growth (cm)", clab = "Phytometer", contour = TRUE, cex.axis = 0.75, cex.lab = 0.75,  zlim=(c(0,30))) #phi and theta change angle0
scatter3D(x=lacdatasm$xcord, y=lacdatasm$ycord, z=lacdatasm$Early.Growth..cm., add = TRUE, type = "h", pch = 19, lwd = 0.2, colkey = FALSE)

persp3D(x = xaxis, y = yaxis, z = z_axisbf, phi = 45, theta = 55, xlab = "x-coord", ylab = "y-coord", zlab = "early growth (cm)", clab = "Phytometer", contour = TRUE, cex.axis = 0.75, cex.lab = 0.75,  zlim=(c(0,30))) #phi and theta change angle
scatter3D(x=lacdatasm$xcord, y=lacdatasm$ycord, z=lacdatasm$Early.Growth..cm., add = TRUE, type = "h", pch = 19, lwd = 0.2, colkey = FALSE)


#out of curiosity, estimate fitness function for all points, including lines:
for (i in 1:nrow(masterdata)){
  xci <- lettodist(masterdata$X[i])[1]
  yci <- lettodist(masterdata$X[i])[2]
  xcord[i]<-as.numeric(xci)
  ycord[i]<-as.numeric(yci)
}
masterdata<-cbind(masterdata,xcord, ycord)


gmapall=lm(masterdata$Early.Growth..cm.~masterdata$xcord*masterdata$ycord+I(masterdata$xcord^2)+I(masterdata$ycord^2))
gmap3=lm(masterdata$Early.Growth..cm.~masterdata$xcord*masterdata$ycord+I(masterdata$xcord^2)+I(masterdata$ycord^2)+I(masterdata$xcord^3)+I(masterdata$ycord^3)+masterdata$ycord*I(masterdata$xcord^2)+masterdata$xcord*I(masterdata$ycord^2))
gmap4=lm(masterdata$Early.Growth..cm.~masterdata$xcord*masterdata$ycord+I(masterdata$xcord^2)+I(masterdata$ycord^2)
         +I(masterdata$xcord^3)+I(masterdata$ycord^3)+masterdata$ycord*I(masterdata$xcord^2)+masterdata$xcord*I(masterdata$ycord^2)
         +I(masterdata$xcord^4)+I(masterdata$ycord^4)+masterdata$ycord*I(masterdata$xcord^3)+masterdata$xcord*I(masterdata$ycord^3)
         +I(masterdata$ycord^2)*I(masterdata$xcord^2))
gmap0<-lm(masterdata$Early.Growth..cm.~1)
gstepall<-step(gmap4, direction=c("both"), scope=formula(gmap0))
gmapbf<- lm(masterdata$Early.Growth..cm.~ masterdata$ycord + 
              I(masterdata$xcord^3) + I(masterdata$ycord^2) + I(masterdata$ycord^3) + 
              I(masterdata$ycord^4))#from output in step func


galre= function(x,y)
{coef(gmapall)[1] + coef(gmapall)[2]*x + coef(gmapall)[3]*y + coef(gmapall)[4]*x^2 + coef(gmapall)[5]*y^2 + coef(gmapall)[6]*x*y}
galre3= function(x,y)
{coef(gmap3)[1] + coef(gmap3)[2]*x + coef(gmap3)[3]*y + coef(gmap3)[4]*x^2 + coef(gmap3)[5]*y^2 + coef(gmap3)[6]*x^3+coef(gmap3)[7]*y^3+coef(gmap3)[8]*x*y+coef(gmap3)[9]*y*x^2+coef(gmap3)[10]*x*y^2}
galre4=function(x,y)
{coef(gmap4)[1] + coef(gmap4)[2]*x + coef(gmap4)[3]*y + coef(gmap4)[4]*x^2 + coef(gmap4)[5]*y^2 + coef(gmap4)[6]*x^3+coef(gmap4)[7]*y^3+coef(gmap4)[8]*x^4+coef(gmap4)[9]*y^4+coef(gmap4)[10]*x*y+coef(gmap4)[11]*y*x^2+coef(gmap4)[12]*x*y^2+coef(gmap4)[13]*y*x^3+coef(gmap4)[14]*x*y^3+coef(gmap4)[15]*(x^2)*y^2}
galbf=function(x,y)
{coef(gmapbf)[1] + coef(gmapbf)[2]*y + coef(gmapbf)[3]*x^3 + coef(gmapbf)[4]*y^2 + coef(gmapbf)[5]*y^3+coef(gmapbf)[6]*y^4}

xaxis=seq(0,186, 1)
yaxis=seq(0,95, 1)

z_axis = outer(X = xaxis,Y = yaxis, galre)
z_axis3 = outer(X = xaxis,Y = yaxis, galre3)
z_axis4 = outer(X = xaxis,Y = yaxis, galre4)
z_axisbf= outer(X = xaxis,Y = yaxis, galbf)


contour2D(x = xaxis,y = yaxis, z = z_axis, xlab = "x-coord", ylab = "y-coord", clab = "phyt func")
colfunc <- colorRampPalette(c("cyan", "red"))
ourcolpal2<-cbind(sort(masterdata$Early.Growth..cm),colfunc(length(masterdata$Early.Growth..cm)))
points(x=masterdata$xcord, y=masterdata$ycord, col=ourcolpal2[match(masterdata$Early.Growth..cm, ourcolpal2),2], cex=0.5, pch=19)

contour2D(x = xaxis,y = yaxis, z = z_axis3, xlab = "x-coord", ylab = "y-coord", clab = "phyt func")
persp3D(x = xaxis, y = yaxis, z = z_axis, phi = 45, theta = 55, xlab = "x-coord", ylab = "y-coord", zlab = "early growth (cm)", clab = "Phytometer", contour = TRUE, cex.axis = 0.75, cex.lab = 0.75,  zlim=(c(0,30))) #phi and theta change angle
scatter3D(x=lacdata$xcord, y=lacdata$ycord, z=lacdata$Early.Growth..cm., add = TRUE, type = "h", pch = 19, lwd = 0.2, colkey = FALSE)

persp3D(x = xaxis, y = yaxis, z = z_axis3, phi = 45, theta = 55, xlab = "x-coord", ylab = "y-coord", zlab = "early growth (cm)", clab = "Phytometer", contour = TRUE, cex.axis = 0.75, cex.lab = 0.75,  zlim=(c(0,30))) #phi and theta change angle
scatter3D(x=lacdata$xcord, y=lacdata$ycord, z=lacdata$Early.Growth..cm., add = TRUE, type = "h", pch = 19, lwd = 0.2, colkey = FALSE)

persp3D(x = xaxis, y = yaxis, z = z_axis4, phi = 45, theta = 55, xlab = "x-coord", ylab = "y-coord", zlab = "early growth (cm)", clab = "Phytometer", contour = TRUE, cex.axis = 0.75, cex.lab = 0.75,  zlim=(c(0,30))) #phi and theta change angle
scatter3D(x=lacdata$xcord, y=lacdata$ycord, z=lacdata$Early.Growth..cm., add = TRUE, type = "h", pch = 19, lwd = 0.2, colkey = FALSE)

persp3D(x = xaxis, y = yaxis, z = z_axisbf, phi = 45, theta = 55, xlab = "x-coord", ylab = "y-coord", zlab = "early growth (cm)", clab = "Phytometer", contour = TRUE, cex.axis = 0.75, cex.lab = 0.75,  zlim=(c(0,30))) #phi and theta change angle
scatter3D(x=lacdata$xcord, y=lacdata$ycord, z=lacdata$Early.Growth..cm., add = TRUE, type = "h", pch = 19, lwd = 0.2, colkey = FALSE)

scatter3D(x=masterdata$xcord, y=masterdata$ycord, z=masterdata$Early.Growth..cm., cex=1, pch=20, phi = 45, theta = 55, type = "h", xlab = "x-coord", ylab = "y-coord", zlab = "early growth (cm)", clab = "Phytometer")


#now for smaller part:
datasm<-subset(masterdata, xcord<91)
datasm<-subset(datasm, ycord<46)

gmapsm=lm(datasm$Early.Growth..cm.~datasm$xcord*datasm$ycord+I(datasm$xcord^2)+I(datasm$ycord^2))
gresm= function(x,y)
{coef(gmapsm)[1] + coef(gmapsm)[2]*x + coef(gmapsm)[3]*y + coef(gmapsm)[4]*x^2 + coef(gmapsm)[5]*y^2 + coef(gmapsm)[6]*x*y}
xaxis=seq(0,91, 1)
yaxis=seq(0,46, 1)

z_axis = outer(X = xaxis,Y = yaxis, gresm)
contour2D(x = xaxis,y = yaxis, z = z_axis, xlab = "x-coord", ylab = "y-coord", clab = "phyt func")
colfunc <- colorRampPalette(c("cyan", "red"))
ourcolpal<-cbind(sort(lacdatasm$Early.Growth..cm),colfunc(length(lacdatasm$Early.Growth..cm)))
points(x=lacdatasm$xcord, y=lacdatasm$ycord, col=ourcolpal[match(lacdatasm$Early.Growth..cm, ourcolpal),2], cex=0.5, pch=19)
colfunc2 <- colorRampPalette(c("white", "black"))
ourcolpal2<-cbind(sort(datasm$Early.Growth..cm),colfunc2(length(datasm$Early.Growth..cm)))
points(x=datasm$xcord, y=datasm$ycord, col=ourcolpal2[match(datasm$Early.Growth..cm, ourcolpal2),2], cex=0.5, pch=19)

persp3D(x = xaxis, y = yaxis, z = z_axis, phi = 45, theta = 55, xlab = "x-coord", ylab = "y-coord", zlab = "early growth (cm)", clab = "Phytometer", contour = TRUE, cex.axis = 0.75, cex.lab = 0.75,  zlim=(c(0,30))) #phi and theta change angle
scatter3D(x=datasm$xcord, y=datasm$ycord, z=datasm$Early.Growth..cm., add = TRUE, type = "h", pch = 19, lwd = 0.2, colkey = FALSE)

#now with more complex models
gmapsm3=lm(datasm$Early.Growth..cm.~datasm$xcord*datasm$ycord+I(datasm$xcord^2)+I(datasm$ycord^2)+I(datasm$xcord^3)+I(datasm$ycord^3)+datasm$ycord*I(datasm$xcord^2)+datasm$xcord*I(datasm$ycord^2))
gmapsm4=lm(datasm$Early.Growth..cm.~datasm$xcord*datasm$ycord+I(datasm$xcord^2)+I(datasm$ycord^2)
           +I(datasm$xcord^3)+I(datasm$ycord^3)+datasm$ycord*I(datasm$xcord^2)+datasm$xcord*I(datasm$ycord^2)
           +I(datasm$xcord^4)+I(datasm$ycord^4)+datasm$ycord*I(datasm$xcord^3)
           +I(datasm$ycord^2)*I(datasm$xcord^3)+I(datasm$ycord^2)*I(datasm$xcord^4)+I(datasm$ycord^3)*I(datasm$xcord^4)+datasm$xcord*I(datasm$ycord^3)+I(datasm$ycord^3)*I(datasm$xcord^2)+I(datasm$ycord^4)*I(datasm$xcord^2)+I(datasm$ycord^4)*I(datasm$xcord^3)
           +I(datasm$ycord^2)*I(datasm$xcord^2)+I(datasm$xcord^3)*I(datasm$ycord^3)+I(datasm$xcord^4)*I(datasm$ycord^4))
gm0<-lm(datasm$Early.Growth..cm.~1)
gstep<-step(gmapsm4, direction=c("both"), scope=formula(gm0))
gmapsmbf<- lm(datasm$Early.Growth..cm.~datasm$xcord + datasm$ycord + 
                I(datasm$xcord^2) + I(datasm$ycord^2) +I(datasm$ycord^3) + 
                I(datasm$ycord^4) + datasm$ycord*I(datasm$xcord^2) + 
                datasm$xcord*I(datasm$ycord^3))#from output in step func

gresm3= function(x,y)
{coef(gmapsm3)[1] + coef(gmapsm3)[2]*x + coef(gmapsm3)[3]*y + coef(gmapsm3)[4]*x^2 + coef(gmapsm3)[5]*y^2 + coef(gmapsm3)[6]*x^3+coef(gmapsm3)[7]*y^3+coef(gmapsm3)[8]*x*y+coef(gmapsm3)[9]*y*x^2+coef(gmapsm3)[10]*x*y^2}
gresm4=function(x,y)
{coef(gmapsm4)[1] + coef(gmapsm4)[2]*x + coef(gmapsm4)[3]*y + coef(gmapsm4)[4]*x^2 + coef(gmapsm4)[5]*y^2 + coef(gmapsm4)[6]*x^3+coef(gmapsm4)[7]*y^3+coef(gmapsm4)[8]*x^4+coef(gmapsm4)[9]*y^4+coef(gmapsm4)[10]*x*y+coef(gmapsm4)[11]*y*x^2+coef(gmapsm4)[12]*x*y^2+coef(gmapsm4)[13]*y*x^3+coef(gmapsm4)[14]*(y^2)*x^3+coef(gmapsm4)[15]*(y^2)*x^4+coef(gmapsm4)[16]*(y^3)*x^4+coef(gmapsm4)[17]*x*y^3+coef(gmapsm4)[18]*(x^2)*y^3+coef(gmapsm4)[19]*(x^2)*y^4+coef(gmapsm4)[20]*(x^3)*y^4+coef(gmapsm4)[21]*(x^2)*y^2+coef(gmapsm4)[22]*(x^3)*y^3+coef(gmapsm4)[23]*(x^4)*y^4}
grebfsm=function(x,y)
{coef(gmapsmbf)[1] + coef(gmapsmbf)[2]*x + coef(gmapsmbf)[3]*y + coef(gmapsmbf)[4]*x^2 + coef(gmapsmbf)[5]*y^2+ coef(gmapsmbf)[6]*y^3 + coef(gmapsmbf)[7]*y^4+
    +coef(gmapsmbf)[8]*y*x^2+coef(gmapsmbf)[9]*x*y^3}


z_axis3 = outer(X = xaxis,Y = yaxis, gresm3)
z_axis4 = outer(X = xaxis,Y = yaxis, gresm4)
z_axisbf = outer(X = xaxis,Y = yaxis, grebfsm)
contour2D(x = xaxis,y = yaxis, z = z_axis3, xlab = "x-coord", ylab = "y-coord", clab = "phyt func")
contour2D(x = xaxis,y = yaxis, z = z_axis4, xlab = "x-coord", ylab = "y-coord", clab = "phyt func")
contour2D(x = xaxis,y = yaxis, z = z_axisbf, xlab = "x-coord", ylab = "y-coord", clab = "phyt func")

points(x=datasm$xcord, y=datasm$ycord, col=ourcolpal2[match(datasm$Early.Growth..cm, ourcolpal2),2], cex=0.5, pch=19)
points(x=lacdatasm$xcord, y=lacdatasm$ycord, col=ourcolpal[match(lacdatasm$Early.Growth..cm, ourcolpal),2], cex=0.5, pch=19)

persp3D(x = xaxis, y = yaxis, z = z_axis3, phi = 45, theta = 55, xlab = "x-coord", ylab = "y-coord", zlab = "early growth (cm)", clab = "Phytometer", contour = TRUE, cex.axis = 0.75, cex.lab = 0.75,  zlim=(c(0,30))) #phi and theta change angle
scatter3D(x=lacdatasm$xcord, y=lacdatasm$ycord, z=lacdatasm$Early.Growth..cm., add = TRUE, type = "h", pch = 19, lwd = 0.2, colkey = FALSE)

persp3D(x = xaxis, y = yaxis, z = z_axis4, phi = 45, theta = 55, xlab = "x-coord", ylab = "y-coord", zlab = "early growth (cm)", clab = "Phytometer", contour = TRUE, cex.axis = 0.75, cex.lab = 0.75,  zlim=(c(0,30))) #phi and theta change angle0
scatter3D(x=lacdatasm$xcord, y=lacdatasm$ycord, z=lacdatasm$Early.Growth..cm., add = TRUE, type = "h", pch = 19, lwd = 0.2, colkey = FALSE)

persp3D(x = xaxis, y = yaxis, z = z_axisbf, phi = 45, theta = 55, xlab = "x-coord", ylab = "y-coord", zlab = "early growth (cm)", clab = "Phytometer", contour = TRUE, cex.axis = 0.75, cex.lab = 0.75,  zlim=(c(0,30))) #phi and theta change angle
scatter3D(x=lacdatasm$xcord, y=lacdatasm$ycord, z=lacdatasm$Early.Growth..cm., add = TRUE, type = "h", pch = 19, lwd = 0.2, colkey = FALSE)

#let's try a spline-based approach:

library(MBA)
spline <- mba.surf(as.data.frame(cbind(lacdatasm$xcord,lacdatasm$ycord,lacdatasm$Early.Growth..cm.)),31,10)
zcol<-spline$xyz.est
xaxis=seq(0,90, 3)
yaxis=seq(0,45, 5)
persp3D(x = yaxis, y = xaxis, z = zcol$z, phi = 45, theta = 55, xlab = "y-coord", ylab = "x-coord", zlab = "early growth (cm)", clab = "Phytometer", contour = TRUE, cex.axis = 0.75, cex.lab = 0.75,  zlim=(c(0,30))) #phi and theta change angle
scatter3D(y=lacdatasm$xcord, x=lacdatasm$ycord, z=lacdatasm$Early.Growth..cm., add = TRUE, type = "h", pch = 19, lwd = 0.2, colkey = FALSE)

contour2D(x = yaxis, y = xaxis, z = zcol$z, xlab = "y-coord", ylab = "x-coord", clab = "phyt func")
points(y=lacdatasm$xcord, x=lacdatasm$ycord, col=ourcolpal[match(lacdatasm$Early.Growth..cm, ourcolpal),2], cex=0.5, pch=19)


#now for all
splinel <- mba.surf(as.data.frame(cbind(datasm$xcord,datasm$ycord,datasm$Early.Growth..cm.)),31,10)
zcol2<-splinel$xyz.est
persp3D(x = yaxis, y = xaxis, z = zcol2$z, phi = 45, theta = 55, xlab = "y-coord", ylab = "x-coord", zlab = "early growth (cm)", clab = "Phytometer", contour = TRUE, cex.axis = 0.75, cex.lab = 0.75,  zlim=(c(0,30))) #phi and theta change angle
scatter3D(y=datasm$xcord, x=datasm$ycord, z=datasm$Early.Growth..cm., add = TRUE, type = "h", pch = 19, lwd = 0.2, colkey = FALSE)

contour2D(x = yaxis, y = xaxis, z = zcol2$z, xlab = "y-coord", ylab = "x-coord", clab = "phyt func")
points(y=datasm$xcord, x=datasm$ycord, col=ourcolpal2[match(datasm$Early.Growth..cm, ourcolpal2),2], cex=0.5, pch=19)



#now estimating environmental values for the on-lac data points:

xcord<-ycord<-rep(NA, nrow(nolac))

for (i in 1:nrow(nolac)){
  xci <- lettodist(nolac$X[i])[1]
  yci <- lettodist(nolac$X[i])[2]
  xcord[i]<-as.numeric(xci)
  ycord[i]<-as.numeric(yci)
}
nolac<-cbind(nolac,xcord, ycord)

longtab<-matrix(NA, length(zcol$x)*length(zcol$y), 3)
for (i in 1:nrow(zcol$z)){
  yval=zcol$y[i]
  for (j in 1:ncol(zcol$z)){
    xval=zcol$x[j]
    longtab[j+ncol(zcol$z)*(i-1),]<-c(xval, yval, zcol$z[i,j])
  }
}
phytv<-rep(NA, nrow(nolac))
for (i in 1:nrow(nolac)){
  for (j in 1:nrow(longtab)){
    if (longtab[j,1]==nolac$xcord[i]&&longtab[j,2]==nolac$ycord[i]){
      phytv[i]<-longtab[j,3]
    }
  }
}
nolac<-cbind(nolac, phytv)

nldasm<-subset(nolac, xcord<91)
nldasm<-subset(nldasm, ycord<46)

####for data w/lines####:
longtab2<-matrix(NA, length(zcol2$x)*length(zcol2$y), 3)
for (i in 1:nrow(zcol2$z)){
  yval=zcol2$y[i]
  for (j in 1:ncol(zcol2$z)){
    xval=zcol2$x[j]
    longtab2[j+ncol(zcol2$z)*(i-1),]<-c(xval, yval, zcol2$z[i,j])
  }
}
phytv2<-rep(NA, nrow(nolac))
for (i in 1:nrow(nolac)){
  for (j in 1:nrow(longtab2)){
    if (longtab2[j,1]==nolac$xcord[i]&&longtab2[j,2]==nolac$ycord[i]){
      phytv2[i]<-longtab[j,3]
      phyhold<-longtab[j,]
    }
  }
  smatx<-c(nolac$xcord[i],nolac$ycord[i],nolac$Early.Growth..cm.[i],nolac$Total.seed.weight[i], phyhold)
  print(smatx)
}
nolace<-cbind(nolac,phytv2)

plot(nldasm$Early.Growth..cm.~nldasm$phytv, xlab="Model 2 Phytometer Values", ylab="Early Growth")
plot(nldasm$Total.seed.weight~nldasm$Early.Growth..cm., xlab="Early Growth", ylab="Seed Set (g)")
abline(lm(nldasm$Total.seed.weight~nldasm$Early.Growth..cm.))

abline(lm(nldasm$Early.Growth..cm.~nldasm$phytv))
summary(lm(nldasm$Early.Growth..cm.~nldasm$phytv))
summary(lmer(Early.Growth..cm.~phytv+(1|RILID),data=nldasm))
resid(lm(nldasm$Early.Growth..cm.~nldasm$phytv)) 
plot(nldasm$Total.seed.weight~nldasm$phytv,xlab="Model 2 Phytometer Values", ylab="Seed Set (g)")
abline(lm(nldasm$Total.seed.weight~nldasm$phytv))
summary(lmer(Total.seed.weight~phytv+(1|RILID),data=nldasm))
invratio<-nldasm$phytv/mean(nldasm$phytv, na.rm=T)
adjweight<-round(nldasm$Total.seed.weight/invratio,5)
nldas<-cbind(nldasm,adjweight)
plot(nldas$adjweight~nldas$phytv)
abline(lm(nldas$adjweight~nldas$phytv))
summary(lmer(adjweight~phytv+(1|RILID),data=nldas))
summary(lmer(adjweight~RILcolor+(1|RILID),data=nldas))
#even this version isn't a good predictor of early growth/seed set


egest00<-growthre(nolac$xcord, nolac$ycord)
egest<-growthrebf(nolac$xcord, nolac$ycord)
egest03<-growthre3(nolac$xcord, nolac$ycord)
egest04<-growthre4(nolac$xcord, nolac$ycord)


egestx<-growthresm(nolac$xcord, nolac$ycord)
egest2<-growthrebfsm(nolac$xcord, nolac$ycord) #rebfsm is from quadratic
egest2m<-growthrebfsm(nldasm$xcord, nldasm$ycord)
egest3<-growthresm3(nolac$xcord, nolac$ycord)
egest4<-growthresm4(nolac$xcord, nolac$ycord)

egestxx<-gresm(nolac$xcord, nolac$ycord)
egestbf<-grebfsm(nolac$xcord, nolac$ycord)
egest33<-gresm3(nolac$xcord, nolac$ycord)
egest44<-gresm4(nolac$xcord, nolac$ycord)

nolacs<-cbind(nolac,egest2, egestbf)
nld<-cbind(nldasm,egest2m)

plot(nld$Early.Growth..cm.~nld$egest2m)
abline(lm(nld$Early.Growth..cm.~nld$egest2m))
summary(lm(nld$Early.Growth..cm.~nld$egest2m))

resid(lm(nld$Early.Growth..cm.~nld$egest2))#appears to be a good predictor (R-squared is 0.08136)


plot(nld$Total.seed.weight~nld$egest2m, xlab="Model 1 Phytometer Values", ylab="Seed set (g)")
abline(lm(nld$Total.seed.weight~nld$egest2m))
summary(lm(nld$Total.seed.weight~nld$egest2m))
plot(nld$Projected..~nld$egest2m, xlim=c(13,23))
abline(lm(nld$Projected..~nld$egest2m))
summary(lm(nld$Projected..~nld$egest2m))
segr<-resid(lm(nld$Total.seed.weight~nld$egest2m)) #less of a good fit, bit estimates do have some predictive power in seed set

nolacwseed<-subset(nld, as.numeric(nld$Total.seed.weight)>=0)
plot(nolacwseed$phytv~nolacwseed$egest2m)
abline(lm(nolacwseed$phytv~nolacwseed$egest2m))
summary(lm(nolacwseed$phytv~nolacwseed$egest2m))

nolacseeds<-cbind(nolacwseed,segr)

#adjusted seed weight
seedweightadj<-nolacseeds$Total.seed.weight*ave(nolacseeds$egest2m)/nolacseeds$egest2m
seednumbadj<-nolacseeds$Projected..*ave(nolacseeds$egest2m)/nolacseeds$egest2m
##
nolacseeds<-cbind(nolacseeds,seedweightadj, seednumbadj)
plot(nolacseeds$seedweightadj~nolacseeds$egest2m)
abline(lm(nolacseeds$seedweightadj~nolacseeds$egest2m))
plot(nolacseeds$seednumbadj~nolacseeds$egest2m)
abline(lm(nolacseeds$seednumbadj~nolacseeds$egest2m))
summary(lm(nolacseeds$seedweightadj~nolacseeds$egest2m))
summary(lm(nolacseeds$seednumbadj~nolacseeds$egest2m))
summary(lmer(seedweightadj~egest2m+RILcolor+(1|RILID),data=nolacseeds))
summary(lmer(seedweightadj~RILcolor+(1|RILID),data=nolacseeds))
summary(lmer(segr~RILcolor+(1|RILID),data=nolacseeds))
summary(lmer(seednumbadj~egest2m+RILcolor+(1|RILID),data=nolacseeds))
summary(lmer(seednumbadj~RILcolor+(1|RILID),data=nolacseeds))
plot(nolacseeds$Early.Growth..cm.~nolacseeds$egest2m)
abline(lm(nolacseeds$Early.Growth..cm.~nolacseeds$egest2m))
summary(lm(nolacseeds$Early.Growth..cm.~nolacseeds$egest2m))
summary(lmer(Early.Growth..cm.~egest2m+(1|RILID),data=nolacseeds))
#strangely, there is a negative correlation with early growth predictions based off of lac and seed set
#early growth however was predicted by seed set


pinds<-subset(nolacseeds, nolac$RILcolor=="P")
winds<-subset(nolacseeds, nolac$RILcolor=="W")

##################################
#make a line based grid for selection analysis by averaging seed weights and resids
linenumb<-nlevels(as.factor(nolacseeds$RILID))
linematrix<-matrix(data=NA, nrow=linenumb, ncol=10)
xystore<-matrix(data=NA, nrow=linenumb, ncol=3)
for (i in 1:linenumb){
  templine<-subset(nolacseeds, as.numeric(as.factor(nolacseeds$RILID))==i)
  stemp<-mean(templine$Total.seed.weight)
  ptemp<-mean(templine$Projected..)
  rtemp<-mean(templine$segr)
  adjtemp<-mean(templine$seedweightadj)
  adjtemp2<-mean(templine$seednumbadj)
  fztemp<-mean(templine$Cor_Len_Var)
  sctemp<-mean(templine$Sugar_UND)
  xcordm<-mean(templine$xcord)
  ycordm<-mean(templine$ycord)
  sline<-1-mean(as.numeric(templine$Dead.), na.rm=T)
  coltemp<-templine$RILcolor[1]
  linehold<-templine$RILID[1]
  linedata<-c(linehold, coltemp,as.numeric(stemp),as.numeric(ptemp),as.numeric(rtemp),as.numeric(adjtemp),as.numeric(adjtemp2),as.numeric(fztemp),as.numeric(sctemp), as.numeric(sline))
  linematrix[i,]<-linedata
  xystore[i,]<-c(as.numeric(linehold),xcordm,ycordm)
}
linematrix<-as.data.frame(linematrix)
labelen<- c("Line", "Color", "Seed_Set","Projected_Number", "Resids", "Adj_seed_weight", "Adj_number", "Flower_Size", "Sugar_Conc", "Survival")
colnames(linematrix)<-labelen
linematrix$Seed_Set<-as.numeric(linematrix$Seed_Set)
linematrix$Projected_Number<-as.numeric(linematrix$Projected_Number)
linematrix$Resids<-as.numeric(linematrix$Resids)
linematrix$Adj_seed_weight<-as.numeric(linematrix$Adj_seed_weight)
linematrix$Adj_number<-as.numeric(linematrix$Adj_number)
linematrix$Flower_Size<-as.numeric(linematrix$Flower_Size)
linematrix$Sugar_Conc<-as.numeric(linematrix$Sugar_Conc)
linematrix$Survival<-as.numeric(linematrix$Survival)

#######now for selection models
#with seed set
#first as a lm, then as a glm
selcsc<-lm(Seed_Set~Color+Color*Sugar_Conc, data=linematrix)
selcfz<-lm(Seed_Set~Color+Color*Flower_Size, data=linematrix)
selscfz<-lm(Seed_Set~Flower_Size*Sugar_Conc, data=linematrix)
selcscp<-lm(Projected_Number~Color+Color*Sugar_Conc, data=linematrix)
selcfzp<-lm(Projected_Number~Color*Flower_Size, data=linematrix)
selscfzp<-lm(Projected_Number~Flower_Size*Sugar_Conc, data=linematrix)
sel3w<-lm(Seed_Set~Color+Color*Flower_Size*Sugar_Conc, data=linematrix)
summary(selscfzp)

miw<-as.matrix(linematrix)
sel3wg<-glm(Seed_Set~Color+Color*Flower_Size+Color*Sugar_Conc, data=linematrix, family=Gamma)
sel3wgp<-glm(round(Projected_Number)~Color+Color*Flower_Size+Color*Sugar_Conc, data=linematrix, family=poisson)

#not working b/c of list read error

#with residuals
sel3rs<-lm(Resids~Color+Color*Flower_Size*Sugar_Conc, data=linematrix)
selcscrs<-lm(Resids~Color+Color*Sugar_Conc, data=linematrix)
selcfzrs<-lm(Resids~Color+Color*Flower_Size, data=linematrix)
summary(selcfzrs)
boxplot(Resids~Color,data=linematrix)
boxplot(Adj_seed_weight~Color,data=linematrix)
boxplot(Adj_number~Color,data=linematrix)

#with adjusted:
asel3<-lm(Adj_seed_weight~Color+Color*Flower_Size*Sugar_Conc, data=linematrix)
aselcsc<-lm(Adj_seed_weight~Color+Color*Sugar_Conc, data=linematrix)
aselcfz<-lm(Adj_seed_weight~Color+Color*Flower_Size, data=linematrix)
summary(aselcfz)
boxplot(Adj_seed_weight~Color,data=linematrix)
summary(lm(Adj_seed_weight~Color, data=linematrix))

asel3<-lm(Adj_number~Color+Color*Flower_Size*Sugar_Conc, data=linematrix)
aselcsc<-lm(Adj_number~Color+Color*Sugar_Conc, data=linematrix)
aselcfz<-lm(Adj_number~Color+Color*Flower_Size, data=linematrix)
summary(aselcfz)
boxplot(Adj_number~Color,data=linematrix)
summary(lm(Adj_number~Color, data=linematrix))

######################no more lines until later

#let's do this piecemeal:
#first, are the possible predictors correlated?
#we can take the most complicated glmmodel and just look at the correlational matrix
#anything over 0.4 should be under consideration as "highly correlated
#we can also assume that certain traits that irene showed to be highly coreelated re still highly correlated
summary(lmer(Sugar_UND~Nectar_UNC+(1|RILID), data=nolac))
plot(Sugar_UND~Nectar_UNC, data=nolac)
abline(lm(Sugar_UND~Nectar_UNC, data=nolac)) #adjusted R-squared is 0.2
#sugar and nectar, for example, are highly correlated
summary(lmer(Sugar_UND~Sugar_D.2.5+(1|RILID), data=nolac))
plot(Sugar_UND~Sugar_D.2.5, data=nolac)
abline(lm(Sugar_UND~Sugar_D.2.5, data=nolac)) #adjusted R-squared is 0.42
#sugar measures are highly correlated
summary(lmer(Seed.Mass~Seed.Length+(1|RILID), data=nolac))
plot(Seed.Mass~Seed.Length, data=nolac)
abline(lm(Seed.Mass~Seed.Length, data=nolac))#adjusted R-squared is 0.41
#seed mass and seed length highly covary

#here is the covariance matrix for most traits:
saxx<-glmer(Total.seed.weight~RILcolor+Cor_Len_Var+Seed.Mass+Sugar_UND+(1|blockID)+(1|RILID), family=Gamma, data=nolac)
saxb<-glmer(Total.seed.weight~blockID+RILcolor+Cor_Len_Var+Seed.Mass+Sugar_UND+(1|RILID), family=Gamma, data=nolac)
summary(saxx)
summary(saxb)
saxy<-glmer(Projected..~RILcolor+Cor_Len_Var+Seed.Mass+Sugar_UND+(1|blockID)+(1|RILID), family=poisson, data=nolac)
summary(saxy)
#all are under 0.4, but let's test some of the covariances
summary(lmer(Seed.Mass~Cor_Len_Var+(1|RILID), data=nolac))
plot(Seed.Mass~Cor_Len_Var, data=nolac)
abline(lm(Seed.Mass~Cor_Len_Var, data=nolac))#adjusted R-squared is 0.0122; probably pretty negligible covariance
summary(lmer(Seed.Mass~Sugar_UND+(1|RILID), data=nolac))
plot(Seed.Mass~Sugar_UND, data=nolac)
abline(lm(Seed.Mass~Sugar_UND, data=nolac))#adjusted R-squared is 0.004; probably pretty negligible covariance
summary(lmer(Seed.Mass~RILcolor+(1|RILID), data=nolac))
boxplot(Seed.Mass~RILcolor, data=nolac)
abline(lm(Seed.Mass~RILcolor, data=nolac))#highly sig, but adj R squared is 0.045
summary(lmer(Sugar_UND~RILcolor+(1|RILID), data=nolac))
boxplot(Sugar_UND~RILcolor, data=nolac)
abline(lm(Sugar_UND~RILcolor, data=nolac))#adj R squared is near 0, -0.00097
summary(lmer(Cor_Len_Var~RILcolor+(1|RILID), data=nolac))
boxplot(Cor_Len_Var~RILcolor, data=nolac)
pnoted<-subset(nolac, nolac$Obs.Flower.Color=="P")
wnoted<-subset(nolac, nolac$Obs.Flower.Color=="W")
tnoted<-rbind(pnoted,wnoted)
boxplot(Cor_Len_Var~Obs.Flower.Color, data=tnoted)
summary(lmer(Cor_Len_Var~Obs.Flower.Color+(1|RILID), data=tnoted))
boxplot(Sugar_UND~Obs.Flower.Color, data=tnoted)
summary(lmer(Sugar_UND~Obs.Flower.Color+(1|RILID), data=tnoted))
abline(lm(Cor_Len_Var~RILcolor, data=nolac))#highly sig, but adj R squared is 0.02862
summary(lmer(Cor_Len_Var~Sugar_UND+(1|RILID), data=nolac))
plot(Cor_Len_Var~Sugar_UND, data=nolac, xlab="Sugar Concentration", ylab="Flower Size")
abline(lm(Cor_Len_Var~Sugar_UND, data=nolac))#adjusted R-squared of 0.103
#seems like the only clear correlation is btwn sugar conc and flower size, maybecolor with these traits
#to test this, we can look at glmer models with single effects for each factor and joint effect models
seedwtc<-glmer(Total.seed.weight~RILcolor+(1|blockID)+(1|RILID), family=Gamma, data=nolac)
seedwtco<-glmer(Total.seed.weight~Obs.Flower.Color+(1|blockID)+(1|RILID), family=Gamma, data=tnoted)
seedwtcb<-glmer(Total.seed.weight~blockID+RILcolor+(1|RILID), family=Gamma, data=nolac)
seedwtcbo<-glmer(Total.seed.weight~blockID+Obs.Flower.Color+(1|RILID), family=Gamma, data=tnoted)


seedwtcj<-glmer(seedweightadj~Early.Growth..cm.+RILcolor+(1|blockID)+(1|RILID), family=Gamma, data=nolacseeds)
summary(seedwtcb)
seedwtcp<-glmer(Projected..~RILcolor+(1|blockID)+(1|RILID), family=poisson, data=nolac)
summary(seedwtcp)
seedwtcpb<-glmer(Projected..~blockID+RILcolor+(1|RILID), family=poisson, data=nolac)
summary(seedwtcpb)
oseedwtcp<-glmer(Projected..~Obs.Flower.Color+(1|blockID)+(1|RILID), family=poisson, data=tnoted)
summary(oseedwtcp) #this is throwing an error
oseedwtcpb<-glmer(Projected..~blockID+Obs.Flower.Color+(1|RILID), family=poisson, data=tnoted)
summary(oseedwtcpb)

boxplot(Total.seed.weight~RILcolor,data=nolac, col=c("magenta","white"), xlab="Line Color", ylab="Seed Set(g)")
boxplot(Total.seed.weight~Obs.Flower.Color,data=tnoted, xlab="Line Color", ylab="Seed Set(g)")
boxplot(Projected..~RILcolor,data=nolac, col=c("magenta","white"), xlab="Line Color", ylab="Estimated Number of Seeds")
boxplot(Projected..~Obs.Flower.Color,data=tnoted,  xlab="Line Color", ylab="Estimated Number of Seeds")

Anova(seedwtc)
Anova(seedwtcp)
Anova(seedwtco)
Anova(oseedwtcp) #this is off

Anova(lmer(Total.seed.weight~RILcolor+(1|blockID)+(1|RILID), data=nolac))
Anova(lmer(Projected..~RILcolor+(1|blockID)+(1|RILID), data=nolac))
Anova(lmer(Total.seed.weight~Early.Growth..cm.+(1|RILID), data=nolac))
Anova(lmer(seedweightadj~Early.Growth..cm.+RILcolor+(1|blockID)+(1|RILID), data=nolacseeds))
Anova(lmer(segr~Early.Growth..cm.+RILcolor+(1|blockID)+(1|RILID), data=nolacseeds))

Anova(lmer(Total.seed.weight~Obs.Flower.Color+(1|blockID)+(1|RILID), data=tnoted))
Anova(lmer(Total.seed.weight~Early.Growth..cm.+(1|RILID), data=tnoted))

summary(aov(Total.seed.weight~Error(RILID)+Early.Growth..cm.+Cor_Len_Var*Sugar_UND, data=nolacseeds))
summary(aov(seedweightadj~Error(RILID)+Early.Growth..cm.+Cor_Len_Var*Sugar_UND, data=nolacseeds))
summary(aov(segr~Error(RILID)+Early.Growth..cm.+Cor_Len_Var*Sugar_UND, data=nolacseeds))

summary(aov(Total.seed.weight~Error(RILID)+Cor_Len_Var*Sugar_UND, data=nolacseeds))
summary(aov(seedweightadj~Error(RILID)+Cor_Len_Var*Sugar_UND, data=nolacseeds))
summary(aov(segr~Error(RILID)+Cor_Len_Var*Sugar_UND, data=nolacseeds))
#adjusting for early growth or using residuals both provide more ability to pick up significat trends

#color not sig predictor
seddsg<-glmer(Total.seed.weight~Early.Growth..cm.+Sugar_UND+(1|blockID)+(1|RILID), family=Gamma, data=nolacseeds)
summary(seddsg) 
Anova(seddsg)
seddsgp<-glmer(Projected..~Sugar_UND+(1|blockID)+(1|RILID), family=poisson, data=nolacseeds)
summary(seddsgp) 
Anova(seddsgp)

#p-value is .108 for just sugar for weight
#p-value is .0703 for estimated number
seddsgj<-glmer(seedweightadj~Early.Growth..cm.+Sugar_UND+(1|blockID)+(1|RILID), family=Gamma, data=nolacseeds)
summary(seddsgj) 
Anova(seddsgj) #corrects down to suggestive with adjusted for early growth

seedsz<-glmer(Total.seed.weight~Cor_Len_Var+(1|blockID)+(1|RILID), family=Gamma, data=nolacseeds)
summary(seedsz) 
Anova(seedsz)
seedszp<-glmer(Projected..~Cor_Len_Var+(1|blockID)+(1|RILID), family=poisson, data=nolacseeds)
summary(seedszp) 
Anova(seedszp)
#p-value is .101 for size
#Anova is .08613 for projected

seedszj<-glmer(seedweightadj~Early.Growth..cm.+Cor_Len_Var+(1|blockID)+(1|RILID), family=Gamma, data=nolacseeds)
summary(seedszj) 
Anova(seedszj)
#but p-value goes up with size to 0.135

seedsm<-glmer(Total.seed.weight~Early.Growth..cm.+Seed.Mass+(1|blockID)+(1|RILID), family=Gamma(link = "log"), data=nolac)
summary(seedsm)
Anova(seedsm)
seedsm2<-glmer(Total.seed.weight~Early.Growth..cm.+X50.seed.weight+(1|blockID)+(1|RILID), family=Gamma(link = "log"), data=nolac)
summary(seedsm2)
Anova(seedsm2)

seedsmj<-glmer(seedweightadj~Early.Growth..cm.+Seed.Mass+(1|blockID)+(1|RILID), family=Gamma(link = "log"), data=nolacseeds)

#seed mass not sig predictor
#not let's look at some joint predictions
#first one that we know is covaried 0 size and sugar conc
seedszg<-glmer(Total.seed.weight~Cor_Len_Var*Sugar_UND+(1|blockID)+(1|RILID), family=Gamma(link = "log"), data=nolacseeds)
summary(seedszg)
Anova(seedszg)
seedszgp<-glmer(Projected..~Cor_Len_Var*Sugar_UND+(1|blockID)+(1|RILID), family=poisson, data=nolacseeds)
summary(seedszgp)
Anova(seedszgp)

seedszgj<-glmer(seedweightadj~Early.Growth..cm.+Cor_Len_Var*Sugar_UND+(1|blockID)+(1|RILID), family=Gamma(link = "log"), data=nolacseeds)
summary(seedszgj)
Anova(seedszgj)

seedszng<-glmer(Total.seed.weight~Cor_Len_Var*Sugar_UND+(1|blockID)+(1|RILID), family=Gamma(link = "log"), data=nolacseeds)
seedszngj<-glmer(seedweightadj~Cor_Len_Var*Sugar_UND+(1|blockID)+(1|RILID), family=Gamma(link = "log"), data=nolacseeds)


Anova(lmer(Total.seed.weight~Early.Growth..cm.+Cor_Len_Var*Sugar_UND+(1|blockID)+(1|RILID), data=nolacseeds))
Anova(lmer(seedweightadj~Early.Growth..cm.+Cor_Len_Var*Sugar_UND+(1|blockID)+(1|RILID), data=nolacseeds))
Anova(lmer(segr~Early.Growth..cm.+Cor_Len_Var*Sugar_UND+(1|blockID)+(1|RILID), data=nolacseeds))

Anova(lmer(Total.seed.weight~Cor_Len_Var*Sugar_UND+(1|blockID)+(1|RILID), data=nolacseeds))
Anova(lmer(seedweightadj~Cor_Len_Var*Sugar_UND+(1|blockID)+(1|RILID), data=nolacseeds))
Anova(lmer(segr~Cor_Len_Var*Sugar_UND+(1|blockID)+(1|RILID), data=nolacseeds))

summary(aov(Total.seed.weight~Error(RILID)+Early.Growth..cm.+Cor_Len_Var*Sugar_UND, data=nolacseeds))
summary(aov(seedweightadj~Error(RILID)+Early.Growth..cm.+Cor_Len_Var*Sugar_UND, data=nolacseeds))
summary(aov(segr~Error(RILID)+Early.Growth..cm.+Cor_Len_Var*Sugar_UND, data=nolacseeds))

summary(aov(Total.seed.weight~Error(RILID)+Cor_Len_Var*Sugar_UND, data=nolacseeds))
summary(aov(seedweightadj~Error(RILID)+Cor_Len_Var*Sugar_UND, data=nolacseeds))
summary(aov(segr~Error(RILID)+Cor_Len_Var*Sugar_UND, data=nolacseeds))


#now everything is sig, including predictors that weren't individually
seddsgc<-glmer(Total.seed.weight~Sugar_UND*RILcolor+(1|blockID)+(1|RILID), family=Gamma, data=nolacseeds)
seddsgcj<-glmer(seedweightadj~Sugar_UND*RILcolor+(1|blockID)+(1|RILID), family=Gamma, data=nolacseeds)
summary(seddsgc)
Anova(seddsgc)
seddsgcp<-glmer(Projected..~Sugar_UND*RILcolor+(1|blockID)+(1|RILID), family=poisson, data=nolacseeds)
summary(seddsgcp)
Anova(seddsgcp)

Anova(seddsgcj)
#only sugar is suggestive in both (mixed is close to suggestiveness though)
#for projected number, this pattern is more striking
seedcsz<-glmer(Total.seed.weight~Early.Growth..cm.+RILcolor+Cor_Len_Var*RILcolor+(1|blockID)+(1|RILID), family=Gamma, data=nolacseeds)
seedcszj<-glmer(seedweightadj~Early.Growth..cm.+RILcolor+Cor_Len_Var*RILcolor+(1|blockID)+(1|RILID), family=Gamma, data=nolacseeds)
summary(seedcsz)
Anova(seedcsz)
seddcszp<-glmer(Projected..~Cor_Len_Var*RILcolor+(1|blockID)+(1|RILID), family=poisson, data=nolacseeds)
summary(seddcszp)
Anova(seddcszp)

Anova(seedcszj)
#nothing is sig nor sugg with correction?
#for Poisson size itself is suggestive
seedcsm<-glmer(Total.seed.weight~Early.Growth..cm.+Seed.Mass*RILcolor+(1|blockID)+(1|RILID), family=Gamma(link = "log"), data=nolacseeds)
seedcsmj<-glmer(seedweightadj~Early.Growth..cm.+Seed.Mass*RILcolor+(1|blockID)+(1|RILID), family=Gamma(link = "log"), data=nolacseeds)
summary(seedcsm)
Anova(seedcsm)
#seed mass no sig predictors
#now let's do a 3 way model
seed3<-glmer(Total.seed.weight~Early.Growth..cm.+Cor_Len_Var*Sugar_UND*RILcolor+(1|blockID)+(1|RILID), family=Gamma(link = "log"), data=nolacseeds)
seed32<-glmer(Total.seed.weight~Cor_Len_Var*Sugar_UND+Sugar_UND*RILcolor+Cor_Len_Var*RILcolor+(1|blockID)+(1|RILID), family=Gamma(link = "log"), data=nolacseeds)
seed32b<-glmer(Total.seed.weight~blockID+Cor_Len_Var*Sugar_UND+Sugar_UND*RILcolor+Cor_Len_Var*RILcolor+(1|RILID), family=Gamma(link = "log"), data=nolacseeds)

seed32p<-glmer(Projected..~Cor_Len_Var*Sugar_UND+Sugar_UND*RILcolor+Cor_Len_Var*RILcolor+(1|blockID)+(1|RILID), family=poisson, data=nolacseeds)
seed32pb<-glmer(Projected..~blockID+Cor_Len_Var*Sugar_UND+Sugar_UND*RILcolor+Cor_Len_Var*RILcolor+(1|RILID), family=poisson, data=nolacseeds)



#without early growth, this model seems broken (all highly sig, p<10^-16)
seed3j<-glmer(seedweightadj~Early.Growth..cm.+Cor_Len_Var*Sugar_UND*RILcolor+(1|blockID)+(1|RILID), family=Gamma(link = "log"), data=nolacseeds)
seed3r<-lmer(segr~Early.Growth..cm.+Cor_Len_Var*Sugar_UND*RILcolor+(1|blockID)+(1|RILID), data=nolacseeds)

summary(seed3)
summary(seed3j)
summary(seed3r)
Anova(seed3)
Anova(seed3j)
Anova(seed3r)

#interestingly with adjusted seed weights, only sugar is sig in 3-way model

#seed3 ancova
seed3anc<-aov(Total.seed.weight~Cor_Len_Var*Sugar_UND*RILcolor+Error(RILID), data=nolac)
seed3ancp<-aov(Projected..~Cor_Len_Var*Sugar_UND*RILcolor+Error(RILID), data=nolac)


#now the 3 way corr, 2 way sug and color are sig, just color is sug as well as colorxsize
seedxplct<-glmer(Total.seed.weight~Early.Growth..cm.+Herbiv+Cor_Len_Var+Sugar_UND+RILcolor+Seed.Mass+Seed.Mass*RILcolor+Cor_Len_Var*Sugar_UND+Sugar_UND*RILcolor+Cor_Len_Var*RILcolor+Cor_Len_Var*Sugar_UND*RILcolor+(1|blockID)+(1|RILID), family=Gamma(link = "log"), data=nolac)
seed0<-glmer(Total.seed.weight~1+(1|blockID)+(1|RILID), family=Gamma(link = "log"), data=nolac)
seedeg<-glmer(Total.seed.weight~Early.Growth..cm.+(1|blockID)+(1|RILID), family=Gamma(link = "log"), data=nolac)
seedsmax<-glmer(Total.seed.weight~Early.Growth..cm.+Seed.Mass*RILcolor+Cor_Len_Var*Sugar_UND+(1|blockID)+(1|RILID), family=Gamma(link = "log"), data=nolac)

seedxplctj<-glmer(seedweightadj~Early.Growth..cm.+Herbiv+Cor_Len_Var+Sugar_UND+RILcolor+Seed.Mass+Seed.Mass*RILcolor+Cor_Len_Var*Sugar_UND+Sugar_UND*RILcolor+Cor_Len_Var*RILcolor+Cor_Len_Var*Sugar_UND*RILcolor+(1|blockID)+(1|RILID), family=Gamma(link = "log"), data=nolacseeds)
seed0j<-glmer(seedweightadj~1+(1|blockID)+(1|RILID), family=Gamma(link = "log"), data=nolacseeds)
seedegj<-glmer(seedweightadj~Early.Growth..cm.+(1|blockID)+(1|RILID), family=Gamma(link = "log"), data=nolacseeds)
seedsmaxj<-glmer(seedweightadj~Early.Growth..cm.+Seed.Mass*RILcolor+Cor_Len_Var*Sugar_UND+(1|blockID)+(1|RILID), family=Gamma(link = "log"), data=nolacseeds)


library(bbmle)
bbmle::ICtab(seedxplct,seedsmax,seedcsm,seedwtc,seedsm, seedeg, seed0, base = T)#model selection
#seed color x mass
bbmle::ICtab(seedxplct,seed3,seedcsz,seedwtc,seedsz, seedeg,seed0, base = T)#model selection
#seed color x size
bbmle::ICtab(seedxplct,seed3,seddsgc,seedwtc,seddsg, seedeg,seed0, base = T)#model selection
#seed color X sugar
bbmle::ICtab(seedxplct,seed3,seedszg,seedsz,seddsg,seedeg, seed0, base = T)#model selection

bbmle::ICtab(seedxplct,seed3,seedsmax,seedcsm, seedsm, seedcsz,seedwtc,seedszg,seedsz,seddsg,seedeg, seed0, base = T)#model selection

#for adjusted
bbmle::ICtab(seedxplctj,seedsmaxj,seedcsmj,seedwtcj,seedsmj, seedegj, seed0j, base = T)#model selection
#seed color x mass
bbmle::ICtab(seedxplctj,seed3j,seedcszj,seedwtcj,seedszj, seedegj,seed0j, base = T)#model selection
#seed color x size
bbmle::ICtab(seedxplctj,seed3j,seddsgcj,seedwtcj,seddsgj, seedegj,seed0j, base = T)#model selection
#seed color X sugar
bbmle::ICtab(seedxplctj,seed3j,seedszgj,seedszj,seddsgj,seedegj, seed0j, base = T)#model selection

bbmle::ICtab(seedxplct,seed3,seedsmax,seedcsm, seedsm, seedcsz,seedwtc,seedszg,seedsz,seddsg,seedeg, seed0, base = T)#model selection


summary(aov(nolac$Herbiv~nolac$RILcolor))
summary(aov(nolac$Herbiv~Error(nolac$RILID)))
herbm<-glmer(Herbiv~RILcolor+Cor_Len_Var*RILcolor+Seed.Mass*RILcolor+Sugar_UND*RILcolor+(1|blockID)+(1|RILID), family=poisson, data=nolac)
herbc<-glmer(Herbiv~RILcolor+(1|blockID)+(1|RILID), family=poisson, data=nolac)
herbcsz<-glmer(Herbiv~Cor_Len_Var*RILcolor+(1|blockID)+(1|RILID), family=poisson, data=nolac)
herbsz<-glmer(Herbiv~Cor_Len_Var+(1|blockID)+(1|RILID), family=poisson, data=nolac)
herbcsg<-glmer(Herbiv~Sugar_UND*RILcolor+(1|blockID)+(1|RILID), family=poisson, data=nolac)
#maybe followup with color x sugar effects on herbivory?

#now for early growth
summary(aov(nolacseeds$Early.Growth..cm.~nolacseeds$RILcolor))
summary(aov(nolacseeds$Early.Growth..cm.~Error(nolacseeds$RILID)))

egm<-glmer(Early.Growth..cm.~RILcolor+Cor_Len_Var*RILcolor+Seed.Mass*RILcolor+Sugar_UND*RILcolor+(1|blockID)+(1|RILID), family=Gamma(link = "log"), data=nolacseeds)
egc<-glmer(Early.Growth..cm.~RILcolor+(1|blockID)+(1|RILID), family=Gamma(link = "log"), data=nolacseeds)
egcsz<-glmer(Early.Growth..cm.~Cor_Len_Var*RILcolor+(1|blockID)+(1|RILID), family=Gamma(link = "log"), data=nolacseeds)
egsz<-glmer(Early.Growth..cm.~Cor_Len_Var+(1|blockID)+(1|RILID), family=Gamma(link = "log"), data=nolacseeds)
egcsg<-glmer(Early.Growth..cm.~Sugar_UND*RILcolor+(1|blockID)+(1|RILID), family=Gamma(link = "log"), data=nolacseeds)
#no predictive power of color or size on early growth

#testing seed weight differences w/line
boxplot(X50.seed.weight~RILcolor, data=nolacseeds, ylab="Mass of 50 seeds (g)")
summary(lm(X50.seed.weight~RILcolor,data=nolacseeds))
Anova(lmer(X50.seed.weight~RILcolor+(1|RILID)+(1|blockID),data=nolacseeds))
summary(lmer(X50.seed.weight~RILcolor+(RILcolor|RILID)+(1|blockID),data=nolacseeds))
summary(lmer(X50.seed.weight~blockID+RILcolor+(1|RILID),data=nolacseeds))

#random effect of RIL makes result n.s.
shapiro.test(resid(aov(X50.seed.weight~RILcolor,data=nolacseeds))) #normal enough to avoid glm


plot(X50.seed.weight~Seed.Mass, data=nolacseeds, ylab="")
abline(lm(X50.seed.weight~Seed.Mass, data=nolacseeds))
summary(lm(X50.seed.weight~Seed.Mass, data=nolacseeds))

#adjusting seed set to exclude 0 entry- change to .0000001

seedwt<-glmer(Total.seed.weight~RILcolor+Cor_Len_Var*RILcolor+Seed.Mass*RILcolor+Sugar_UND*RILcolor+(1|blockID)+(1|RILID), family=Gamma, data=nolac)
seddwtx<-glmer(Total.seed.weight~Early.Growth..cm.+RILcolor+Cor_Len_Var*RILcolor+Seed.Mass*RILcolor+Sugar_UND*RILcolor+(1|blockID)+(1|RILID), family=Gamma, data=nolac)
#need to test if gamma dist is appropriate over multiplying x1000, rounding and testing as poisson
#plotting residuals:
reseed<-resid(seddwtx)
plot(fitted(seddwtx,reseed))
qqnorm(reseed)
qqline(reseed)
plot(density(reseed))

#################################
#I think this next section is just repetitive code
#now seed weight



#both seem to indicate that color and sugar conc are sig predictors,let's run individual models:
seedwtc<-glmer(Total.seed.weight~Early.Growth..cm.+RILcolor+(1|blockID)+(1|RILID), family=Gamma, data=nolac)
seedwtcegr<-glmer(Total.seed.weight~(1|Early.Growth..cm.)+RILcolor+(1|blockID)+(1|RILID), family=Gamma, data=nolac)
#note that color not a sign predictor here
seddsgc<-glmer(Total.seed.weight~Early.Growth..cm.+Sugar_UND*RILcolor+(1|blockID)+(1|RILID), family=Gamma, data=nolac)
seddsgcegr<-glmer(Total.seed.weight~(1|Early.Growth..cm.)+Sugar_UND*RILcolor+(1|blockID)+(1|RILID), family=Gamma, data=nolac)
#sugar sig, color and sugxcolor hang around 0.1-0.15
seddsg<-glmer(Total.seed.weight~Early.Growth..cm.+Sugar_UND+(1|blockID)+(1|RILID), family=Gamma, data=nolac)
seddsgegr<-glmer(Total.seed.weight~(1|Early.Growth..cm.)+Sugar_UND+(1|blockID)+(1|RILID), family=Gamma, data=nolac)
#sugar just suggestive without color
seedcsz<-glmer(Total.seed.weight~Early.Growth..cm.+RILcolor+Cor_Len_Var*RILcolor+(1|blockID)+(1|RILID), family=Gamma, data=nolac)
seedcszegr<-glmer(Total.seed.weight~(1|Early.Growth..cm.)+RILcolor+Cor_Len_Var*RILcolor+(1|blockID)+(1|RILID), family=Gamma, data=nolac)
#size is sig, and color, color*size is sug
seedsz<-glmer(Total.seed.weight~Early.Growth..cm.+Cor_Len_Var+(1|blockID)+(1|RILID), family=Gamma, data=nolac)
seedszegr<-glmer(Total.seed.weight~(1|Early.Growth..cm.)+Cor_Len_Var+(1|blockID)+(1|RILID), family=Gamma, data=nolac)

#3wayinteractions
seed3w<-glmer(Total.seed.weight~Early.Growth..cm.+RILcolor+Cor_Len_Var*RILcolor+Sugar_UND*RILcolor+(1|blockID)+(1|RILID), family=Gamma, data=nolac)
seed3wegr<-glmer(Total.seed.weight~(1|Early.Growth..cm.)+RILcolor+Cor_Len_Var*RILcolor+Sugar_UND*RILcolor+(1|blockID)+(1|RILID), family=Gamma, data=nolac)
#all 3are sugg now - what does all of this mean


#lets try with offset for early growth
toff<-glmer(Total.seed.weight~RILcolor+Cor_Len_Var*RILcolor+Seed.Mass*RILcolor+Sugar_UND*RILcolor+(1|blockID)+(1|RILID)+offset(log(Early.Growth..cm.)), family=Gamma, data=nolac)



#need to test if gamma dist is appropriate over multiplying x1000, rounding and testing as poisson
#plotting residuals:
reseed<-resid(seddwtx)
plot(fitted(seddwtx,reseed))
qqnorm(reseed)
qqline(reseed)
plot(density(reseed))

#now as a lambda approximation:
nolac2<-nolac
nolac2$Total.seed.weight<-round(nolac2$Total.seed.weight*10) #weight in mg
seddwtx2<-glmer(Total.seed.weight~Early.Growth..cm.+RILcolor+Cor_Len_Var*RILcolor+Seed.Mass*RILcolor+Sugar_UND*RILcolor+(1|blockID)+(1|RILID), family=poisson, data=nolac2)
seddwtx2egr<-glmer(Total.seed.weight~(1|Early.Growth..cm.)+RILcolor+Cor_Len_Var*RILcolor+Seed.Mass*RILcolor+Sugar_UND*RILcolor+(1|blockID)+(1|RILID), family=poisson, data=nolac2)
seddnsz2egr<-glmer(Total.seed.weight~(1|Early.Growth..cm.)+RILcolor+Seed.Mass*RILcolor+Sugar_UND*RILcolor+(1|blockID)+(1|RILID), family=poisson, data=nolac2)
seddszegr<-glmer(Total.seed.weight~(1|Early.Growth..cm.)+Cor_Len_Var+(1|blockID)+(1|RILID), family=poisson, data=nolac2)
#none of this really makes sense

#let's try e.g. as offset
fset<-glmer(Total.seed.weight~RILcolor+Cor_Len_Var*RILcolor+Seed.Mass*RILcolor+Sugar_UND*RILcolor+(1|blockID)+(1|RILID)+offset(log(Early.Growth..cm.)), family=poisson, data=nolac2)
fsetnsz<-glmer(Total.seed.weight~RILcolor+Seed.Mass*RILcolor+Sugar_UND*RILcolor+(1|blockID)+(1|RILID)+offset(log(Early.Growth..cm.)), family=poisson, data=nolac2)
fsetcs<-glmer(Total.seed.weight~RILcolor+Sugar_UND*RILcolor+(1|blockID)+(1|RILID)+offset(log(Early.Growth..cm.)), family=poisson, data=nolac2)
fsets<-glmer(Total.seed.weight~Sugar_UND+(1|blockID)+(1|RILID)+offset(log(Early.Growth..cm.)), family=poisson, data=nolac2)
fsetc<-glmer(Total.seed.weight~RILcolor+(1|blockID)+(1|RILID)+offset(log(Early.Growth..cm.)), family=poisson, data=nolac2)



reseed2<-resid(seddwtx2)
plot(fitted(seddwtx2,reseed2))
qqnorm(reseed2)
qqline(reseed2)
plot(density(reseed2))

seedwtc2<-glmer(Total.seed.weight~Early.Growth..cm.+RILcolor+(1|blockID)+(1|RILID), family=poisson, data=nolac2)
seed2m<-glmer(Total.seed.weight~Early.Growth..cm.+Cor_Len_Var*Sugar_UND*RILcolor+(1|blockID)+(1|RILID), family=poisson, data=nolac2)
seed2megr<-glmer(Total.seed.weight~(1|Early.Growth..cm.)+Cor_Len_Var*Sugar_UND*RILcolor+(1|blockID)+(1|RILID), family=poisson, data=nolac2)

seedwtcegr2<-glmer(Total.seed.weight~(1|Early.Growth..cm.)+RILcolor+(1|blockID)+(1|RILID), family=poisson, data=nolac2)


lnest<-lmer(as.numeric(Total.seed.weight)~(1|blockID)+RILcolor+(1|RILID)+(1|RILID:RILcolor), data=nolac)
Anova(lmer(as.numeric(Total.seed.weight)~(1|blockID)+RILcolor+(1|RILID), data=nolac))

lnestsize<-lmer(as.numeric(Total.seed.weight)~(1|blockID)+RILcolor+Cor_Len_Var+Cor_Len_Var*RILcolor+(1|RILID)+(1|RILID:RILcolor) +(1|RILID:Cor_Len_Var), data=nolac)
lnestsug<-lmer(as.numeric(Total.seed.weight)~(1|blockID)+RILcolor+Sugar_UND+Sugar_UND*RILcolor+(1|RILID)+(1|RILID:RILcolor) +(1|RILID:Sugar_UND), data=nolac)
Anova(lmer(as.numeric(Total.seed.weight)~(1|blockID)+Sugar_UND+(1|RILID), data=nolac))

lnestsm<-lmer(as.numeric(Total.seed.weight)~(1|blockID)+RILcolor+Seed.Mass+Seed.Mass*RILcolor+(1|RILID)+(1|RILID:RILcolor) +(1|RILID:Seed.Mass), data=nolac)
Anova(lmer(as.numeric(Total.seed.weight)~(1|blockID)+Seed.Mass+(1|RILID), data=nolac))
lnestnec<-lmer(as.numeric(Total.seed.weight)~(1|blockID)+RILcolor+Nectar_UNC+Nectar_UNC*RILcolor+(1|RILID)+(1|RILID:RILcolor) +(1|RILID:Nectar_UNC), data=nolac)
lnestsl<-lmer(as.numeric(Total.seed.weight)~(1|blockID)+RILcolor+Seed.Length+Seed.Length*RILcolor+(1|RILID)+(1|RILID:RILcolor) +(1|RILID:Seed.Length), data=nolac)
Anova(lnest)
Anova(lnestsug)
Anova(lnestsm) #ignore, not good predictor
Anova(lnestnec) #ignore, not good predictor

Anova(lmer(as.numeric(Total.seed.weight)~(1|blockID)+RILcolor+Early.Growth..cm.+Early.Growth..cm.*RILcolor+(1|RILID)+(1|RILID:RILcolor), data=nolac))
Anova(lmer(as.numeric(Total.seed.weight)~(1|blockID)+RILcolor+Shoot..+Shoot..*RILcolor+(1|RILID)+(1|RILID:RILcolor), data=nolac))
Anova(lmer(as.numeric(Total.seed.weight)~(1|blockID)+Herbiv+(1|RILID), data=nolac)) #ignore herbiv
Anova(lnestsl)#ignore, not good predictor

nestgpegf<-lmer(as.numeric(Total.seed.weight)~(1|blockID)+Early.Growth..cm.+Shoot..+RILcolor+Sugar_UND+Cor_Len_Var+Sugar_UND*RILcolor+Cor_Len_Var*RILcolor+(1|RILID)+(RILcolor|RILID)+(Cor_Len_Var|RILID)+(Sugar_UND|RILID), data=nolac)
nestgpegr<-lmer(as.numeric(Total.seed.weight)~(1|blockID)+(1|Early.Growth..cm.)+(1|Shoot..)+RILcolor+Sugar_UND+Cor_Len_Var+Sugar_UND*RILcolor+Cor_Len_Var*RILcolor+(1|RILID)+(RILcolor|RILID)+(Cor_Len_Var|RILID)+(Sugar_UND|RILID), data=nolac)
nestgpegrnc<-lmer(as.numeric(Total.seed.weight)~(1|blockID)+(1|Early.Growth..cm.)+(1|Shoot..)+Sugar_UND+Cor_Len_Var+Sugar_UND*RILcolor+Cor_Len_Var*RILcolor+(1|RILID)+(RILcolor|RILID)+(Cor_Len_Var|RILID)+(Sugar_UND|RILID), data=nolac)
nestmegrnc<-lmer(as.numeric(Total.seed.weight)~(1|blockID)+(1|Early.Growth..cm.)+(1|Shoot..)+Cor_Len_Var*Sugar_UND*RILcolor+(1|RILID), data=nolac)



lnestall0<-lmer(as.numeric(Total.seed.weight)~(1|blockID)+Early.Growth..cm.+Shoot..+Herbiv +Cor_Len_Var+Nectar_UNC+Sugar_UND+RILcolor+Seed.Mass+Seed.Length+Seed.Length*RILcolor+Seed.Mass*RILcolor+Early.Growth..cm.*RILcolor+Sugar_UND*RILcolor+Nectar_UNC*RILcolor+Shoot..*RILcolor+Cor_Len_Var*RILcolor+Herbiv*RILcolor+(1|RILID)+(1|RILID:RILcolor)+(1|RILID:Cor_Len_Var)+(1|RILID:Sugar_UND) +(1|RILID:Seed.Mass)+(1|RILID:Herbiv), data=nolac)
lnestall6<-lmer(as.numeric(Total.seed.weight)~(1|blockID)+Early.Growth..cm.+Herbiv +Cor_Len_Var+RILcolor+Seed.Mass+Seed.Mass*RILcolor+Early.Growth..cm.*RILcolor+Sugar_UND*RILcolor+Cor_Len_Var*RILcolor+Herbiv*RILcolor+(1|RILID)+(1|RILID:RILcolor)+(1|RILID:Cor_Len_Var)+(1|RILID:Sugar_UND) +(1|RILID:Seed.Mass)+(1|RILID:Herbiv), data=nolac)
lnestall5egr<-lmer(as.numeric(Total.seed.weight)~(1|blockID)+(1|Early.Growth..cm.)+Herbiv +Cor_Len_Var+RILcolor+Seed.Mass+Seed.Mass*RILcolor+Sugar_UND*RILcolor+Cor_Len_Var*RILcolor+(1|RILID)+(1|RILID:RILcolor)+(1|RILID:Cor_Len_Var)+(1|RILID:Sugar_UND) +(1|RILID:Seed.Mass)+(1|RILID:Herbiv), data=nolac)
simpre<-lmer(as.numeric(Total.seed.weight)~(1|blockID)+(1|Early.Growth..cm.)+Herbiv +Cor_Len_Var+RILcolor+Seed.Mass+Seed.Mass*RILcolor+Sugar_UND*RILcolor+Cor_Len_Var*RILcolor+(1|RILID), data=nolac)


nest<- aov(nolac$Herbiv~nolac$blockID+ nolac$RILcolor + Error(1/factor(nolac$RILID)))

nest <- aov(nolac$Herbiv~nolac$RILcolor/ factor(nolac$RILID))#w/o block
summary(nest)#suggests that while color is a significant predictor of herbivory (white=more herbiv),line is a much stronger predictor
nestb<-aov(nolac$Herbiv~factor(nolac$blockID)+nolac$RILcolor/factor(nolac$RILID))#w/o block
summary(nestb) #color is a predictor, but block is a much stronger predictor
plot(nolac$Herbiv, nolac$Total.seed.weight, col=as.numeric(as.factor(nolac$RILcolor)))

#for seed set (via weight)
summary(lm(as.numeric(nolac$Total.seed.weight)~nolac$RILcolor))
summary(aov(as.numeric(nolac$Total.seed.weight)~nolac$RILcolor))
nestw<- aov(as.numeric(nolac$Total.seed.weight)~nolac$RILcolor/factor(nolac$RILID)) #RIL ID is highly determinant of seed set, color merely suggestive
summary(nestw)
#w/block
nestwb<- aov(as.numeric(nolac$Total.seed.weight)~factor(nolac$blockID)+nolac$RILcolor/factor(nolac$RILID), na.rm=T) #RIL ID is highly determinant of seed set, color merely suggestive
summary(nestwb) #block and color are both suggestive, again biggest effect is by line
boxplot(as.numeric(nolac$Total.seed.weight)~nolac$RILcolor, xlab="Color", ylab="Seed Weight (g)")


#let's look at the other purport ns factors
nesteg <- aov(nolac$Early.Growth..cm.~nolac$RILcolor/ factor(nolac$RILID))#w/o block
summary(nesteg)#color is not sig predictor
nestegb<-aov(nolac$Early.Growth..cm.~factor(nolac$blockID)+nolac$RILcolor/factor(nolac$RILID))#w/o block
summary(nestegb)#but not when factoring in block

nestsn <- aov(nolac$Shoot..~nolac$RILcolor/ factor(nolac$RILID))#w/o block
summary(nestsn)#color is now sig predictor
nestsnb<-aov(nolac$Shoot..~factor(nolac$blockID)+nolac$RILcolor/factor(nolac$RILID))#w/o block
summary(nestsnb)

head(nolac)


##################################################################
#calculating survival success based off of flower color
##################################################################
#estimting survival
mdata1<-masterdata <-masterdata[1:1199,]
n <-nrow(mdata1)
dvec <- rep(NA,n)

for(i in 1:n){
  if (masterdata$Dead.[i]==0 | is.na(masterdata$Dead.[i])){
    mdata1$Dead.[i] <-0
    dvec[i] <-0
  } else {
    mdata1$Dead.[i] <-1
    dvec[i] <-1
  }
}

sum(dvec)
mdata1<-cbind(mdata1, dvec)
blockID <- gsub('\\b(\\pL)|.','\\U\\1',mdata1$X,perl = TRUE)#adapted from stackoverflow
mdata1<-cbind(mdata1,blockID)
nldata<-subset(mdata1, mdata1$RILcolor!="L")
npdata<-subset(mdata1, mdata1$RILcolor!="P")
nwdata<-subset(mdata1, mdata1$RILcolor!="W")
pdatax<-subset(mdata1, mdata1$RILcolor=="P")
wdatax<-subset(mdata1, mdata1$RILcolor=="W")
pwdatax<-rbind(pdatax,wdatax)
pw1<-pwdatax[,2:3]
pwdataxx<-cbind(pw1,pwdatax$Dead.)
colnames(pwdataxx)=c("RILID", "RILcolor", "Dead.")

nxn<-nlevels(as.factor(nldata$RILID))

dlmx<-matrix(data=NA, nrow=nxn, ncol=4)
group_order <- names(sort(tapply(as.numeric(as.factor(pwdataxx$RILcolor)), pwdataxx$RILID, median)))
pwdataxx$RILID<- factor(pwdataxx$RILID, levels=group_order)

coline<-rep(NA,nxn)
for (i in 1:nxn){
  tempmtx<-subset(pwdataxx, as.numeric(as.factor(pwdataxx$RILID))==i)
  cID<-tempmtx$RILID[1]
  ccol<-tempmtx$RILcolor[1]
  dmean<-mean(as.numeric(tempmtx$Dead.))
  dstd<-sqrt(var(tempmtx$Dead.))
  dlmx[i,]<- c(cID,ccol, dmean,dstd)
  if (tempmtx$RILcolor[1]=="P"){
    coline[i]="magenta"
    print("pink")
  } else if (tempmtx$RILcolor[1]=="W") {
    coline[i]="black"
    print("white")
  } else{
    print("error")
  }
}

colnames(dlmx)=c("RILID", "Color", "ave", "sd")
dlmx[,3]<-as.numeric(dlmx[,3])
dlmx[,4]<-as.numeric(dlmx[,4])
pdlmx<-subset(dlmx, dlmx[,2]=="P")
wdlmx<-subset(dlmx, dlmx[,2]=="W")
mean(1-as.numeric(pdlmx[,3])) #average pink survival/line
sdpmx<-sqrt(var(1-as.numeric(pdlmx[,3])))#sd
sdpmx/sqrt(nrow(pdlmx))
mean(1-as.numeric(wdlmx[,3]))#average white survival/line
sdwmx<-sqrt(var(1-as.numeric(wdlmx[,3]))) #sd
sdwmx/sqrt(nrow(pdlmx))
pontius<-1-as.numeric(dlmx[,3])

bx<-boxplot((rep(1, 798)-as.numeric(pwdataxx$Dead.))~pwdataxx$RILID, xlab="Line", ylab="Average Survival Rate", border=coline,xaxt="n")
points(pontius, col=coline, cex=1.5, pch=18)
tick <- seq_along(bx$names)
axis(1, at = tick, labels = FALSE)
text(tick, par("usr")[3] - 0.075, bx$names, srt = 90, xpd = TRUE)

br<-barplot(dlmx[,3]~x)



ptabs<-table(1-as.numeric(pdatax$Dead.))
wtabs<-table(1-as.numeric(wdatax$Dead.))
shalla<-rbind(ptabs,wtabs)
row.names(shalla)<-c("P","W")
colnames(shalla)<-c("Dead", "Alive")
barplot((shalla), xlab="Survival", ylab="Number of Plants", col=c("magenta", "white"), beside=T)

#better graph would have counts on y-axis of survival?
mdataw<-cbind(mdata1$RILID,mdata1$RILcolor,mdata1$Dead.)
mdead<-subset(mdataw,as.numeric(mdata1$Dead.)>0)
pdead<-nrow(subset(mdead,mdead[,2]=="P"))
wdead<-nrow(subset(mdead,mdead[,2]=="W"))
ldead<-nrow(mdead)-pdead-wdead
mlive<-subset(mdataw,as.numeric(mdata1$Dead.)==0)
plive<-nrow(subset(mlive,mlive[,2]=="P"))
wlive<-nrow(subset(mlive,mlive[,2]=="W"))
llive<-nrow(mlive)-plive-wlive #these measurments form a table for data vis

summary(lm(mdata1$dvec~mdata1$RILcolor))
d2all<-glmer(dvec~RILcolor+Cor_Len_Var*RILcolor+Seed.Mass*RILcolor+Sugar_UND*RILcolor+(1|blockID)+(1|RILID), family=binomial, data=mdata1)
d2c<- glmer(dvec~RILcolor+(1|blockID)+(1|RILID), family=binomial, data=mdata1)
summary(d2c)
Anova(d2c, type=2) #this is reported p-values/chi-squared
d2csg<-glmer(dvec~Sugar_UND*RILcolor+(1|blockID)+(1|RILID), family=binomial, data=mdata1) #error here
d2cfz<-glmer(dvec~Cor_Len_Var*RILcolor+(1|blockID)+(1|RILID), family=binomial, data=mdata1)

#glm with RIL as random effect clearly shows RIL accounts for most of the variance.
#no sig diff with urvivval across colors/sizes/ etc.

summary(aov(mdata1$dvec~mdata1$RILcolor))
nestdt <- aov(mdata1$dvec~mdata1$RILcolor/factor(mdata1$RILID))#w/o block
nestdtwb <- aov(mdata1$dvec~mdata1$blockID+mdata1$RILcolor/factor(mdata1$RILID))#w/o block
summary(nestdtwb) #block and lineare best indicators of survival



#w/o lac
summary(aov(nldata$dvec~nldata$RILcolor))

nestdt2 <- aov(nldata$dvec~nldata$RILcolor/factor(nldata$RILID))#w/o block
nestdtwb2 <- aov(nldata$dvec~nldata$blockID+nldata$RILcolor+factor(nldata$RILID))#w/o block
summary(nestdtwb2) #block and lineare best indicators of survival
d2nl<- glmer(dvec~RILcolor+(1|blockID)+(1|RILID), family=binomial, data=nldata)
summary(d2nl)
#w/o pink
summary(aov(npdata$dvec~npdata$RILcolor))
nestdt3 <- aov(npdata$dvec~npdata$RILcolor/factor(npdata$RILID))#w/o block
nestdtwb3 <- aov(npdata$dvec~npdata$blockID+npdata$RILcolor/factor(npdata$RILID))#w/o block
summary(nestdtwb3)
d2np<- glmer(dvec~RILcolor+(1|blockID)+(1|RILID), family=binomial, data=npdata)
summary(d2np)
d2nw<- glmer(dvec~RILcolor+(1|blockID)+(1|RILID), family=binomial, data=nwdata)
summary(d2nw)

#great examples of where linear model is wrong approach
#when considering as a binomial glm, nothing predicts death (RILs are responsible for nearly all variance)

#early growth, same approach
boxplot(as.numeric(mdata1$Early.Growth..cm.)*2.54~mdata1$RILcolor, xlab="Line Color", ylab="Shoot Length at 2 weeks (cm)")
boxplot(as.numeric(nldata$Early.Growth..cm.)*2.54~nldata$RILcolor, xlab="Line Color", ylab="Shoot Length at 2 weeks (cm)")


#all3
eg2<- lmer(Early.Growth..cm.~RILcolor+(1|blockID)+(1|RILID), data=mdata1)
summary(eg2)
#w/o lac
eg2nl<- lmer(Early.Growth..cm.~RILcolor+(1|blockID)+(1|RILID), data=nldata)
summary(eg2nl)
#w/o pink
eg2np<- lmer(Early.Growth..cm.~RILcolor+(1|blockID)+(1|RILID), data=npdata)
summary(eg2np)
eg2nw<- lmer(Early.Growth..cm.~RILcolor+(1|blockID)+(1|RILID), data=nwdata)
summary(eg2nw)

#line based analysis, to incorporate flower size and s.c.
linenumbs<-nlevels(as.factor(nldata$RILID))
linematrixs<-matrix(data=NA, nrow=linenumbs, ncol=7)
for (i in 1:linenumbs){
  templine<-subset(nldata, as.numeric(as.factor(nldata$RILID))==i)
  fztemp<-mean(templine$Cor_Len_Var, na.rm=T)
  sctemp<-mean(templine$Sugar_UND, na.rm=T)
  sline<-1-mean(as.numeric(templine$Dead.), na.rm=T)
  smass<-mean(templine$Seed.Mass, na.rm=T)
  coltemp<-templine$RILcolor[1]
  linehold<-templine$RILID[1]
  nlines<-nrow(templine)
  linedatas<-c(linehold, coltemp,as.numeric(fztemp),as.numeric(sctemp),as.numeric(smass), as.numeric(sline), as.numeric(nlines))
  linematrixs[i,]<-linedatas
}
linematrixs<-as.data.frame(linematrixs)
labelens<- c("Line", "Color", "Flower_Size", "Sugar_Conc", "Seed_Mass","Survival","LineWeight")
colnames(linematrixs)<-labelens
linematrixs$Flower_Size<-as.numeric(linematrixs$Flower_Size)
linematrixs$Sugar_Conc<-as.numeric(linematrixs$Sugar_Conc)
linematrixs$Seed_Mass<-as.numeric(linematrixs$Seed_Mass)
linematrixs$Survival<-as.numeric(linematrixs$Survival)
linematrixs$LineWeight<-as.numeric(linematrixs$LineWeight)

boxplot(linematrixs$Survival~linematrixs$Color)
plot(linematrixs$Survival~linematrixs$Sugar_Conc)
abline(lm(linematrixs$Survival~linematrixs$Sugar_Conc))
plot(linematrixs$Survival~linematrixs$Flower_Size)
abline(lm(linematrixs$Survival~linematrixs$Flower_Size))
plot(linematrixs$Survival~linematrixs$Seed_Mass)
abline(lm(linematrixs$Survival~linematrixs$Seed_Mass))

pines<-subset(linematrixs, linematrix$Color=="P")
wines<-subset(linematrixs, linematrix$Color=="W")
plot(linematrixs$Survival~linematrixs$Flower_Size, xlab="Corolla Length (mm)", ylab="Survival(g)", col=as.numeric(as.factor(linematrix$Color))%%2*5+1, pch=as.numeric(as.factor(linematrix$Color))*2+17%%20)
abline(lm(pines$Survival~pines$Flower_Size), col="magenta", lwd=2)
abline(lm(wines$Survival~wines$Flower_Size), col="black", lwd=2)
summary(lm(linematrixs$Survival~linematrixs$Flower_Size))
summary(glm(linematrixs$Survival~linematrixs$Flower_Size, family=binomial, weights=linematrixs$LineWeight))
summary(lm(linematrixs$Survival~linematrixs$Flower_Size*linematrixs$Color))
summary(glm(linematrixs$Survival~linematrixs$Flower_Size*linematrixs$Color, family=binomial, weights=linematrixs$LineWeight))
Anova(glm(linematrixs$Survival~linematrixs$Flower_Size*linematrixs$Color, family=binomial, weights=linematrixs$LineWeight), type=3, contrasts=list(topic=contr.sum, sys=contr.sum))

plot(linematrixs$Survival~linematrixs$Flower_Size, xlab="Corolla Length (mm)", ylab="Survival Rate", col=as.numeric(as.factor(linematrix$Color))%%2*5+1, pch=as.numeric(as.factor(linematrix$Color))*2+17%%20)
pfit<-(glm(Survival~Flower_Size, data=pines, family=binomial, weights=LineWeight)) #
wfit<-(glm(Survival~Flower_Size, data=wines, family=binomial, weights=LineWeight)) #
pblankz <- data.frame(Flower_Size=seq(0, max(pines$Flower_Size)+5,len=100))
wblankz <- data.frame(Flower_Size=seq(0, max(pines$Flower_Size)+5,len=100))
pblankz$Survival<-predict(pfit, newdata=pblankz,type="response")
wblankz$Survival<-predict(wfit, newdata=wblankz,type="response")
lines(Survival ~ Flower_Size, pblankz, col="magenta", lwd=2)
lines(Survival ~ Flower_Size, wblankz, col="black", lwd=2)

plot(linematrixs$Survival~linematrixs$Sugar_Conc, xlab="Sugar Conc. (mg/ml)", ylab="Survival(g)", col=as.numeric(as.factor(linematrix$Color))%%2*5+1, pch=as.numeric(as.factor(linematrix$Color))*2+17%%20)
abline(lm(pines$Survival~pines$Sugar_Conc), col="magenta", lwd=2)
abline(lm(wines$Survival~wines$Sugar_Conc), col="black", lwd=2)
summary(lm(linematrixs$Survival~linematrixs$Sugar_Conc))
summary(glm(linematrixs$Survival~linematrixs$Sugar_Conc, family=binomial, weights=linematrixs$LineWeight))
summary(lm(linematrixs$Survival~linematrixs$Sugar_Conc*linematrixs$Color))
summary(glm(linematrixs$Survival~linematrixs$Sugar_Conc*linematrixs$Color, family=binomial, weights=linematrixs$LineWeight))
Anova(glm(linematrixs$Survival~linematrixs$Sugar_Conc*linematrixs$Color, family=binomial, weights=linematrixs$LineWeight), type=3, contrasts=list(topic=contr.sum, sys=contr.sum))

plot(linematrixs$Survival~linematrixs$Sugar_Conc, xlab="Sugar Conc. (mg/ml)", ylab="Survival Rate", col=as.numeric(as.factor(linematrix$Color))%%2*5+1, pch=as.numeric(as.factor(linematrix$Color))*2+17%%20)
pfit<-(glm(Survival~Sugar_Conc, data=pines, family=binomial, weights=LineWeight)) #
wfit<-(glm(Survival~Sugar_Conc, data=wines, family=binomial, weights=LineWeight)) #

pblankz <- data.frame(Sugar_Conc=seq(0, max(pines$Sugar_Conc)+5,len=100))
wblankz <- data.frame(Sugar_Conc=seq(0, max(pines$Sugar_Conc)+5,len=100))
pblankz$Survival<-predict(pfit, newdata=pblankz,type="response")
wblankz$Survival<-predict(wfit, newdata=wblankz,type="response")
lines(Survival ~ Sugar_Conc, pblankz, col="magenta", lwd=2)
lines(Survival ~ Sugar_Conc, wblankz, col="black", lwd=2)

plot(linematrixs$Survival~linematrixs$Seed_Mass, xlab="Mass/Seed(g)", ylab="Survival", col=as.numeric(as.factor(linematrix$Color))%%2*5+1, pch=as.numeric(as.factor(linematrix$Color))*2+17%%20)
abline(lm(pines$Survival~pines$Seed_Mass), col="magenta", lwd=2)
abline(lm(wines$Survival~wines$Seed_Mass), col="black", lwd=2)
summary(lm(linematrixs$Survival~linematrixs$Seed_Mass))
summary(glm(linematrixs$Survival~linematrixs$Seed_Mass, family=binomial, , weights=linematrixs$LineWeight))
summary(lm(linematrixs$Survival~linematrixs$Seed_Mass*linematrixs$Color))
summary(glm(linematrixs$Survival~linematrixs$Seed_Mass*linematrixs$Color, family=binomial, weights=linematrixs$LineWeight))
Anova(glm(linematrixs$Survival~linematrixs$Seed_Mass*linematrixs$Color, family=binomial, weights=linematrixs$LineWeight), type=3, contrasts=list(topic=contr.sum, sys=contr.sum))

plot(linematrixs$Survival~linematrixs$Seed_Mass, xlab="Mass/Seed(g)", ylab="Survival Rate", col=as.numeric(as.factor(linematrix$Color))%%2*5+1, pch=as.numeric(as.factor(linematrix$Color))*2+17%%20)
pfit<-(glm(Survival~Seed_Mass, data=pines, family=binomial, weights=LineWeight)) #
wfit<-(glm(Survival~Seed_Mass, data=wines, family=binomial, weights=LineWeight)) #
pblankz <- data.frame(Seed_Mass=seq(0, max(pines$Seed_Mass)+5,len=100))
wblankz <- data.frame(Seed_Mass=seq(0, max(pines$Seed_Mass)+5,len=100))
pblankz$Survival<-predict(pfit, newdata=pblankz,type="response")
wblankz$Survival<-predict(wfit, newdata=wblankz,type="response")
lines(Survival ~ Seed_Mass, pblankz, col="magenta", lwd=2)
lines(Survival ~ Seed_Mass, wblankz, col="black", lwd=2)
#binomials complicate things

########incorporating germination data###############
germcsv<-read.csv("C:/Users/jonat/OneDrive/Documents/Grad School Research Files/germdata.csv")
germd<-germcsv[1:800,1:26]

group_order2 <- names(sort(tapply(as.numeric(as.factor(germd$Color)), germd$Line, median)))
germd$Line<- factor(germd$Line, levels=group_order2)

coline<-rep(NA,nlevels(as.factor(germd$Line)))
for (i in 1:nlevels(as.factor(germd$Line))){
  tempmtx<-subset(germd, as.numeric(as.factor(germd$Line))==i)
  if (tempmtx$Color[1]=="P"){
    coline[i]="magenta"
    print("pink")
  } else if (tempmtx$Color[1]=="W") {
    coline[i]="white"
    print("white")
  } else{
    print("error")
  }
}

boxplot(as.numeric(germd$Total)~germd$Color, xlab="Line Color", ylab="Total Germinant Number")
boxplot(as.numeric(germd$Germination.Fraction)~germd$Color, xlab="Line", ylab="Germination Fraction")

bt<-boxplot(as.numeric(germd$Germination.Fraction)~germd$Line,xlab="Line", ylab="Germination Fraction",col=coline, xaxt="n")
tick <- seq_along(bt$names)
axis(1, at = tick, labels = FALSE)
text(tick, par("usr")[3] - 0.075, bt$names, srt = 90, xpd = TRUE)


germbb<- glmer(Total~Color+(1|Line), family=poisson, data=germd)
#seperating out row and column (space)
rc<-matrix(unlist(strsplit(germd$DATES,"R", "C")),ncol=2,byrow =T)[,2]
rc2<-matrix(as.numeric(unlist(strsplit(rc,"C"))),ncol=2,byrow =T)
colnames(rc2)<-c("row", "column")
germxy<-cbind(germd[,2:ncol(germd)],rc2)
germxy[is.na(germxy)] <- 0
germbxy<- glmer(Total~Color+(1|Line)+(1|row)+(1|column), family=poisson, data=germxy)
Anova(germbxy, type=2)
germbxy1w<- glmer((X17.Jun+X24.Jun+X30.Jun+X8.Jul)~Color+(1|Line)+(1|row)+(1|column), family=poisson, data=germxy)
boxplot((X17.Jun+X24.Jun+X30.Jun+X8.Jul)~Color,data=germxy)
#no sig dif by color

#over time?
install.packages("tidyr")
library(tidyr)
colnames(germxy)<-c("Color", "Line", 7,14,20,28,33,38,41,47,52,59,63,70,77,89,97,103,115,119,123,130,160, "Total", "Germination.Fraction","row","column") #days after planting (June 10)
gxylong<-pivot_longer(germxy,col=3:23,names_to ="days" , values_to="germinants",values_drop_na = FALSE)
gxylong$days<-as.numeric(gxylong$days)
gxylong$Germination.Fraction<-gxylong$germinants/10
pgxy<-subset(gxylong, gxylong$Color=="P")
wgxy<-subset(gxylong, gxylong$Color=="W")
plot(jitter(germinants,0.2)~days, data=gxylong, ylab="germinants", col=as.numeric(as.factor(Color))%%2*5+1, pch=as.numeric(as.factor(Color))*2+17%%20)
summary(lmer(germinants~Color*days+(1|row)+(1|column)+(1|Line:Color),data=gxylong))
Anova(lmer(germinants~Color*days+(1|row)+(1|column)+(1|Line:Color),data=gxylong), type=3, contrasts=list(topic=contr.sum, sys=contr.sum))

abline(lm(germinants~days, data=pgxy), col="magenta", lwd=2)
abline(lm(germinants~days, data=wgxy), col="black", lwd=2)

summary(glmer(Germination.Fraction~Color*days+(1|row)+(1|column)+(1|Line:Color),data=gxylong, family=binomial, weights=rep(10, nrow(gxylong))))
Anova(glmer(Germination.Fraction~Color*days+(1|row)+(1|column)+(1|Line:Color),data=gxylong, family=binomial, weights=rep(10, nrow(gxylong))),  type=3, contrasts=list(topic=contr.sum, sys=contr.sum))
plot(jitter(Germination.Fraction,0.2)~days, data=gxylong, ylab="germinantion fraction", col=as.numeric(as.factor(Color))%%2*5+1, pch=as.numeric(as.factor(Color))*2+17%%20)
pfit<-(glm(as.numeric(Germination.Fraction)~as.numeric(days), data=pgxy, family=binomial, weights=rep(10, nrow(pgxy)))) #
wfit<-(glm(as.numeric(Germination.Fraction)~as.numeric(days), data=wgxy, family=binomial, weights=rep(10, nrow(wgxy)))) #
pblankz <- data.frame(days=seq(0, max(as.numeric(plong$days))+5,len=100))
wblankz <- data.frame(days=seq(0, max(as.numeric(plong$days))+5,len=100))
pblankz$Germination.Fraction<-predict(pfit, newdata=pblankz,type="response")
wblankz$Germination.Fraction<-predict(wfit, newdata=wblankz,type="response")
lines(Germination.Fraction ~ days, pblankz, col="magenta", lwd=2)
lines(Germination.Fraction ~ days, wblankz, col="black", lwd=2)
oldprop<-rep(NA,6)
for (i in 1:nlevels(as.factor(pgxy$days))){
  vprop<-table(subset(pgxy,as.numeric(as.factor(pgxy$days))==i)$germinants)
  vhalla<-rbind(oldprop,vprop)
  oldprop<-vhalla
}
vhalla<-vhalla[2:22,]
rownames(vhalla)<-levels(as.factor(pgxy$days))
colnames(vhalla)<- c(0,0.1,0.2,0.3,0.4,0.5)
vhalla[,6]<-rep(0,21)
vhalla[c(11,16,18),2:5]<-0
vhalla[c(1,6:7,9:21),3:5]<-0
vhalla[c(3,5,8),4:5]<-0
barplot(prop.table(t(vhalla[])), xlab="Date", ylab="Frequency", col=c("black","red1", "orange1", "green1", "blue1","purple1"), beside=T)
barplot(prop.table(t(vhalla[,2:6])), xlab="Date", ylab="Frequency", col=c("red1", "orange1", "green1", "blue1","purple1"), beside=T)
tdist<-(t(vhalla[,2:6]))[1,]+(t(vhalla[,2:6]))[2,]*2+(t(vhalla[,2:6]))[3,]*3+(t(vhalla[,2:6]))[4,]*4+(t(vhalla[,2:6]))[5,]*5
tdist<-rbind(as.numeric(levels(as.factor(pgxy$days))),tdist)
weighted.mean(tdist[1,],tdist[2,])
df.freq= as.vector(rep(tdist[1,], tdist[2,]))
hist(df.freq, xlab="Date",ylab="Number of Germinants", col=c("magenta"),breaks=seq(7,160,l=11),xlim=c(0,160))
abline(v=weighted.mean(tdist[1,],tdist[2,]),lty=3,lwd=3)

oldprop<-rep(NA,6)
for (i in 1:nlevels(as.factor(wgxy$days))){
  vprop<-table(subset(wgxy,as.numeric(as.factor(wgxy$days))==i)$Germination.Fraction)
  vhalla2<-rbind(oldprop,vprop)
  oldprop<-vhalla2
}
vhalla2<-vhalla2[2:22,]
rownames(vhalla2)<-levels(as.factor(wgxy$days))
colnames(vhalla2)<- c(0,0.1,0.2,0.3,0.4,0.5)
vhalla2[,6]<-rep(0,21)
vhalla2[2,4:6]<-c(0,1,1)
vhalla2[c(1,9,10,17,21),2:5]<-0
vhalla2[c(3,6:20),3:5]<-0
vhalla2[5,4:5]<-0
barplot(prop.table(t(vhalla2[])), xlab="Date", ylab="Frequency", col=c("black","red1", "orange1", "green1", "blue1","purple1"), beside=T)
barplot(prop.table(t(vhalla2[,2:6])), xlab="Date", ylab="Frequency", col=c("red1", "orange1", "green1", "blue1","purple1"), beside=T)
tdist2<-(t(vhalla2[,2:6]))[1,]+(t(vhalla2[,2:6]))[2,]*2+(t(vhalla2[,2:6]))[3,]*3+(t(vhalla2[,2:6]))[4,]*4+(t(vhalla2[,2:6]))[5,]*5
tdist2<-rbind(as.numeric(levels(as.factor(wgxy$days))),tdist2)
weighted.mean(tdist2[1,],tdist2[2,])
df.freq2= as.vector(rep(tdist2[1,], tdist2[2,]))
hist(df.freq2, xlab="Date",ylab="Number of Germinants", col=c("white"), breaks=seq(7,160,l=11),xlim=c(0,160))
abline(v=weighted.mean(tdist2[1,],tdist2[2,]),lty=3,lwd=3)

#abline(lm(Germination.Fraction~days, data=pgxy), col="magenta", lwd=2)
#abline(lm(Germination.Fraction~days, data=wgxy), col="black", lwd=2)
#sig w/time*color but estimates are miniscule for joint effects (estimate for pink favored by rate ~0.00699 germ fraction/day)

#now to make line estimates of germiantion rate
glineco<-nlevels(as.factor(as.character(germxy$Line)))
linestoreg<-matrix(data=NA, nrow=glineco,ncol=ncol(germxy)+1)
#calculates averages across lines, if needed
for (i in 1:glineco){
  wogerm<-subset(germxy,as.numeric(as.factor(germd$Line))==i)
  colhold<-wogerm[1,1]
  ghld<-as.data.frame(wogerm[,2:ncol(wogerm)])
  ghld$Line<-as.numeric(as.character(ghld$Line))
  gerlineave<-colMeans(ghld,na.rm=T)
  gweights<-nrow(wogerm)
  linestoreg[i,]<-c(gerlineave,colhold,as.numeric(gweights))
}
linestoreg<-as.data.frame(linestoreg)
colnames(linestoreg)<-c(colnames(germxy[2:27]),"Color","gweight")
#for color gweights
linestoreg$gweight<-as.numeric(linestoreg$gweight)
linestoreg$Germination.Fraction<-as.numeric(linestoreg$Germination.Fraction)
summary(lm(Germination.Fraction~Color,data=linestoreg)) #evenly balanced
summary(glm(Germination.Fraction~as.factor(Color),data=linestoreg, family=quasibinomial, weights=gweight)) #needs to be as.factor for model run

gerline<- lmer(as.numeric(Germination.Fraction)~Color+(1|row)+(1|column), data=linestoreg)
summary(gerline)
Anova(gerline, type=2)
gerlineg<- glmer(as.numeric(Germination.Fraction)~Color+(1|row)+(1|column), data=linestoreg,family=binomial, weights=gweight)
summary(gerlineg)

group_order <- names(sort(tapply(germd$Total, germd$Line, mean)))
germd$Line <- factor(germd$Line, levels=group_order)

boxplot(as.numeric(germd$Germination.Fraction)~germd$Line,xlab="Line",col=as.numeric(as.factor(linestoreg$Color))+5, ylab="Germination Fraction")

boxplot(as.numeric(linestoreg$Germination.Fraction)~linestoreg$Color, xlab="Line Color", ylab="Germination Fraction")

#with using line estimates, germination rate is stil same regardless of color

######should probably incorporate line estimates of sugar conc, flower size, and seed mass(?) here?

#with time
colnames(linestoreg)<-c("Line", 7,14,20,28,33,38,41,47,52,59,63,70,77,89,97,103,115,119,123,130,160, "Total", "Germination.Fraction","row","column","Color","gweight") #days after planting (June 10)
llong<-pivot_longer(linestoreg,cols=2:22, names_to ="days" , values_to="germinants",values_drop_na = FALSE)
plong<-subset(llong, llong$Color=="P")
wlong<-subset(llong, llong$Color=="W")
plot(germinants~days, data=llong)
plot(jitter(as.numeric(germinants),10)~days, data=llong, ylab="germinants", col=as.numeric(as.factor(Color))%%2*5+1, pch=as.numeric(as.factor(Color))*2+17%%20)
abline(lm(as.numeric(germinants)~as.numeric(days), data=plong), col="magenta", lwd=2)
abline(lm(as.numeric(germinants)~as.numeric(days), data=llong), col="black", lwd=2)

plot(jitter(as.numeric(germinants),10)~days, data=llong, ylab="germinants", col=as.numeric(as.factor(Color))%%2*5+1, pch=as.numeric(as.factor(Color))*2+17%%20)
pfit<-(glm(as.numeric(germinants)~as.numeric(days), data=plong, family=binomial, weights=gweight)) #
wfit<-(glm(as.numeric(germinants)~as.numeric(days), data=wlong, family=binomial, weights=gweight)) #
pblankz <- data.frame(days=seq(0, max(as.numeric(plong$days))+5,len=100))
wblankz <- data.frame(days=seq(0, max(as.numeric(plong$days))+5,len=100))
pblankz$germinants<-predict(pfit, newdata=pblankz,type="response")
wblankz$germinants<-predict(wfit, newdata=wblankz,type="response")
lines(germinants ~ days, pblankz, col="magenta", lwd=2)
lines(germinants ~ days, wblankz, col="black", lwd=2)

summary(lm(as.numeric(germinants)~as.numeric(days), data=llong))
Anova(lm(as.numeric(germinants)~as.numeric(days), data=llong))
Anova(lm(as.numeric(germinants)~as.numeric(days)*Color, data=llong))
Anova(glm(as.numeric(germinants)~as.numeric(days)*Color, data=llong),family=binomial, weights=gweight)
Anova(glmer(as.numeric(germinants)~as.numeric(days)*Color), data=llong ,family=binomial, weights=gweight)


#make first flower data

#split character along dash
#special condition for before dates "<9/30; <9/22"
#check for NAs, then <sign
#if less than sign, take last 2 characcters, convert to number, and add 242
#if August = 7*30+1
#if September=8*30+2
#if October=9*30+2
#add numbers together for Julian date



#running a linear model

library(lubridate)


#quantifying flower number
#converting dates to Julian dates
masterdate <- mdy(mdata1$Flower.Date)
masterdate1 <- mdy(mdata1$Flower.Date.1)
masterdate2 <- mdy(mdata1$Flower.Date.2)
masterdate3 <- mdy(mdata1$Flower.Date.3)

#testconversionplots
plot (x = masterdata$Early.Growth..cm, y = as.numeric(masterdate1))
plot (x=masterdate1, y=masterdata$Flower...1)
points(x=masterdate2, y=masterdata$Flower...2)
points(x=masterdate3, y=masterdata$Flower...3)

mdate1<-as.numeric(masterdate1)
mdate12<-as.numeric(masterdate1)^2
flowercurve <- lm(masterdata$Flower...1~ mdate1 + mdate12, data=masterdata)

#now we have to run this through a for loop for each RIL line, while holding each RIL-total flower numer
flowernumbhold <- rep (NA, nlevels(masterdata$RILID))
#inits
#plot (x=0, xlim=c(18850,18920), ylim=c(0,120))

#we need to reformat the data to iterate out this process
tempid <- masterdata[,1:3]
tempfl <- masterdata[,11:16]

tempfl$Flower.Date.1 <- masterdate1
tempfl$Flower.Date.2 <- masterdate2
tempfl$Flower.Date.3 <- masterdate3
tempa <- cbind(tempid, tempfl)
temp1 <-tempa[,1:5]
temp2a <-tempa[,6:7]
temp2 <-cbind(tempid,temp2a)
temp3a <-tempa[,8:9]
temp3 <-cbind(tempid,temp3a)
tempx <- rep(NA, ncol(temp1))

names(temp1) <-names(temp2)<- names(temp3)<- c('Location', 'RILID', 'RILcolor', 'FlowerNumber', 'FlowerDate')

for (i in 1:n){
  tempr <- temp1[i,]
  if (is.na(temp3[i,4])){
    if (is.na(temp2[i,4])) {
    } else {
      tempr <-rbind(tempr, temp2[i,])
    }
  } else {
    tempr <-rbind(tempr,temp3[i,])
  }
  tempx <-rbind(tempx,tempr)
  rm(tempr)
}
tempfinal <-tempx[2:2195,]

md1<-mdy('8/11/2021')
md2<-mdy('11/3/2021')

fcurveint <-rep(NA, 3)
plot (x=0, xlim=c(md1-12,md2+28), ylim=c(0,120), xlab="Flower Date", ylab="Flower Number")
dayValues <- seq(18800, 18990, 1)
fnumatrix <- matrix(data=NA, nrow=nlevels(as.factor(tempfinal$RILID)), ncol=6)

#now we can calculate curves, plot point, and estimate total flower numebr for each RIL line


######nls estimates for curves#####

#from ft2022 code
#creating non-linear models - Gaussian,

#following model is adapted from stack-overflow
#model of normal distribution 
Gaussmodel <- function(x, theta)  { 
  m <- theta[1]; sd <- theta[2]; a <- theta[3]; b <- theta[4];
  a*exp(-0.5*((x-m)/sd)^2) + b
}

#putting into a 4-loop/iterating across all flowers


plot (x=0, xlim=c(md1-12,md2+28), ylim=c(0,120), xlab="Flower Date", ylab="Flower Number")
dayValues <- seq(18800, 18990, 1)
fnumatrix <- matrix(data=NA, nrow=nlevels(as.factor(tempfinal$RILID)), ncol=6)
meanRILftlist <- rep(NA, nlevels(as.factor(tempfinal$RILID)))
gmhold<- matrix(data=NA, nrow=nlevels(as.factor(tempfinal$RILID)), ncol=191)
control1 <- nls.control(maxiter= 100,tol=1e-02, warnOnly=TRUE)
earlyzs<-data.frame(as.numeric(c(md1-2, md1-1, md1)) , c(0,0,0))
names(earlyzs) <- c('mdate', 'fn')
latezs<-data.frame(as.numeric(c(md2, md2+1, md2+2)) , c(0,0,0))
names(latezs) <- c('mdate', 'fn')

#for loop done by hand b/c of zero gradient error w/line 24
#BE VERY CAREFUL WITH FLOWER NUMBER - estimates are very high for late flowering plants due to how model is fitted
for (i in 1:nlevels(as.factor(tempfinal$RILID))) {
  RILhold<-subset(tempfinal, as.numeric(as.factor(RILID))==i)
  if (is.na(RILhold$FlowerDate[1])){
    next
  } else {
    RILnn <- as.character(RILhold$RILID[1])
    RILcol <-as.character(RILhold$RILcolor[1])
    if (RILcol=="P"){
      coloring=3
    } else if (RILcol=="W"){
      coloring=2
    } else if (RILcol=="L"){
      coloring=1
    }
    mdate<-fn<-rep(NA, nlevels(as.factor(RILhold$FlowerDate)))
    print(nlevels(as.factor(RILhold$FlowerDate)))
    for (j in 1:nlevels(as.factor(RILhold$FlowerDate))){
      tmp1<-subset(RILhold, as.numeric(as.factor(FlowerDate))==j)
      mdate[j]<-as.numeric(tmp1$FlowerDate[1])
      fn[j]<-ave(as.numeric(tmp1$FlowerNumber[1]))
    }
    minit <- mdate[which.max(fn)] #initial mean which may move
    sinit<-(max(mdate, na.rm=T)-min(mdate, na.rm=T))/4 #rough estimate of sd
    binit <- 0; #minimizing the data start point
    ainit <- (max(fn, na.rm=T)-min(fn, na.rm=T)) #largest possible range
    gfmtx<- rbind(earlyzs, data.frame(mdate,fn), latezs)
    Gaussflowerfit <- nls(fn ~ Gaussmodel(mdate,c(m,sd,a,b)), gfmtx, start=list(m=minit, sd=sinit, a=ainit, b=binit), trace=TRUE, control=control1)    
    if (isTRUE(Gaussflowerfit$convInfo$isConv)){
      meanlft<- summary(Gaussflowerfit)$parameters["m", 1]
    } else {
      meanlft<- weighted.mean(mdate, fn) #accounts for lines that don't fit
      print("THIS IS TEMPORARY!!!!!!!!!!!!!!!!!!!!!!!!")
    }
    meanRILftlist[i]<-meanlft
    points(gfmtx$mdate, gfmtx$fn, pch=16, col=coloring, ylim=c(0,120))
    gm <-Gaussmodel(dayValues, coef(Gaussflowerfit))
    gmhold[i,]<-gm
    lines(dayValues, gm, col=coloring, lwd=2)
    mdate<- rep(NA, nrow(RILhold))
    fnumt <-sum(gm)
    tta <-append(RILnn,RILcol)
    ttb <-append(tta, fcurveint)
    ttt <- append(ttb,fnumt)
    fnumatrix[i,] <- ttt
  } 
  print(i)
  print(as.factor(RILhold$RILcolor[1]))
}
gmcv<-matrix(rep(gmhold[,1],ncol(gmhold)),nrow=nrow(gmhold))
gmnew<-gmhold-gmcv
gmnew[gmnew>65] <- 65 #this cleans the data and corrects for overestimates late
fnumatrix[,6]<-rowSums(gmnew)

#let's try a couple of other versions with different md bounds
#for highly restrictive bounds- first flower was Aug 16:
md1<-mdy('8/15/2021')
#plants likely had few to 0 flowers by October 25th
md2<-mdy('10/25/2021')
#add this code above when testing differences setting the bounds makes

flowermatrix<-as.data.frame(cbind(fnumatrix[,1:2], meanRILftlist,fnumatrix[,6]))
colnames(flowermatrix)<-c('RIL', 'color', 'meanft', 'flowersum')
greps<- cbind(flowermatrix, gmhold)

plot(as.numeric(flowermatrix$flowersum)~as.numeric(flowermatrix$meanft), col=as.numeric(as.factor(flowermatrix$color)), xlab="Peak Flower Time", ylab="Estimated Flower Number")
boxplot(as.numeric(flowermatrix$meanft)~flowermatrix$color, xlab="Line Color", ylab="Julian Date")
boxplot(as.numeric(flowermatrix$flowersum)~flowermatrix$color, ylim=c(0, 2000), xlab="Line Color", ylab="Estimated Flower Number")
plot(as.numeric(flowermatrix$flowersum)~as.numeric(flowermatrix$meanft), ylim=c(0, 1500))
abline(lm(as.numeric(flowermatrix$flowersum)~as.numeric(flowermatrix$meanft)))
summary(lm(as.numeric(flowermatrix$flowersum)~as.numeric(flowermatrix$meanft)))
summary(lm(as.numeric(flowermatrix$flowersum)~flowermatrix$color))
summary(lm(as.numeric(flowermatrix$meanft)~flowermatrix$color))
nltest<-subset(flowermatrix, flowermatrix$color!="L")

summary(lm(as.numeric(nltest$meanft)~nltest$color))
summary(lm(as.numeric(nltest$flowersum)~nltest$color))
#seems like color may be sig predictor of flower number?? 
#note low p-value though, and the fact that model has hard time predicting averages for right-biased data)

#now let's add this data to the line average data
match(linematrix$Line, nltest$RIL)#these match perfectly, minus 23
nltest2<-rbind(nltest[1:22,], nltest[24:52,]) #excl 23 (prob error in notetaking)
match(linematrix$Line, nltest2$RIL)

linematrix2<-cbind(linematrix,nltest2[,3:4])
colnames(xystore)<-c("Line","x","y")
match(linematrix2$Line, xystore[,1])

#for quick analyses, import matrices
linematrix2<-read.csv("C:/Users/jonat/OneDrive/Documents/Grad School Research Files/linematrix2.csv")

linematrix3<-cbind(linematrix2,xystore[,2:3])

mftwfe<-lmer(as.numeric(meanft)~Color+(1|x)+(1|y), data=linematrix3)
fnwfe<-lmer(as.numeric(flowersum)~Color+(1|x)+(1|y), data=linematrix3)
###

plot(as.numeric(linematrix2$Seed_Set)~as.numeric(linematrix2$Flower_Size), xlab="Flower Size", ylab="Seed Set (g)")
abline(lm(as.numeric(linematrix2$Seed_Set)~as.numeric(linematrix2$Flower_Size)))
plot(as.numeric(linematrix2$Seed_Set)~as.numeric(linematrix2$Sugar_Conc), xlab="Sugar Concentration", ylab="Seed Set(g)")
abline(lm(as.numeric(linematrix2$Seed_Set)~as.numeric(linematrix2$Sugar_Conc)))

pine<-subset(linematrix2, linematrix2$Color=="P")
wine<-subset(linematrix2, linematrix2$Color=="W")

plot(as.numeric(linematrix2$Seed_Set)~as.numeric(linematrix2$Flower_Size), xlab="Corolla Length (mm)", ylab="Seed Set (g)", col=as.numeric(as.factor(linematrix2$Color))%%2*5+1, pch=as.numeric(as.factor(linematrix2$Color))*2+17%%20)
abline(lm(as.numeric(pine$Seed_Set)~as.numeric(pine$Flower_Size)), col="magenta", lwd=2)
abline(lm(as.numeric(wine$Seed_Set)~as.numeric(wine$Flower_Size)), col="black", lwd=2)
summary(lm(as.numeric(linematrix2$Seed_Set)~as.numeric(linematrix2$Flower_Size)))

plot(as.numeric(linematrix2$Seed_Set)~as.numeric(linematrix2$Sugar_Conc), xlab="Sugar Concentration (mg/ml)", ylab="Seed Set (g)", col=as.numeric(as.factor(linematrix2$Color))%%2*5+1, pch=as.numeric(as.factor(linematrix2$Color))*2+17%%20)
abline(lm(as.numeric(pine$Seed_Set)~as.numeric(pine$Sugar_Conc)), col="magenta", lwd=2)
abline(lm(as.numeric(wine$Seed_Set)~as.numeric(wine$Sugar_Conc)), col="black", lwd=2)
summary(lm(as.numeric(linematrix2$Seed_Set)~as.numeric(linematrix2$Sugar_Conc)))

plot(as.numeric(linematrix2$Projected_Number)~as.numeric(linematrix2$Flower_Size), xlab="Corolla Length (mm)", ylab="Estimated Seed Number", col=as.numeric(as.factor(linematrix2$Color))%%2*5+1, pch=as.numeric(as.factor(linematrix2$Color))*2+17%%20)
abline(lm(as.numeric(pine$Projected_Number)~as.numeric(pine$Flower_Size)), col="magenta", lwd=2)
abline(lm(as.numeric(wine$Projected_Number)~as.numeric(wine$Flower_Size)), col="black", lwd=2)
summary(lm(as.numeric(linematrix2$Projected_Number)~as.numeric(linematrix2$Flower_Size)))

plot(as.numeric(linematrix2$Projected_Number)~as.numeric(linematrix2$Sugar_Conc), xlab="Sugar Concentration (mg/ml)", ylab="Estimated Seed Number", col=as.numeric(as.factor(linematrix2$Color))%%2*5+1, pch=as.numeric(as.factor(linematrix2$Color))*2+17%%20)
abline(lm(as.numeric(pine$Projected_Number)~as.numeric(pine$Sugar_Conc)), col="magenta", lwd=2)
abline(lm(as.numeric(wine$Projected_Number)~as.numeric(wine$Sugar_Conc)), col="black", lwd=2)
summary(lm(as.numeric(linematrix2$Projected_Number)~as.numeric(linematrix2$Sugar_Conc)))

plot(as.numeric(linematrix2$Adj_seed_weight)~as.numeric(linematrix2$Flower_Size), xlab="Corolla Length (mm)", ylab="Adjusted Seed Set (g)", col=as.numeric(as.factor(linematrix2$Color))%%2*5+1, pch=as.numeric(as.factor(linematrix2$Color))*2+17%%20)
abline(lm(as.numeric(pine$Adj_seed_weight)~as.numeric(pine$Flower_Size)), col="magenta", lwd=2)
abline(lm(as.numeric(wine$Adj_seed_weight)~as.numeric(wine$Flower_Size)), col="black", lwd=2)
summary(lm(as.numeric(linematrix2$Adj_seed_weight)~as.numeric(linematrix2$Flower_Size)))

plot(as.numeric(linematrix2$Adj_seed_weight)~as.numeric(linematrix2$Sugar_Conc), xlab="Sugar Concentration (mg/ml)", ylab="Adjusted Seed Set (g)", col=as.numeric(as.factor(linematrix2$Color))%%2*5+1, pch=as.numeric(as.factor(linematrix2$Color))*2+17%%20)
abline(lm(as.numeric(pine$Adj_seed_weight)~as.numeric(pine$Sugar_Conc)), col="magenta", lwd=2)
abline(lm(as.numeric(wine$Adj_seed_weight)~as.numeric(wine$Sugar_Conc)), col="black", lwd=2)
summary(lm(as.numeric(linematrix2$Adj_seed_weight)~as.numeric(linematrix2$Sugar_Conc)))

plot(as.numeric(linematrix2$Adj_number)~as.numeric(linematrix2$Flower_Size), xlab="Corolla Length (mm)", ylab="Adjusted Projected Seed Number", col=as.numeric(as.factor(linematrix2$Color))%%2*5+1, pch=as.numeric(as.factor(linematrix2$Color))*2+17%%20)
abline(lm(as.numeric(pine$Adj_number)~as.numeric(pine$Flower_Size)), col="magenta", lwd=2)
abline(lm(as.numeric(wine$Adj_number)~as.numeric(wine$Flower_Size)), col="black", lwd=2)
summary(lm(as.numeric(linematrix2$Adj_number)~as.numeric(linematrix2$Flower_Size)))

plot(as.numeric(linematrix2$Adj_number)~as.numeric(linematrix2$Sugar_Conc), xlab="Sugar Concentration (mg/ml)", ylab="Adjusted Projected Seed Number", col=as.numeric(as.factor(linematrix2$Color))%%2*5+1, pch=as.numeric(as.factor(linematrix2$Color))*2+17%%20)
abline(lm(as.numeric(pine$Adj_number)~as.numeric(pine$Sugar_Conc)), col="magenta", lwd=2)
abline(lm(as.numeric(wine$Adj_number)~as.numeric(wine$Sugar_Conc)), col="black", lwd=2)
summary(lm(as.numeric(linematrix2$Adj_number)~as.numeric(linematrix2$Sugar_Conc)))

plot(as.numeric(linematrix2$Resids)~as.numeric(linematrix2$Flower_Size), xlab="Corolla Length (mm)", ylab="Residual difference",col=as.numeric(as.factor(linematrix2$Color))%%2*5+1, pch=as.numeric(as.factor(linematrix2$Color))*2+17%%20)
abline(lm(as.numeric(pine$Resids)~as.numeric(pine$Flower_Size)), col="magenta", lwd=2)
abline(lm(as.numeric(wine$Resids)~as.numeric(wine$Flower_Size)), col="black", lwd=2)
summary(lm(as.numeric(linematrix2$Resids)~as.numeric(linematrix2$Flower_Size)))

plot(as.numeric(linematrix2$Resids)~as.numeric(linematrix2$Sugar_Conc), xlab="Sugar Concentration (mg/ml)", ylab="Residual difference", col=as.numeric(as.factor(linematrix2$Color))%%2*5+1, pch=as.numeric(as.factor(linematrix2$Color))*2+17%%20)
abline(lm(as.numeric(pine$Resids)~as.numeric(pine$Sugar_Conc)), col="magenta", lwd=2)
abline(lm(as.numeric(wine$Resids)~as.numeric(wine$Sugar_Conc)), col="black", lwd=2)
summary(lm(as.numeric(linematrix2$Resids)~as.numeric(linematrix2$Sugar_Conc)))


write.csv(linematrix2, "C:/Users/jonat/OneDrive/Documents/Grad School Research Files/linematrix2.csv")
pline3<-subset(linematrix3, linematrix$Color=="P")
wline3<-subset(linematrix3, linematrix$Color=="W")

Anova(lmer(Seed_Set~Color*Sugar_Conc+(1|x)+(1|y), data=linematrix3), type="II")
Anova(lmer(Projected_Number~Color*Sugar_Conc+(1|x)+(1|y), data=linematrix3), type="II")
Anova(lmer(Adj_seed_weight~Color*Sugar_Conc+(1|x)+(1|y), data=linematrix3), type="II")
Anova(lmer(Adj_number~Color*Sugar_Conc+(1|x)+(1|y), data=linematrix3), type="II")
Anova(lmer(Resids~Color*Sugar_Conc+(1|x)+(1|y), data=linematrix3), type="II")

#for individual colors
Anova(lmer(Seed_Set~Sugar_Conc+(1|x)+(1|y), data=pline3), type="II")
Anova(lmer(Projected_Number~Sugar_Conc+(1|x)+(1|y), data=pline3), type="II")
Anova(lmer(Adj_seed_weight~Sugar_Conc+(1|x)+(1|y), data=pline3), type="II")
Anova(lmer(Adj_number~Sugar_Conc+(1|x)+(1|y), data=pline3), type="II")
Anova(lmer(Resids~Sugar_Conc+(1|x)+(1|y), data=pline3), type="II")
#sig (seed set) for pink lines
Anova(lmer(Seed_Set~Sugar_Conc+(1|x)+(1|y), data=wline3), type="II")
Anova(lmer(Projected_Number~Sugar_Conc+(1|x)+(1|y), data=wline3), type="II")
Anova(lmer(Adj_seed_weight~Sugar_Conc+(1|x)+(1|y), data=wline3), type="II")
Anova(lmer(Adj_number~Sugar_Conc+(1|x)+(1|y), data=wline3), type="II")
Anova(lmer(Resids~Sugar_Conc+(1|x)+(1|y), data=wline3), type="II")
#n.s.
Anova(lmer(Seed_Set~Color*Sugar_Conc+(1|x)+(1|y), data=linematrix3), type=3, contrasts=list(topic=contr.sum, sys=contr.sum))
Anova(lmer(Projected_Number~Color*Sugar_Conc+(1|x)+(1|y), data=linematrix3), type=3, contrasts=list(topic=contr.sum, sys=contr.sum))
Anova(lmer(Adj_seed_weight~Color*Sugar_Conc+(1|x)+(1|y), data=linematrix3), type=3, contrasts=list(topic=contr.sum, sys=contr.sum))
Anova(lmer(Adj_number~Color*Sugar_Conc+(1|x)+(1|y), data=linematrix3), type=3, contrasts=list(topic=contr.sum, sys=contr.sum))
Anova(lmer(Resids~Color*Sugar_Conc+(1|x)+(1|y), data=linematrix3), type=3, contrasts=list(topic=contr.sum, sys=contr.sum))
##
##
##
Anova(lmer(Seed_Set~Color*Flower_Size+(1|x)+(1|y), data=linematrix3), type="II")
Anova(lmer(Projected_Number~Color*Flower_Size+(1|x)+(1|y), data=linematrix3), type="II")
Anova(lmer(Adj_seed_weight~Color*Flower_Size+(1|x)+(1|y), data=linematrix3), type="II")
Anova(lmer(Adj_number~Color*Flower_Size+(1|x)+(1|y), data=linematrix3), type="II")
Anova(lmer(Resids~Color*Flower_Size+(1|x)+(1|y), data=linematrix3), type="II")
#for individual colors
Anova(lmer(Seed_Set~Flower_Size+(1|x)+(1|y), data=pline3), type="II")
Anova(lmer(Projected_Number~Flower_Size+(1|x)+(1|y), data=pline3), type="II")
Anova(lmer(Adj_seed_weight~Flower_Size+(1|x)+(1|y), data=pline3), type="II")
Anova(lmer(Adj_number~Flower_Size+(1|x)+(1|y), data=pline3), type="II")
Anova(lmer(Resids~Flower_Size+(1|x)+(1|y), data=pline3), type="II")
#sig (seed set) or sugg (flower number) for pink lines
Anova(lmer(Seed_Set~Flower_Size+(1|x)+(1|y), data=wline3), type="II")
Anova(lmer(Projected_Number~Flower_Size+(1|x)+(1|y), data=wline3), type="II")
Anova(lmer(Adj_seed_weight~Flower_Size+(1|x)+(1|y), data=wline3), type="II")
Anova(lmer(Adj_number~Flower_Size+(1|x)+(1|y), data=wline3), type="II")
Anova(lmer(Resids~Flower_Size+(1|x)+(1|y), data=wline3), type="II")
#n.s. for white
Anova(lmer(Seed_Set~Color*Flower_Size+(1|x)+(1|y), data=linematrix3), type=3, contrasts=list(topic=contr.sum, sys=contr.sum))
Anova(lmer(Projected_Number~Color*Flower_Size+(1|x)+(1|y), data=linematrix3), type=3, contrasts=list(topic=contr.sum, sys=contr.sum))
Anova(lmer(Adj_seed_weight~Color*Flower_Size+(1|x)+(1|y), data=linematrix3), type=3, contrasts=list(topic=contr.sum, sys=contr.sum))
Anova(lmer(Adj_number~Color*Flower_Size+(1|x)+(1|y), data=linematrix3), type=3, contrasts=list(topic=contr.sum, sys=contr.sum))
Anova(lmer(Resids~Color*Flower_Size+(1|x)+(1|y), data=linematrix3), type=3, contrasts=list(topic=contr.sum, sys=contr.sum))
#adjusted seed number or projected seed number way to go?

plot(as.numeric(linematrix2$Projected_Number)~as.numeric(linematrix2$Flower_Size), xlab="Flower Size", ylab="Projected Seed Number", col=as.numeric(as.factor(linematrix2$Color)))
abline(lm(as.numeric(pine$Projected_Number)~as.numeric(pine$Flower_Size)), col="black")
abline(lm(as.numeric(wine$Projected_Number)~as.numeric(wine$Flower_Size)), col="red")
plot(as.numeric(linematrix2$Projected_Number)~as.numeric(linematrix2$Sugar_Conc), xlab="Sugar Conc", ylab="Projected Seed Number", col=as.numeric(as.factor(linematrix2$Color)))
abline(lm(as.numeric(pine$Projected_Number)~as.numeric(pine$Sugar_Conc)), col="black")
abline(lm(as.numeric(wine$Projected_Number)~as.numeric(wine$Sugar_Conc)), col="red")

linematrix2$meanft<-as.numeric(linematrix2$meanft)
linematrix2$flowersum<-as.numeric(linematrix2$flowersum)

#excluding lines with inflated flower number:
slinem<-subset(linematrix2, linematrix2$flowersum<1100)

plot(linematrix2$Seed_Set~linematrix2$flowersum, xlab="Estimated Flower Number", ylab="Seed set", col=as.numeric(as.factor(linematrix2$Color)))
summary(lm(linematrix2$Seed_Set~linematrix2$flowersum))
abline(lm(linematrix2$Seed_Set~linematrix2$flowersum))

plot(linematrix2$Projected_Number~linematrix2$flowersum, xlab="Estimated Flower Number", ylab="Seed set", col=as.numeric(as.factor(linematrix2$Color)))
summary(lm(linematrix2$Projected_Number~linematrix2$flowersum))
abline(lm(linematrix2$Projected_Number~linematrix2$flowersum))


plot(slinem$Seed_Set~slinem$flowersum,  xlab="Estimated Flower Number", ylab="Seed set", col=as.numeric(as.factor(linematrix2$Color)))
summary(lm(slinem$Seed_Set~slinem$flowersum))
abline(lm(slinem$Seed_Set~slinem$flowersum))

plot(slinem$Projected_Number~slinem$flowersum,  xlab="Estimated Flower Number", ylab="Seed set", col=as.numeric(as.factor(linematrix2$Color)))
summary(lm(slinem$Projected_Number~slinem$flowersum))
abline(lm(slinem$Projected_Number~slinem$flowersum))
#taking out outliers that appear to be overestimating model fit shows sig predic power

plot(linematrix2$Seed_Set~linematrix2$meanft, xlab="Julian Date", ylab="Seed set", col=as.numeric(as.factor(linematrix2$Color)))
abline(lm(as.numeric(pine$Seed_Set)~as.numeric(pine$meanft)), col="black")
abline(lm(as.numeric(wine$Seed_Set)~as.numeric(wine$meanft)), col="red")
Anova(aov(Seed_Set~Color*meanft, data=linematrix2), type="III")


plot(linematrix2$Projected_Number~linematrix2$meanft, xlab="Julian Date", ylab="Seed set", col=as.numeric(as.factor(linematrix2$Color)))
abline(lm(as.numeric(pine$Projected_Number)~as.numeric(pine$meanft)), col="black")
abline(lm(as.numeric(wine$Projected_Number)~as.numeric(wine$meanft)), col="red")

Anova(aov(Projected_Number~Color*meanft, data=linematrix2), type="III")

summary(lm(linematrix2$Projected_Number~linematrix2$meanft))
summary(lm(slinem$Projected_Number~slinem$meanft))

plot(slinem$Seed_Set~slinem$meanft)
summary(lm(slinem$Seed_Set~slinem$meanft))
abline(lm(slinem$Seed_Set~slinem$meanft))
#but it doesn't efect n.s. result from ft

plot(linematrix2$flowersum~linematrix2$Flower_Size)
plot(slinem$flowersum~slinem$Flower_Size)
summary(lm(linematrix2$flowersum~linematrix2$Flower_Size))
summary(lm(slinem$flowersum~slinem$Flower_Size))
abline(lm(linematrix2$flowersum~linematrix2$Flower_Size))
plot(linematrix2$flowersum~linematrix2$Sugar_Conc)
plot(slinem$flowersum~slinem$Sugar_Conc)
summary(lm(linematrix2$flowersum~linematrix2$Sugar_Conc))
summary(lm(slinem$flowersum~slinem$Sugar_Conc))
abline(lm(linematrix2$flowersum~linematrix2$Sugar_Conc))

plot(linematrix2$meanft~linematrix2$Flower_Size)
summary(lm(linematrix2$meanft~linematrix2$Flower_Size))
abline(lm(linematrix2$meanft~linematrix2$Flower_Size))
plot(linematrix2$meanft~linematrix2$Sugar_Conc)
summary(lm(linematrix2$meanft~linematrix2$Sugar_Conc))
abline(lm(linematrix2$meanft~linematrix2$Sugar_Conc))

boxplot(slinem$flowersum~slinem$Color)
summary(lm(slinem$flowersum~slinem$Color))
#now for selection models
#with seed set
sscsc<-lmer(as.numeric(Seed_Set)~as.numeric(as.factor(Color))*as.numeric(Sugar_Conc)+(1|x)+(1|y), data=linematrix3) #error in processing color
sscfz<-lmer(as.numeric(Seed_Set)~as.numeric(as.factor(Color))*as.numeric(Flower_Size)+(1|x)+(1|y), data=linematrix3) #error in processing color


#first as a lm, then as a glm
sscft<-lmer(as.numeric(Seed_Set)~as.numeric(as.factor(Color))*as.numeric(meanft)+(1|x)+(1|y), data=linematrix3) #error in processing color
sscfm<-lmer(as.numeric(Seed_Set)~as.numeric(as.factor(Color))*as.numeric(flowersum)+(1|x)+(1|y), data=linematrix3) #error in processing color
plot(as.numeric(Seed_Set)~as.numeric(flowersum),data=linematrix3, col=as.numeric(as.factor(linematrix2$Color))%%2*5+1, pch=as.numeric(as.factor(linematrix2$Color))*2+17%%20)
abline(lm(as.numeric(pine$Seed_Set)~as.numeric(pine$flowersum)), col="magenta", lwd=2)
abline(lm(as.numeric(wine$Seed_Set)~as.numeric(wine$flowersum)), col="black", lwd=2)
#for pink limbed lines, no effect of flower prediction on seed set, but white limbed lines there are effects?
#slscfm<-lm(Seed_Set~Color*flowersum, data=slinem)
ssfzft<-lmer(as.numeric(Seed_Set)~as.numeric(Flower_Size)*as.numeric(meanft)+(1|x)+(1|y), data=linematrix3) #error in processing color
ssfzfm<-lmer(as.numeric(Seed_Set)~as.numeric(Flower_Size)*as.numeric(flowersum)+(1|x)+(1|y), data=linematrix3) #error in processing color
#flower size by flower number may be a +joint effect?
#plot(as.numeric(Seed_Set)~as.numeric(flowersum),data=linematrix3, col=round((as.numeric(Flower_Size)-20)/3),pch=as.numeric(round((as.numeric(Flower_Size)-20)/3)*2+17%%20))
#abline(lm(as.numeric(Seed_Set)~as.numeric(flowersum), data=linematrix3))

ssscft<-lmer(as.numeric(Seed_Set)~as.numeric(Sugar_Conc)*as.numeric(meanft)+(1|x)+(1|y), data=linematrix3)
ssscfm<-lmer(as.numeric(Seed_Set)~as.numeric(Sugar_Conc)*as.numeric(flowersum)+(1|x)+(1|y), data=linematrix3)
ssftfm<-lmer(as.numeric(Seed_Set)~as.numeric(meanft)*as.numeric(flowersum)+(1|x)+(1|y), data=linematrix3)
#this is a significant (-) joint prediction of small effect size, both independently have positive effects

ssall<-lmer(as.numeric(Seed_Set)~as.numeric(as.factor(Color))*as.numeric(meanft)+as.numeric(as.factor(Color))*as.numeric(flowersum)+as.numeric(Flower_Size)*as.numeric(meanft)+as.numeric(Flower_Size)*as.numeric(flowersum)
          +as.numeric(Sugar_Conc)*as.numeric(meanft)+as.numeric(Sugar_Conc)*as.numeric(flowersum)+as.numeric(flowersum)*as.numeric(meanft)+as.numeric(as.factor(Color))*as.numeric(Flower_Size)+
            as.numeric(as.factor(Color))*as.numeric(Sugar_Conc)+as.numeric(Sugar_Conc)*as.numeric(Flower_Size)+(1|x)+(1|y), data=linematrix3)
ssallnfs<-lmer(as.numeric(Seed_Set)~as.numeric(as.factor(Color))*as.numeric(meanft)+as.numeric(as.factor(Color))+as.numeric(Flower_Size)*as.numeric(meanft)
             +as.numeric(Sugar_Conc)*as.numeric(meanft)+as.numeric(as.factor(Color))*as.numeric(Flower_Size)+
               as.numeric(as.factor(Color))*as.numeric(Sugar_Conc)+as.numeric(Sugar_Conc)*as.numeric(Flower_Size)+(1|x)+(1|y), data=linematrix3)
sscnfs<-lmer(as.numeric(Seed_Set)~as.numeric(as.factor(Color))*as.numeric(meanft)+as.numeric(as.factor(Color))*as.numeric(Flower_Size)+as.numeric(as.factor(Color))*as.numeric(Sugar_Conc)+(1|x)+(1|y), data=linematrix3)
sscnft<-lmer(as.numeric(Seed_Set)~as.numeric(as.factor(Color))*as.numeric(Flower_Size)+as.numeric(as.factor(Color))*as.numeric(Sugar_Conc)+(1|x)+(1|y), data=linematrix3)

ss5w<-lm(as.numeric(Seed_Set)~Color*meanft*flowersum*Flower_Size*Sugar_Conc, data=linematrix2)
ss4w<-lm(as.numeric(Seed_Set)~Color*meanft*Flower_Size*Sugar_Conc, data=linematrix2)

summary(ssall)
nullt<-lmer(as.numeric(Seed_Set)~(1|x)+(1|y), data=linematrix3)
ftstep<-step(ssall, direction=c("both"), scope=formula(nullt))
ftstep2<-step(ssall, direction=c("backward"), scope=formula(nullt))
ftstep3<-step(ssallnfs, direction=c("forward"), scope=formula(nullt))
ftstep33<-step(nullt, direction=c("forward"), scope=formula(ssallnfs))
ftstep4<-step(sscnfs, direction=c("both"), scope=formula(nullt))
ftstep5<-step(nullt, direction=c("both"), scope=formula(sscnfs))
ftstep6<-step(sscnft, direction=c("both"), scope=formula(nullt))
ftstep7<-step(nullt, direction=c("both"), scope=formula(sscnft))
#suggest that no effects are significant with random effects of average x and y coordinates!

psallnfs<-lm(Projected_Number~Color*meanft+Flower_Size*meanft +Sugar_Conc*meanft+Color*Flower_Size+
               Color*Sugar_Conc+Sugar_Conc*Flower_Size, data=linematrix2)
pscnfs<-lm(Projected_Number~Color+Flower_Size+Sugar_Conc+meanft+Color*meanft+Color*Flower_Size+Color*Sugar_Conc, data=linematrix2)
nullp<-lm(Projected_Number~0, data=linematrix2)
fstep<-step(psallnfs, direction=c("both"), scope=formula(nullp))
fstep2<-step(nullp, direction=c("both"), scope=formula(psallnfs))
fstep4<-step(pscnfs, direction=c("both"), scope=formula(nullp))
fstep5<-step(nullp, direction=c("both"), scope=formula(pscnfs))

#glm step function?
resft<-as.numeric(linematrix3$meanft)-18878

gssall<-glmer(as.numeric(Seed_Set)~as.numeric(as.factor(Color))*resft+
                as.numeric(as.factor(Color))*log(as.numeric(flowersum))+
                as.numeric(Flower_Size)*resft+
              as.numeric(Sugar_Conc)*resft+
              as.numeric(as.factor(Color))*as.numeric(Flower_Size)+
              as.numeric(as.factor(Color))*as.numeric(Sugar_Conc)+
              as.numeric(Sugar_Conc)*as.numeric(Flower_Size)
            +(1|x)+(1|y), data=linematrix3, family=Gamma(link=log), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), verbose = 100)
nullg<-glmer(as.numeric(Seed_Set)~(1|x)+(1|y), family=Gamma(link=log),data=linematrix3)
gstep<-step(gssall, direction=c("both"), scope=formula(nullg))
gstep2<-step(nullg, direction=c("both"), scope=formula(gssall))



#for glm
fts3wg<-glmer(as.numeric(Seed_Set)~as.numeric(as.factor(Color))*as.numeric(Flower_Size)+as.numeric(Sugar_Conc)+(1|x)+(1|y), data=linematrix3, family=Gamma)
Anova(fts3wg, type="III")
fts3wgp<-glmer(as.integer(Projected_Number)~as.numeric(as.factor(Color))*as.numeric(Flower_Size)+as.numeric(Sugar_Conc)+(1|x)+(1|y), data=linematrix3, family=poisson)
Anova(fts3wg, type="III")
f3wgg<-glmer(as.integer(Seed_Set)~as.numeric(as.factor(Color))*as.numeric(Flower_Size)+as.numeric(as.factor(Color))*as.numeric(Sugar_Conc)+(1|x)+(1|y), data=linematrix3, family=Gamma)
Anova(f3wgg, type="III") #also overconfident
f3wgl<-lmer(as.integer(Seed_Set)~as.numeric(as.factor(Color))*as.numeric(Flower_Size)+as.numeric(as.factor(Color))*as.numeric(Sugar_Conc)+(1|x)+(1|y), data=linematrix3)
Anova(f3wgl, type="III")
f2sg<-lmer(as.integer(Seed_Set)~as.numeric(as.factor(Color))*as.numeric(Sugar_Conc)+(1|x)+(1|y), data=linematrix3)
Anova(f2sg, type="III")
f2fz<-lmer(as.integer(Seed_Set)~as.numeric(as.factor(Color))*as.numeric(Flower_Size)+(1|x)+(1|y), data=linematrix3)
Anova(f2fz, type="III")
f3wall<-lmer(as.integer(Seed_Set)~as.numeric(as.factor(Color))*as.numeric(Flower_Size)*as.numeric(Sugar_Conc)+(1|x)+(1|y), data=linematrix3)
Anova(f3wall, type="III")
f3w<-lmer(as.integer(Seed_Set)~as.numeric(as.factor(Color))+(1|x)+(1|y), data=linematrix3)
Anova(f3w, type="III")
boxplot(as.numeric(linematrix3$Seed_Set)~linematrix3$Color, xlab="Line Color", ylab="Line Average of Seed Set(g)")

f3wgp<-glmer(as.integer(Projected_Number)~as.numeric(as.factor(Color))*as.numeric(Flower_Size)+as.numeric(as.factor(Color))*as.numeric(Sugar_Conc)+(1|x)+(1|y), data=linematrix3, family=poisson)
Anova(fts3wgp, type="III")
f3wgl<-lmer(as.integer(Projected_Number)~as.numeric(as.factor(Color))*as.numeric(Flower_Size)+as.numeric(as.factor(Color))*as.numeric(Sugar_Conc)+(1|x)+(1|y), data=linematrix3)
Anova(f3wgl, type="III")
f2sg<-lmer(as.integer(Projected_Number)~as.numeric(as.factor(Color))*as.numeric(Sugar_Conc)+(1|x)+(1|y), data=linematrix3)
Anova(f2sg, type="III")
f2fz<-lmer(as.integer(Projected_Number)~as.numeric(as.factor(Color))*as.numeric(Flower_Size)+(1|x)+(1|y), data=linematrix3)
Anova(f2fz, type="III")
f3wall<-lmer(as.integer(Projected_Number)~as.numeric(as.factor(Color))*as.numeric(Flower_Size)*as.numeric(Sugar_Conc)+(1|x)+(1|y), data=linematrix3)
Anova(f3wall, type="III")
f3wallp<-glmer(as.integer(Projected_Number)~as.numeric(as.factor(Color))*as.numeric(Flower_Size)*as.numeric(Sugar_Conc)+(1|x)+(1|y), data=linematrix3, family=poisson)
Anova(f3wallp, type="III")
#poisson glm models seem to be highly overconfident

Anova(aov(Seed_Set ~ Color*Flower_Size*Sugar_Conc+Color, data =linematrix2), type="III")
Anova(aov(Seed_Set ~ Color*Flower_Size+Sugar_Conc+Color, data =linematrix2), type="III")
Anova(aov(Projected_Number ~ Color*Flower_Size+Sugar_Conc+Color, data =linematrix2), type="III")
Anova(aov(ftstepx))
Anova(aov(Seed_Set ~ Color*Flower_Size+Color, data =linematrix2), type="III")
Anova(aov(Seed_Set ~ Color*Sugar_Conc+Color, data =linematrix2), type="III")


Anova(aov(Seed_Set ~ meanft*Sugar_Conc*Flower_Size*Color, data =linematrix2), type="III")
Anova(aov(Seed_Set ~ meanft*Sugar_Conc*Flower_Size+Color, data =linematrix2), type="III")
Anova(aov(Seed_Set ~ Sugar_Conc*Flower_Size+Color, data =linematrix2), type="III")


bartlett.test(Seed_Set~ Color, data =linematrix2)
shapiro.test(resid(aov(Seed_Set~ Color, data =linematrix2)))
bartlett.test(Projected_Number~ Color, data =linematrix2)
shapiro.test(resid(aov(Projected_Number~ Color, data =linematrix2)))

#emmeans_test(data = linematrix2, formula = Seed_Set~ meanft+Sugar_Conc+Flower_Size, covariate = Color, p.adjust.method = "fdr")



#seems like flower number by color is only sig predictor in 5 way comparison
#susgar conc is sig in ANOVA, but interestingly color x others disappear

#with residuals
rsscft<-lm(Resids~Color*meanft, data=linematrix2)
rsscfm<-lm(Resids~Color*flowersum, data=linematrix2)
rssfzft<-lm(Resids~Flower_Size*meanft, data=linematrix2)
rssfzfm<-lm(Resids~Flower_Size*flowersum, data=linematrix2)
rssscft<-lm(Resids~Sugar_Conc*meanft, data=linematrix2)
rssscfm<-lm(Resids~Sugar_Conc*flowersum, data=linematrix2)
rssftfm<-lm(Resids~flowersum*meanft, data=linematrix2)
rssall<-lm(Resids~Color*meanft+Color*flowersum+Flower_Size*meanft+Flower_Size*flowersum
           +Sugar_Conc*meanft+Sugar_Conc*flowersum+flowersum*meanft+Color*Flower_Size+
             Color*Sugar_Conc+Sugar_Conc*Flower_Size, data=linematrix2)
rss5w<-lm(Resids~Color*meanft*flowersum*Flower_Size*Sugar_Conc, data=linematrix2)
summary(rssall)


#with adjusted:
asscft<-lm(Adj_seed_weight~Color*meanft, data=linematrix2)
asscfm<-lm(Adj_seed_weight~Color*flowersum, data=linematrix2)
assfzft<-lm(Adj_seed_weight~Flower_Size*meanft, data=linematrix2)
assfzfm<-lm(Adj_seed_weight~Flower_Size*flowersum, data=linematrix2)
assscft<-lm(Adj_seed_weight~Sugar_Conc*meanft, data=linematrix2)
assscfm<-lm(Adj_seed_weight~Sugar_Conc*flowersum, data=linematrix2)
assftfm<-lm(Adj_seed_weight~flowersum*meanft, data=linematrix2)
assall<-lm(Adj_seed_weight~Color*meanft+Color*flowersum+Flower_Size*meanft+Flower_Size*flowersum
           +Sugar_Conc*meanft+Sugar_Conc*flowersum+flowersum*meanft+Color*Flower_Size+
             Color*Sugar_Conc+Sugar_Conc*Flower_Size, data=linematrix2)
ass5w<-lm(Adj_seed_weight~Color*meanft*flowersum*Flower_Size*Sugar_Conc, data=linematrix2)
summary(assall)
#all 3 effectively show the same patterns
#whatabout including line as a random effect?

sel3t<-lmer(Seed_Set~Color+Color*Flower_Size*Sugar_Conc+(1|Line), data=linematrix2)



################
#Again, this next section is likely repetitive code, can probably skip
################


#now let's consider other predictors of seed set: flower size
summary(lm(as.numeric(nolac$Total.seed.weight)~nolac$RILcolor+nolac$Cor_Len_Var))

summary(aov(as.numeric(nolac$Total.seed.weight)~nolac$RILcolor+nolac$Cor_Len_Var))
nestws<- aov(as.numeric(nolac$Total.seed.weight)~(nolac$RILcolor+nolac$Cor_Len_Var)/factor(nolac$RILID)) #RIL ID is highly determinant of seed set, color merely suggestive
summary(nestws)
#w/block
nestwbs<- aov(as.numeric(nolac$Total.seed.weight)~factor(nolac$blockID)+(nolac$RILcolor+nolac$Cor_Len_Var)/factor(nolac$RILID), na.rm=T) #RIL ID is highly determinant of seed set, color merely suggestive
summary(nestwbs)#in a nested model, flower size seems to be a strong predictor of seed set 
plot(nolac$Total.seed.weight~nolac$Cor_Len_Var, col=as.numeric(as.factor(nolac$RILcolor))-1, xlab="Corolla Length (mm)", ylab="Seed set (g)")
psizelm<-lm(as.numeric(pdata$Total.seed.weight)~as.numeric(pdata$Cor_Len_Var))
wsizelm<-lm(as.numeric(wdata$Total.seed.weight)~as.numeric(wdata$Cor_Len_Var))
abline(psizelm)
abline(wsizelm, col="red")
#joint
nestwbs<- aov(as.numeric(nolac$Total.seed.weight)~factor(nolac$blockID)+(nolac$RILcolor*nolac$Cor_Len_Var)/factor(nolac$RILID), na.rm=T) #RIL ID is highly determinant of seed set, color merely suggestive
summary(nestwbs)#in a nested model, flower size seems to be a strong predictor of seed set 
#joint highly predictive


#let's try a different predictor: sugar concentration
summary(lm(as.numeric(nolac$Total.seed.weight)~nolac$RILcolor+nolac$Sugar_UND))
summary(aov(as.numeric(nolac$Total.seed.weight)~nolac$RILcolor+nolac$Sugar_UND))
nestwsg<- aov(as.numeric(nolac$Total.seed.weight)~(nolac$RILcolor+nolac$Sugar_UND)/factor(nolac$RILID)) #RIL ID is highly determinant of seed set, color merely suggestive
summary(nestwsg)
nestwbsg<- aov(as.numeric(nolac$Total.seed.weight)~factor(nolac$blockID)+(nolac$RILcolor+nolac$Sugar_UND)/factor(nolac$RILID), na.rm=T) #RIL ID is highly determinant of seed set, color merely suggestive
summary(nestwbsg) #sugar concentration also seems to be a strong predictor of seed set
plot(nolac$Total.seed.weight~nolac$Sugar_UND, col=as.numeric(as.factor(nolac$RILcolor))-1, xlab="Sugar concentration", ylab="Seed set")#positive trend seems a lot less clear here...
psuglm<-lm(as.numeric(pdata$Total.seed.weight)~as.numeric(pdata$Sugar_UND))
wsuglm<-lm(as.numeric(wdata$Total.seed.weight)~as.numeric(wdata$Sugar_UND))
abline(psuglm)
abline(wsuglm, col="red")

#let's try interactive effect
nestwbsgint<- aov(as.numeric(nolac$Total.seed.weight)~factor(nolac$blockID)+(nolac$RILcolor*nolac$Sugar_UND)/factor(nolac$RILID), na.rm=T) #RIL ID is highly determinant of seed set, color merely suggestive
summary(nestwbsgint) #interactive effects appears significant, but unsure how to parse this


#one that shouldn't be predictive: seed length
summary(lm(as.numeric(nolac$Total.seed.weight)~nolac$RILcolor+nolac$Seed.Length))
summary(aov(as.numeric(nolac$Total.seed.weight)~nolac$RILcolor+nolac$Seed.Length))
nestwl<- aov(as.numeric(nolac$Total.seed.weight)~(nolac$RILcolor+nolac$Seed.Length)/factor(nolac$RILID)) #RIL ID is highly determinant of seed set, color merely suggestive
summary(nestwl)
nestwlb<- aov(as.numeric(nolac$Total.seed.weight)~factor(nolac$blockID)+(nolac$RILcolor+nolac$Seed.Length)/factor(nolac$RILID), na.rm=T) #RIL ID is highly determinant of seed set, color merely suggestive
summary(nestwlb) 
plot(nolac$Total.seed.weight~nolac$Seed.Length, col=as.numeric(as.factor(nolac$RILcolor))-1)#as expected, no predictive power
psdllm<-lm(as.numeric(pdata$Total.seed.weight)~as.numeric(pdata$Seed.Length))
wsdllm<-lm(as.numeric(wdata$Total.seed.weight)~as.numeric(wdata$Seed.Length))
abline(psdllm)
abline(wsdllm, col="red")

#now let's try making a model that incorporates all relavent predictors
#predictors: color, line, early growth, shoot number, herbiv, corlen, nectar volume, sugar conc, seed mass, length, width, growth first second third
summary(aov(as.numeric(nolac$Total.seed.weight)~factor(nolac$blockID)+(nolac$RILcolor*nolac$Early.Growth..cm.)/factor(nolac$RILID), na.rm=T)) #RIL ID is highly determinant of seed set, color merely suggestive
#eg predcitive, not joint
summary(aov(as.numeric(nolac$Total.seed.weight)~factor(nolac$blockID)+(nolac$RILcolor*nolac$Shoot..)/factor(nolac$RILID), na.rm=T)) #RIL ID is highly determinant of seed set, color merely suggestive
#shoot numb predict, not joint
summary(aov(as.numeric(nolac$Total.seed.weight)~factor(nolac$blockID)+(nolac$RILcolor*nolac$Herbiv)/factor(nolac$RILID), na.rm=T)) #RIL ID is highly determinant of seed set, color merely suggestive
#herbiv joint predictive, not indiv
summary(aov(as.numeric(nolac$Total.seed.weight)~factor(nolac$blockID)+(nolac$RILcolor*nolac$Nectar_UNC)/factor(nolac$RILID), na.rm=T)) #RIL ID is highly determinant of seed set, color merely suggestive
#nectar vol predict, not joint
summary(aov(as.numeric(nolac$Total.seed.weight)~factor(nolac$blockID)+(nolac$RILcolor*nolac$Seed.Mass)/factor(nolac$RILID), na.rm=T)) #RIL ID is highly determinant of seed set, color merely suggestive
plot(nolac$Total.seed.weight~nolac$Seed.Mass, col=as.numeric(as.factor(nolac$RILcolor))-1, xlab="Averge Mass per seed within line(g)", ylab="Field Plant's Seed Weight (g)")
psdmlm<-lm(as.numeric(pdata$Total.seed.weight)~as.numeric(pdata$Seed.Mass))
wsdmlm<-lm(as.numeric(wdata$Total.seed.weight)~as.numeric(wdata$Seed.Mass))
abline(psdmlm)
abline(wsdmlm, col="red")
#joint highly predictive, not indiv seed mass
summary(aov(as.numeric(nolac$Total.seed.weight)~factor(nolac$blockID)+(nolac$RILcolor*nolac$Seed.Length)/factor(nolac$RILID), na.rm=T)) #RIL ID is highly determinant of seed set, color merely suggestive
#joint predictive, not indiv seed length
summary(aov(as.numeric(nolac$Total.seed.weight)~factor(nolac$blockID)+(nolac$RILcolor*nolac$First.Internode)/factor(nolac$RILID), na.rm=T)) #RIL ID is highly determinant of seed set, color merely suggestive
#color sig, FI sug, joint sig
summary(aov(as.numeric(nolac$Total.seed.weight)~factor(nolac$blockID)+(nolac$RILcolor*nolac$Second.Internode)/factor(nolac$RILID), na.rm=T)) #RIL ID is highly determinant of seed set, color merely suggestive
#color + joint sig, not indiv 2nd
summary(aov(as.numeric(nolac$Total.seed.weight)~factor(nolac$blockID)+(nolac$RILcolor*nolac$Third.Internode)/factor(nolac$RILID), na.rm=T)) #RIL ID is highly determinant of seed set, color merely suggestive
#all sig;note internode measures are sparse
plot(nolac$Total.seed.weight~nolac$Third.Internode, col=as.numeric(as.factor(nolac$RILcolor))-1)

#combined w/internodes
summary(aov(as.numeric(nolac$Total.seed.weight)~factor(nolac$blockID)+(nolac$Shoot..+nolac$RILcolor*nolac$Early.Growth..cm.+nolac$RILcolor*nolac$Herbiv+nolac$RILcolor*nolac$Cor_Len_Var+nolac$RILcolor*nolac$Sugar_UND+nolac$RILcolor*nolac$Seed.Length+nolac$Nectar_UNC+nolac$RILcolor*nolac$Seed.Mass+nolac$RILcolor*nolac$Third.Internode+nolac$RILcolor*nolac$Second.Internode+nolac$RILcolor*nolac$First.Internode)/factor(nolac$RILID), na.rm=T)) #RIL ID is highly determinant of seed set, color merely suggestive
#with all factors, following is sig, in order: Shoot Number(eg), color x sugar conc, Flower size, seed length, color; sugg:first in, 2nd in, color x flower size, early growth, seed mass

#now w/o internodes, as they limit the data used:
summary(aov(as.numeric(nolac$Total.seed.weight)~factor(nolac$blockID)+(nolac$Shoot..+nolac$RILcolor*nolac$Early.Growth..cm.+nolac$RILcolor*nolac$Herbiv+nolac$RILcolor*nolac$Cor_Len_Var+nolac$RILcolor*nolac$Sugar_UND+nolac$RILcolor*nolac$Seed.Length+nolac$Nectar_UNC+nolac$RILcolor*nolac$Seed.Mass)/factor(nolac$RILID), na.rm=T)) #RIL ID is highly determinant of seed set, color merely suggestive
#now sig is as follows: Shoot number (1.8e-09),flower colorxsize (0.000511),sugar conc(0.004380), Flower size (0.009491), flower color x seed mass (0.011910), flower color x seed length (0.045310);herbivory(0.059202)

#this appears to give a picture of what is driving seed production - early growth is most important, as shown through shoot numberat the 2 week measuring period
#next, is the interactive effect of flower color and size, with the effect pronounced for large pink flowers
#next are the inidividual effects of flower size and sugar concentration
#lastly are the sig effects of pink flowered plants generally having lighter and smaller seeds than white limbed one (colorx seed traits)
#sugg nature of herbivory indicates that the ones that are less eaten likely produce more seeds, although this effect may be minimal

#let's do a 3-way anova for size and sugar
nest3w<- aov(as.numeric(nolac$Total.seed.weight)~factor(nolac$blockID)+(nolac$RILcolor*nolac$Sugar_UND*nolac$Cor_Len_Var)/factor(nolac$RILID), na.rm=T) #RIL ID is highly determinant of seed set, color merely suggestive
summary(nest3w) #colorx traits are sig, but not colorxsizexsugar or sizexsugar

nestxw<- aov(as.numeric(nolac$Total.seed.weight)~factor(nolac$blockID)+(nolac$Shoot..*nolac$RILcolor*nolac$Early.Growth..cm.*nolac$Herbiv*nolac$Cor_Len_Var*nolac$Sugar_UND*nolac$Seed.Length*nolac$Nectar_UNC*nolac$Seed.Mass)/factor(nolac$RILID), na.rm=T)
summary(nestxw) #too big for f-statistics

nestxw2<- aov(as.numeric(nolac$Total.seed.weight)~factor(nolac$blockID)+(nolac$RILcolor*nolac$Early.Growth..cm.*nolac$Herbiv*nolac$Cor_Len_Var*nolac$Sugar_UND*nolac$Seed.Length*nolac$Nectar_UNC*nolac$Seed.Mass)/factor(nolac$RILID), na.rm=T)
summary(nestxw2) #nothing really coherent, would needto test for multiple hypothesis correction on p-values i think

nest5w<-aov(as.numeric(nolac$Total.seed.weight)~factor(nolac$blockID)+(nolac$RILcolor*nolac$Early.Growth..cm.*nolac$Cor_Len_Var*nolac$Sugar_UND*nolac$Seed.Mass)/factor(nolac$RILID), na.rm=T)
summary(nest5w)#5 most rep categories

nest6w<-aov(as.numeric(nolac$Total.seed.weight)~factor(nolac$blockID)+(nolac$RILcolor*nolac$Early.Growth..cm.*nolac$Cor_Len_Var*nolac$Herbiv*nolac$Sugar_UND*nolac$Seed.Mass)/factor(nolac$RILID), na.rm=T)
summary(nest6w)
#let's do a fitness trait now with no likely predictor outside of line with the big model
summary(aov(as.numeric(nolac$dvec)~factor(nolac$blockID)+(nolac$Shoot..+nolac$RILcolor*nolac$Early.Growth..cm.+nolac$RILcolor*nolac$Herbiv+nolac$RILcolor*nolac$Cor_Len_Var+nolac$RILcolor*nolac$Sugar_UND+nolac$RILcolor*nolac$Seed.Length+nolac$Nectar_UNC+nolac$RILcolor*nolac$Seed.Mass)/factor(nolac$RILID), na.rm=T)) #RIL ID is highly determinant of seed set, color merely suggestive
#line is sig (p=0.018), nectar volume is sig (p=0.0056), everything else is n.s.
summary(aov(as.numeric(nolac$dvec)~factor(nolac$blockID)+(nolac$RILcolor*nolac$Nectar_UNC)/factor(nolac$RILID), na.rm=T)) #RIL ID is highly determinant of seed set, color merely suggestive
#however, a simpler model shows that nearly all of this is due to the effect of line
plot(nolac$dvec, nolac$Nectar_UNC)
####################
#garbage to stop

for (i in 1: nrow(nh)) {
  for (j in 1:ncol(nh)){
    if (is.numeric(nh[i,j])){
      if (isTRUE(nh[i,j]>0) ||isTRUE(nh[i,j]==0)) {
        rec1a[j] <- jdates[j]
      }
    }
  }
  if (isTRUE(sum(rec1a, na.rm=T)>0)){
    fn <-as.numeric(nh[i, ])
    minit <- rec1a[which.max(fn)] #initial mean which may move
    sinit<-(max(rec1a, na.rm=T)-min(rec1a, na.rm=T))/16 #rough estimate of sd
    binit <- min(fn, na.rm=T); #minimizing the data start point
    ainit <- (max(fn, na.rm=T)-min(fn, na.rm=T)) #largest possible range
    if (sum(fn, na.rm=T)==1){
      meanft <- minit
      ftsd <- NA
    } else {
      #setting initial values that are bounded
      # conducting the fit for a line
      Gaussflowerfit <- nls(fn ~ Gaussmodel(rec1a,c(m,sd,a,b)), data.frame(rec1a,fn), start=list(m=minit, sd=sinit, a=ainit, b=binit), trace=TRUE, control=control1)    
      meanft<- summary(Gaussflowerfit)$parameters["m", 1]
      ftsd<-summary(Gaussflowerfit)$parameters["m", 2]
    }
    meanftlist[i]<-meanft
    sdlist[i] <-ftsd
    points(rec1a, as.numeric(nh[i, ]), pch=16, col=(groups[i]+1), xlim=c(200,300))
    lines(dayValues, Gaussmodel(dayValues, coef(Gaussflowerfit)), col=(groups[i]+1), lwd=2)
    rec1a <- rep(NA, ncol(nh))
  }
  print(i)
}

###




#intercepts(flowerPredict)
#quadratic solver adapted from dkMathstats
quadraticRoots <- function(x) { #x needs to be of length 3
  c <- x[1]
  b <- x[2]
  a <- x[3]
  discriminant <- (b^2) - (4*a*c)
  if(discriminant < 0) {
  } else if(discriminant > 0) {
    x_int_plus <- (-b + sqrt(discriminant)) / (2*a)
    x_int_neg <- (-b - sqrt(discriminant)) / (2*a)
    xc <-c(x_int_neg, x_int_plus)
    return(xc)
  }
  else #discriminant = 0  case
    x_int <- (-b) / (2*a)
  return(x_int)
}
qr <-quadraticRoots(flowercurve$coefficients)

#create sequence of day values
dayValues <- seq(18850, 18920, 1)

#18850-18920
#create list of predicted flower number using quadratic model
flowerPredict <- predict(flowercurve,list(mdate1=dayValues, mdate12=dayValues^2))

#create scatterplot of original data values
plot (x=0, xlim=c(18850,18920), ylim=c(0,120))
points(masterdate1, masterdata$Flower...1, pch=16, xlim=c(18850,18920))
#add predicted lines based on quadratic regression model
lines(dayValues, flowerPredict, col=4)
points(qr[1],0)
points(qr[2],0)
qrl <-as.numeric(round(qr[2]))-18849
qrh <-as.numeric(round(qr[1]))-18849

#integrate
flint1 <- flowerPredict[qrl:qrh]
fnumt1 <-sum(flint1)

#end of base code



zerofix1<- c(NA,NA,NA,0, as.character.Date(md1))
zerofix2<-c(NA,NA,NA,0, as.character.Date(md2) )
zerofix<-rbind(zerofix1,zerofix2)
zerofix <-as.data.frame(rbind(zerofix,zerofix))
names(zerofix) <-c('Location', 'RILID', 'RILcolor', 'FlowerNumber', 'FlowerDate')

for (i in 1:nlevels(tempfinal$RILID)) {
  RILhold<-subset(tempfinal, as.numeric(RILID)==i)
  if (is.na(RILhold$FlowerDate[1])){
    next
  } else {
    RILhold <- rbind(RILhold, zerofix)
    RILnn <- as.character(RILhold$RILID[1])
    RILcol <-as.character(RILhold$RILcolor[1])
    if (RILcol=="P"){
      coloring=3
    } else if (RILcol=="W"){
      coloring=2
    } else if (RILcol=="L"){
      coloring=1
    }
    mdate<-as.numeric(RILhold$FlowerDate)
    mdatet <- as.numeric(mdate)^2
    flowercurvei <- lm(RILhold$FlowerNumber~ mdate + mdatet, data=RILhold)
    fcurveint<- flowercurvei$coefficients
    points(x=RILhold$FlowerDate, y=RILhold$FlowerNumber, pch=i%%26, col=as.numeric(RILhold$RILcolor[1]))
    flowerPredict <- predict(flowercurvei,list(mdate=dayValues, mdatet=dayValues^2))
    lines(dayValues, flowerPredict, col=as.numeric(RILhold$RILcolor[1]))
    qr <-quadraticRoots(fcurveint)
    qrl <-as.numeric(round(qr[2]))-18809
    qrh <-as.numeric(round(qr[1]))-18809
    flint <- flowerPredict[qrl:qrh]
    fnumt <-sum(flint)
    tta <-append(RILnn,RILcol)
    ttb <-append(tta, fcurveint)
    ttt <- append(ttb,fnumt)
    fnumatrix[i,] <- ttt
  }
}
fnumatrix <-as.data.frame(fnumatrix)
names(fnumatrix) <- c('RILID', 'RILcolor', 'A1', 'A2', 'A3' ,'FlowerTotal')
legend(18937.55, 125, legend=c("Lac", "White RIL", "Pink RIL"), col=c("red", "blue", "green"), lty=1, cex=0.8)

#calculate midlines for each curve
quadraticMax <- function(b,a) { 
  root<- -b/(2*a)
  return(root)
}

curvedate <-rep(NA, nrow(fnumatrix))
qint1 <- as.numeric(as.character(fnumatrix[,4]))
qint2 <- as.numeric(as.character(fnumatrix[,5]))
for (i in 1:nrow(fnumatrix)){
  curvedate[i] <-quadraticMax(qint1[i], qint2[i])
}

fnumatrix <-cbind(fnumatrix, curvedate)

pvols <-subset(fnumatrix, RILcolor=="P")
ave(as.numeric(as.character(pvols$FlowerTotal)))
wvols<-subset(fnumatrix, RILcolor=="W")
ave(as.numeric(as.character(wvols$FlowerTotal)))
lvols <- subset(fnumatrix, RILcolor=="L")
ave(as.numeric(as.character(lvols$FlowerTotal)))

#testing if counts are inflated because of combination of assoc. btwn/fc and ft and sample bias
taa <-as.numeric(as.character(fnumatrix$curvedate))
tag <-as.numeric(as.character(fnumatrix$FlowerTotal))
plot (x=taa, y=tag,col=(as.numeric(fnumatrix$RILcolor)+1), xlab="Flower Date", ylab="Total Flower Nmber")
legend(18879, 1860, legend=c("Lac", "White RIL", "Pink RIL"), col=c("red", "blue", "green"), lty=1, cex=0.8)


fline<- lm(tag~taa)
summary(fline)
linePredict <- predict(fline,list(taa=dayValues))
lines(dayValues, linePredict)
#abline(v=ave(as.numeric(as.character(pvols$curvedate[1]))), col=3)
abline(v=18895, col=3)
#abline(v=ave(as.numeric(as.character(wvols$curvedate[1]))), col=4)
abline(v=18892, col=4)
#abline(v=ave(as.numeric(as.character(lvols$curvedate[1]))), col=2)
abline(v=18878, col=2)

abline(h=953, col=3, lty=2)
abline(h=765, col=4,lty=2)
abline(h=970, col=2, lty=2)

dxnanova<- anova(fline)

dxcolaov <- aov(formula=curvedate~RILcolor, data=fnumatrix)
nxcolaov <-aov(formula=as.numeric(as.character(FlowerTotal))~RILcolor, data=fnumatrix)

colormax <-rbind(pvols,wvols)

dxcolaov2<- aov(formula=curvedate~RILcolor, data=colormax)
nxcolaov2 <-aov(formula=as.numeric(as.character(FlowerTotal))~RILcolor, data=colormax)


#test plots of data
colors1 = c("pink","grey","brown")
colors2= c("black", "green", "black", "black", "black", "black", "black", "black", "black", "black")
plot(x = masterdata$RILID, y = masterdata$Early.Growth..cm)
plot(x = masterdata$RILcolor, y = masterdata$Early.Growth..cm, xlab="RIL color", ylab= "Size of largest shoot (in)")

plot(x = masterdata$RILID, y = masterdata$Dead., col=colors2, xlab="RIL line ID", ylab="Proportion Alive")
plot(x = masterdata$RILcolor, y = masterdata$Dead., col=colors2, xlab="RIL color", ylab="Proportion Alive")


dddd <-as.data.frame(cbind(as.numeric(mdata1$dvec), as.character(masterdata$RILID)))
dddd <-as.data.frame(rbind(as.numeric(mdata1$dvec), as.character(masterdata$RILID)))
barplot(dddd)
barplot(x = masterdata$RILcolor, y = mdata1$dvec)

plot(x = as.factor(nolac$RILID), y = nolac$Herbiv)
plot(x = as.factor(nolac$RILcolor), y = nolac$Herbiv)

plot(x = nolac$Herbiv, y = nolac$dvec)


plot(x = masterdata$Early.Growth..cm, y = masterdata$Flower...1)
plot(x = masterdata$Early.Growth..cm, y = mdata1$dvec)

lacdata =subset(masterdata, RILcolor=="L")
plot(x = lacdata$Early.Growth..cm, y = lacdata$Flower...1, xlab="Shoot Length (in)", ylab="Observed Flower Number")
lm(lacdata$Flower...1~lacdata$Early.Growth..cm)

plot(x = mdata1$RILID, y = mdata1$Dead.)
plot(x = mdata1$RILcolor, y = mdata1$Dead.)
