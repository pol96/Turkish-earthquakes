earthquake=read.csv('export_df.csv')
earthquake=earthquake[earthquake$country=='turkey',]
earthquake=earthquake[,c(18,4:14,15,16,17)]
View(earthquake)
eqk=earthquake[,10:14]
attach(eqk)
View(eqk)
#####location indexes####
mean_data=colMeans(eqk,na.rm = T)
mean_data
median_data=apply(eqk, 2, median)
median_data
gmode=function(x){
  u=unique(x)
  u[which.max(tabulate(match(x,u)))]
}
mode_data=apply(eqk,2,gmode)
mode_data

LocInd=rbind(mean_data,median_data,mode_data)
colnames(LocInd)=colnames(eqk)
row.names(LocInd)=c('mean','median','mode')
View(LocInd)

#####dispersion indexes####
sd_data=apply( 2, sd)
sd_data
var_data=var(eqk)
var_data
range_data=apply(eqk,2,range)
range_data

cor_data=cor(eqk)
library(corrplot)
library(RColorBrewer)
  #1. correlation plot
  x11()
  corrplot(cor_data,'numbe',type = 'upper',col = brewer.pal(n=8, name='RdYlBu'))
  graphics.off()
  boxplot(eqk[,c(1,2,4,5)])
  
#####shape indexes####
library(fBasics)
sk_data=apply(eqk,2,skewness) #symmetry of the distribution (0--> perfect symmetry)
sk_data  
kr_data=apply(eqk,2,kurtosis)  #measure of 'flatness' (higher value corresponds to closer data)
kr_data    

#####Depth####
  #1. depth
library(depth)
library(circular)
med(earthquake[,c(2,3)],method = 'Liu')
  #2. plot hexbin
library(hexbin)
b=hexbin(x = xm,y = md,xbins = 20)
x11()
plot(b)
  #3. isodepth
x11()
bagplot(earthquake[,c(3,2)])# contour plot of Tukey's depth function
isodepth(earthquake[,c(2,3)])

#####histograms####
x11()
par(mfrow=c(2,3))
hist(depth, col='light green', main='depth')
hist(mb, col='light blue', main='mb index')
hist(md,col='light blue', main='md index')
hist(ms,col='light blue', main='ms index')
hist(xm,col='light blue', main='mx index')
x11()
hist(mw, col = 'blue', main='Moment magnitude Scale (mw)',ylim=c(0,7000))
boxplot(richter)
#####Richter cluster####
  #3 classes:
      #from 0 to 3.9 || code: 1
      #from 4.0 to 5.9 || code: 2
      #from 6.0 || code: 3

clustering=function(x=richter){
  s=NULL
  for(i in 1:length(x)){
    if(x[i]<4){
      s[i]='Minor'
    } else if (x[i]<6 & x[i]>=4){
      s[i]='Moderate'
    } else{
      s[i]='Major'
    }
    s=c(s,s[i])
  }
  return(s)
}
richter_cluster=clustering(richter)
eqk=cbind(eqk,richter_cluster[1:11850])
earthquake=cbind(earthquake,richter_cluster[1:11850])
View(eqk)



#####Year cluster####
  #1. absolute
e=tapply(earthquake$year, list(earthquake$`richter_cluster[1:11850]`,earthquake$year), length)
et=as.data.frame(t(e))
et=as.data.frame(cbind(as.numeric(rownames(et)),et))
et[is.na(et)]=0
et[sort(et$`as.numeric(rownames(et))`),]
yr_vec=c(seq(1910,2017,10),2017)
xet=matrix(0,length(yr_vec),4)
xet[,1]=yr_vec
for (i in 1:nrow(et)) {
  for(j in 2:nrow(xet)){
    if(et[i,1]<=xet[j,1] & et[i,1]>xet[j-1,1]){
      for(z in 2:ncol(xet)){
        xet[j,z]=xet[j,z]+et[i,z]
      }
    }
  }
}
sum(xet[,2:4])

txet=t(xet)
colnames(txet)=txet[1,]
txet=txet[2:4,]
txet
txet=rbind(txet[2,],txet[3,],txet[1,])
row.names(txet)=c('Minor','Moderate','Major')
r=rbind(txet[3,],txet[2,],txet[1,])
row.names(r)=c('Major','Moderate','Minor')
r=r[,2:11]
r
View(txet)
txet=txet[,2:11]
  #2. relative
myxet=t(r)
c=rowSums(myxet)
for(i in 1:nrow(myxet)){
  for(j in 1:ncol(myxet)){
    myxet[i,j]=myxet[i,j]/c[i]
  }
}
myxet
#####Barplots (absolute and relative)####
  #1. absolute
library(RColorBrewer)
x11()
a=barplot(
  r,
  ylim = c(0,3500),
  main = 'Occurred earthquakes by decade (Absolute)',
  col=brewer.pal(nrow(txet), "Set1"),
  ylab="Total Earthquakes",
  xlab='Years',
  legend.text=TRUE,
  args.legend = list(x='top',bty='n',inset=c(0,-0.065),ncol=3)
)
b=colSums(r)
text(a,b+300,b)
  #2. relative

x11()
barplot(
  t(myxet),
  main = 'Occurred earthquakes by decade (Relative)',
  col=brewer.pal(nrow(txet), "Set1"),
  ylab="Total Earthquakes per Decade",
  xlab='Years',
  legend.text=TRUE,
  args.legend = list(x='top',bty='n',inset=c(0,-0.065),ncol=3)
)
#####Barplots on Moderate and Major####
  #1. Absolute
r=txet[c(3,2),]
r=r[,2:11]
x11()
a=barplot(
  r,
  ylim=c(0,700),
  main = 'Occurred earthquakes by decade (Moderate and Major)',
  col=brewer.pal(nrow(txet), "Set1"),
  ylab="Total Earthquakes",
  xlab='Years',
  legend.text=TRUE,
  args.legend = list(x='top',bty='n',inset=c(0,-0.065),ncol=2)
)
b=colSums(r)
text(a,b+50,b)

  #2. Relative
#####
rel=t(txet[c(2,3),])
c=rowSums(rel)
for(i in 1:nrow(rel)){
  for(j in 1:ncol(rel)){
    rel[i,j]=rel[i,j]/c[i]
  }
}
rel
#####
x11()
barplot(
  t(rel),
  main = 'Occurred earthquakes by decade (Relative)',
  col=brewer.pal(nrow(txet), "Paired"),
  ylab="Total Earthquakes per Decade",
  xlab='Years',
  legend.text=TRUE,
  args.legend = list(x='top',bty='n',inset=c(0,-0.065),ncol=3)
)

#####PCA####
eqk_scale=as.data.frame(scale(eqk))
pc=princomp(eqk_scale, scores = T)
  #1. loadings
l=pc$loadings
x11()
par(mfrow=c(3,2))
for(i in 1:5) barplot(l[,i], main = i)
  #2. variance
x11()
layout(matrix(c(2,3,1,3),2,byrow=T))
barplot(pc$sdev^2, las=2, main='Principal components', ylim=c(0,4), ylab='Variances') 
abline(h=1, col='blue')
barplot(sapply(eqk_scale,sd)^2, las=2, main='Original variables', ylim=c(0,4), ylab='Variances') 
plot(cumsum(pc$sd^2)/sum(pc$sd^2), type='b', axes=F, xlab='number of components', ylab='contribution to the total variace', ylim=c(0,1)) 
abline(h=1, col='blue')
abline(h=0.8, lty=2, col='blue')
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:ncol(s_st),labels=1:ncol(s_st),las=2)
