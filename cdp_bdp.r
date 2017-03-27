args<-commandArgs(T)
dat<-read.table(args[1],header = T)
tags<-5:dim(dat)[2]
dat_values<-dat[,5:dim(dat)[2]]
dat_means<-colMeans(dat_values,na.rm = T)
for(i in seq(length(dat_means))){
  dat_values[,i]<-dat_values[,1]/dat_means[i]
}
dat[,5:dim(dat)[2]]<-dat_values
dat<-na.omit(dat)
bdp<-dat
bdp_values<-bdp[,tags]
bdp_values[bdp_values>0.1]=1
bdp_values[bdp_values<=0.1]=0
bdp[,tags]=bdp_values
bdp_mean=rowMeans(bdp_values,na.rm=T)
#select_line<-na.omit(bdp_mean)
select_line<-bdp_mean>0.05 & bdp_mean<0.95
bdp<-bdp[select_line,]
cdp<-dat[select_line,]
bdp[,5:dim(dat)[2]][dat[,5:dim(dat)[2]]==1]=2
bdp[,5:dim(dat)[2]][dat[,5:dim(dat)[2]]==0]=1
write.table(cdp,args[2],row.names = F,quote=F)
write.table(bdp,args[3],row.names = F,quote=F)
