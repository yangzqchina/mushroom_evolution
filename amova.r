library(poppr)
setwd("/data4/zqyang/pangenome/gatk_results")
ped<-read.table('124_filter_1k_snp.ped',header = F)
map<-read.table('124_filter_1k_snp.map',header=F)
map[,2]<-paste(map[,1],map[,4],seq='_')
dat<-ped[,7:dim(ped)[2]]
nloc<-dim(dat)[2]/2
nsamples<-dim(dat)[1]
dat2<-as.data.frame(matrix(NA,ncol=nloc,nrow=nsamples))
for(i in seq(nloc)){
  dat2[,i]<-paste0(dat[,2*i-1],dat[,2*i])
}
colnames(dat2)<-map[,2]
rownames(dat2)<-ped[,2]
pop<-read.csv('124_pop.csv')
snp<-df2genind(dat2[pop$sample,],ploidy = 2,ncode=1)
pop$'Pop_Subpop'<-paste(pop$Pop,pop$Subpop,sep='_')
pop$Subpop<-as.factor(pop$Subpop)
other(snp)$population_hierarchy<-pop[,c('Pop_Subpop','Pop','Subpop')]
snp$pop<-pop$Pop
strata(snp) <- other(snp)$population_hierarchy[-1]

sink('/data4/zqyang/pangenome/gatk_results/amova_results.txt')
agc <- as.genclone(snp)
agc
amova.result <- poppr.amova(agc, ~Pop/Subpop)
amova.result
amova.test <- randtest(amova.result) # Test for significance
amova.test
sink()
pdf('/data4/zqyang/pangenome/gatk_results/amova_results.pdf',12,10)
plot(amova.test)
dev.off()
results<-cbind(amova.result$results,amova.result$componentsofcovariance)
write.csv(results,'/data4/zqyang/pangenome/gatk_results/amova_results.csv')
