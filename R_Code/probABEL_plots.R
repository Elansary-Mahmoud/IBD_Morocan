library('stringr')
library('qvalue')
regexp = "[[:digit:]]+"

Disease_name = "CD"
if(Disease_name == "IBD"){
	f = file.path("/scratch/new_imputation/probABEL_output",paste("HD",1:184,"_probABEL_add.out.txt",sep=""))
}else{
	f = file.path("/scratch/new_imputation/probABEL_output",paste(Disease_name,"_HD",1:184,"_probABEL_add.out.txt",sep=""))
}
assoc.logistic <- lapply(f[file.exists(f)], read.table,header=T,stringsAsFactors=F,sep="")

HD_regions = as.numeric(unlist(lapply(1:length(f[file.exists(f)]),function(i){ str_extract(strsplit(f[file.exists(f)],"/")[[i]][5],regexp) })))

f = file.path("/scratch/impute2Apr/impute2",paste("HD",HD_regions,"/ibd_rel5_EU_clean_HD",HD_regions, "_HG19.core.allele",sep=""))
HD_positions =  lapply(f[file.exists(f)], read.table,header=F,stringsAsFactors=F,sep="")



for(i in 1:length(assoc.logistic))
{
	print(paste("HD ",HD_regions[i],sep=""))
	png(paste('/scratch/new_imputation/probABEL_plots/',Disease_name,'_mht_HD',HD_regions[i],'.png',sep=''),width=290,height=210,units="mm",res=450)
	T=(assoc.logistic[[i]]$beta_SNP_sex / assoc.logistic[[i]]$sebeta_SNP_sex)^2
	p=pchisq(T,df=1,lower=F)
	assoc.logistic[[i]]$pvalue = p
	assoc.logistic[[i]]$pos = unlist(HD_positions[[i]][2])
	temp = assoc.logistic[[i]]
	#temp = temp[temp$Freq1 >= 0.05,]
	plot(temp$pos,-log10(temp$pvalue),xlab="BP position",main=paste("HD ",HD_regions[i],sep=""))
	dev.off()
	
	png(paste('/scratch/new_imputation/probABEL_plots/',Disease_name,'_qq_HD',HD_regions[i],'.png',sep=''),width=290,height=210,units="mm",res=450)
	qqunif(temp$pvalue,main=paste("HD ",HD_regions[i],sep=""))
	dev.off()
}


system(paste("touch /scratch/new_imputation/probABEL_output/",Disease_name,"_qvalues.txt",sep="" ))
sink(paste("/scratch/new_imputation/probABEL_output/",Disease_name,"_qvalues.txt",sep=""))
for(i in 1:length(assoc.logistic))
{
	print(paste("HD ",HD_regions[i],sep=""))
	T=(assoc.logistic[[i]]$beta_SNP_sex / assoc.logistic[[i]]$sebeta_SNP_sex)^2
	p=pchisq(T,df=1,lower=F)
	summary(qvalue(na.omit(p)))
}
sink()

for(i in 1:length(assoc.logistic))
{
	print(paste("HD ",HD_regions[i],sep=""))
	T=(assoc.logistic[[i]]$beta_SNP_sex / assoc.logistic[[i]]$sebeta_SNP_sex)^2
	p=pchisq(T,df=1,lower=F)
	assoc.logistic[[i]]$pvalue = p
	assoc.logistic[[i]]$pos = unlist(HD_positions[[i]][2])
	
}





plink --noweb --dosage \ ./impute2/HD1/ibd_rel5_EU_clean_HD1_HG19.core.gen format=2 noheader \ skip1=1 --fam IBD.fam --logistic --covar
pcs.covar.v2.txt \ --covar-number 1-5 --out ibd_rel5_EU_clean_HD1_HG19.core.IBD \ --allow-no-sex



