library('stringr')
library('qvalue')


for(j in 1:3)
{
	Disease_name = c("CD","UC","IBD")
	regexp = "[[:digit:]]+"
	if(Disease_name[j] == "IBD"){
		f = file.path("/scratch/new_imputation/probABEL_output",paste("HD",1:184,"_probABEL_add.out.txt",sep=""))
	}else{
		f = file.path("/scratch/new_imputation/probABEL_output",paste(Disease_name[j],"_HD",1:184,"_probABEL_add.out.txt",sep=""))
	}
	assoc.logistic <- lapply(f[file.exists(f)], read.table,header=T,stringsAsFactors=F,sep="")

	HD_regions = as.numeric(unlist(lapply(1:length(f[file.exists(f)]),function(i){ str_extract(strsplit(f[file.exists(f)],"/")[[i]][5],regexp) })))

	f = file.path("/scratch/impute2Apr/impute2",paste("HD",HD_regions,"/ibd_rel5_EU_clean_HD",HD_regions, "_HG19.core.allele",sep=""))
	HD_positions =  lapply(f[file.exists(f)], read.table,header=F,stringsAsFactors=F,sep="")

	for(i in 1:length(assoc.logistic))
	{
		print(paste("HD ",HD_regions[i],sep=""))
		#png(paste('/scratch/new_imputation/probABEL_plots_v2/',Disease_name[j],'_mht_HD',HD_regions[i],'.png',sep=''),width=290,height=210,units="mm",res=450)
		T=(assoc.logistic[[i]]$beta_SNP_sex / assoc.logistic[[i]]$sebeta_SNP_sex)^2
		p=pchisq(T,df=1,lower=F)
		assoc.logistic[[i]]$pvalue = p
		assoc.logistic[[i]]$pos = unlist(HD_positions[[i]][2])
		temp = assoc.logistic[[i]]
		temp = temp[temp$Freq1 >= 0.05,]
		#plot(temp$pos,-log10(temp$pvalue),xlab="BP position",main=paste("HD ",HD_regions[i],sep=""))
		#points(x=head(temp[order(temp$pvalue),],sum(sort(qvalue(na.omit(temp$pvalue))$qvalues) <= 0.1))$pos,y=-log10(head(temp[order(temp$pvalue),],sum(sort(qvalue(na.omit(temp$pvalue))$qvalues) <= 0.1))$pvalue),col="Red",pch = 19)
		#dev.off()
	
		png(paste('/scratch/new_imputation/probABEL_plots_v2/',Disease_name[j],'_qq_HD',HD_regions[i],'.png',sep=''),width=290,height=210,units="mm",res=450)
		qqunif(temp$pvalue,main=paste("HD ",HD_regions[i],sep=""))
		dev.off()
	}
	rm(list=ls())
}

system(paste("rm /scratch/new_imputation/probABEL_output/",Disease_name,"_qvalues_v3.txt",sep="" ))
system(paste("touch /scratch/new_imputation/probABEL_output/",Disease_name,"_qvalues_v3.txt",sep="" ))
sink(paste("/scratch/new_imputation/probABEL_output/",Disease_name,"_qvalues_v3.txt",sep=""))
for(i in 1:length(assoc.logistic))
{
	print(paste("HD ",HD_regions[i],sep=""))
	T=(assoc.logistic[[i]]$beta_SNP_sex / assoc.logistic[[i]]$sebeta_SNP_sex)^2
	p=pchisq(T,df=1,lower=F)
	assoc.logistic[[i]]$pvalue = p
	temp = assoc.logistic[[i]]
	temp = temp[temp$Freq1 >= 0.05,]
	summary(qvalue(na.omit(temp$pvalue)))
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

Disease_name = c("CD","UC","IBD")
for(i in 1:3)
{
	regexp = "[[:digit:]]+"
	if(Disease_name[i] == "IBD"){
		f = file.path("/scratch/new_imputation/probABEL_output",paste("HD",1:184,"_probABEL_add.out.txt",sep=""))
	}else{
		f = file.path("/scratch/new_imputation/probABEL_output",paste(Disease_name[i],"_HD",1:184,"_probABEL_add.out.txt",sep=""))
	}
	HD_regions = as.numeric(unlist(lapply(1:length(f[file.exists(f)]),function(i){ str_extract(strsplit(f[file.exists(f)],"/")[[i]][5],regexp) })))
	system(paste("rm temp.txt",sep="" ))
	system(paste("less /scratch/new_imputation/probABEL_output/",Disease_name[i],"_qvalues_v3.txt | grep q-value > temp.txt",sep="" ))
	qv = read.table("/scratch/new_imputation/probABEL_output/temp.txt",sep="",stringsAsFactors=F)
	qv[,1] = HD_regions
	qv = qv[order(qv[,7]),]
	qv = qv[qv[,7] >= 1,]
	write.table(qv,file=paste("/scratch/new_imputation/probABEL_sign/",Disease_name[i],"/",Disease_name[i],"_significant_v3.txt",sep=""),row.names=F,quote=F)
}



####################################################################################################
system("rm  /scratch/temp_file.txt")
system(paste("grep rs185970107 /scratch/impute2Apr/impute2/HD19/ibd_rel5_EU_clean_HD19_HG19.core.gen > /scratch/temp_file.txt",sep=""))
temp = read.table('/scratch/temp_file.txt',stringsAsFactors=F,sep=" ")
#Genotypes = data.frame(matrix(0,ncol=3,nrow=(ncol(temp)-4)/2))
count <<- 1
insert = function(i)
{
		print(count)
		output = matrix(0,ncol=3,nrow=1)
		output[1,1] <- temp[,i]
		output[1,2] <- temp[,i+1]
		output[1,3] <- 1 - (output[1,1] + output[1,2])
		count <<- count + 1
		output
		
		
}
Genotypes = lapply(seq(from=5,to=ncol(temp),by=2),insert)
countAA=0
countAB=0
countBB=0
count_genotypes = function(i)
{
	#print(i)
	#if(i %in% which(is.na(fam[,2]))) {
	#	print("skip")
	#} else {
		AA = as.character(which.max(Genotypes[[i]]))
		switch(AA, "1"={countAA<<-countAA+1},"2"={countAB<<-countAB+1},"3"={countBB<<-countBB+1})
		if(max(Genotypes[[i]]) < 1) {
			print(max(Genotypes[[i]]))
		}
	#}
}
lapply(1:length(Genotypes),count_genotypes)

probA = ((2*countAA) + countAB) / (2*(countAA+countAB+countBB))
probB = 1 - probA


HD_SNPs = read.table('/scratch/impute2Apr/impute2/HD19/ibd_rel5_EU_clean_HD19_HG19.core.allele',stringsAsFactors=F)
x = as.numeric(system(paste("grep -nr ",HD_SNPs[i,1]," /scratch/impute2Apr/impute2/HD19/ibd_rel5_EU_clean_HD19_HG19.core.allele | gawk '{print $1}' FS=':'",sep=""),intern = TRUE)) + 2
systme("rm /scratch/tmp.mldose")
y = system(paste("cut -f",x," -d$'\t'  /scratch/new_imputation/mldose_files/HD19.mldose > /scratch/tmp.mldose",sep=""),intern = TRUE)

${PLINK} --noweb --dosage ${DOSAGE_DIR}${DOSAGE_FILE} dose1 format=2 noheader skip1=1



######################## ANalysis using glm #############################


HD_SNPs = read.table('/scratch/impute2Apr/impute2/HD19/ibd_rel5_EU_clean_HD19_HG19.core.allele',stringsAsFactors=F)
#P_interaction = matrix(0,nrow=nrow(HD_SNPs),ncol=2)

analyize_HD_region = function(i)
{
	print(i)
	output = matrix(0,ncol=2,nrow=1)
	x = as.numeric(system(paste("grep -nr ",HD_SNPs[i,1]," /scratch/impute2Apr/impute2/HD19/ibd_rel5_EU_clean_HD19_HG19.core.allele | gawk '{print $1}' FS=':'",sep=""),intern = TRUE)) + 2
	system("rm /scratch/tmp.mldose")
	y = system(paste("cut -f",x," -d$'\t'  /scratch/new_imputation/mldose_files/HD19.mldose > /scratch/tmp.mldose",sep=""),intern = TRUE)
	dose = read.table('/scratch/tmp.mldose')
	res = summary(glm(fam[,2] ~ fam[,3] + fam$PCA1 + fam$PCA2 + fam$PCA3 + fam$PCA4 + fam$PCA5 + dose[,1] + fam[,3]*dose[,1]),family=logit)
	output[1,1] = HD_SNPs[i,1]
	output[1,2] = res$coefficients[9,4]
	output
}
print(paste("total SNPs are ",nrow(HD_SNPs),sep=""))
P_interaction = lapply(1:3,analyize_HD_region)


