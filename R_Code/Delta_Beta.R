S_Males = 30396
S_Females = 30397
for(i in 5:50)
{
	print(paste("HD",i,sep=""))
	temp=read.table(paste("/scratch/joel/broad/ImputeDosages/HD",i,"/ibd_rel5_EU_clean_HD",i,"_HG19.core.filtered.dosage",sep=""),header=T,sep="\t")
	Alleles = read.table(paste("/scratch/joel/broad/ImputeDosages/HD",i,"/ibd_rel5_EU_clean_HD",i,"_HG19.core.allele",sep=""),header=F)
	fam = read.table("/scratch/ming/IBD.fam",header=F)
	temp = cbind(sex=fam[,5],temp)
	Males =  temp[temp[,1] == 1,]
	Females =  temp[temp[,1] == 2,]
	Alleles[,1] = gsub(":", ".", Alleles[,1], fixed = T)
	inter = intersect(colnames(Males),as.character(Alleles[,1]))
	subset_SNPs =  Alleles[as.character(Alleles[,1]) %in% inter,]

	x = t(Males[,2:ncol(Males)])
	x = x[2:nrow(x),] 
	Males = cbind(subset_SNPs[,3:4],x)

	y = t(Females[,2:ncol(Females)])
	y = y[2:nrow(y),] 
	Females = cbind(subset_SNPs[,3:4],y)


	Males_fam = fam[fam[,5] == 1,]
	Females_fam = fam[fam[,5] == 2,]
	
	rownames(Males) = as.character(subset_SNPs[,1])
	rownames(Females) = as.character(subset_SNPs[,1])

	write.table(Males_fam,col.names=F,row.names=F,quote=F,file=paste("/scratch/mahmoud/test/HD_Regions/HD",i,"_Males.fam",sep=""))
	write.table(Females_fam,col.names=F,row.names=F,quote=F,file=paste("/scratch/mahmoud/test/HD_Regions/HD",i,"_Females.fam",sep=""))
	write.table(Males,col.names=F,quote=F,file=paste("/scratch/mahmoud/test/HD_Regions/HD",i,"_Males.dosage",sep=""))
	write.table(Females,col.names=F,quote=F,file=paste("/scratch/mahmoud/test/HD_Regions/HD",i,"_Females.dosage",sep=""))
	
	mainDir = "/scratch/mahmoud/test"
	subDir = paste("HD",i,"_simulation",sep="")
	dir.create(file.path(mainDir, subDir))
	
	temp = t(temp)
	for(j in 1:100)
	{
		print(j)
		shuffle = sample(fam[,5],replace=F)
		temp[1,] = shuffle
		temp = temp[,order(temp[1,])]
		
		Males =  temp[2:nrow(temp),1:S_Males]
		Females =  temp[2:nrow(temp),S_Females:ncol(temp)]
		
		Males_fam = fam[as.character(fam[,1]) %in% as.character(Males[1,]),]
		Females_fam = fam[as.character(fam[,1]) %in% as.character(Females[1,]),]
		
		Males = Males[-1,]
		Females = Females[-1,]
		
		Males = cbind(subset_SNPs[,3:4],Males)
		Females = cbind(subset_SNPs[,3:4],Females)
		rownames(Males) = as.character(subset_SNPs[,1])
		rownames(Females) = as.character(subset_SNPs[,1])
		
		write.table(as.matrix(Males_fam),col.names=F,row.names=F,quote=F,file=paste(file.path(mainDir, subDir),"/HD",i,"_s",j,"_Males.fam",sep=""));
		write.table(as.matrix(Females_fam),col.names=F,row.names=F,quote=F,file=paste(file.path(mainDir, subDir),"/HD",i,"_s",j,"_Females.fam",sep=""));
		write.table(as.matrix(Males),col.names=F,quote=F,file=paste(file.path(mainDir, subDir),"/HD",i,"_s",j,"_Males.dosage",sep=""));
		write.table(as.matrix(Females),col.names=F,quote=F,file=paste(file.path(mainDir, subDir),"/HD",i,"_s",j,"_Females.dosage",sep=""))
	}
}





S_Males = 30396
S_Females = 30397
for(i in 1:184)
{
	print(paste("HD",i,sep=""))
	temp=read.table(paste("/scratch/joel/broad/ImputeDosages/HD",i,"/ibd_rel5_EU_clean_HD",i,"_HG19.core.filtered.dosage",sep=""),header=T,sep="\t",check.names=F)
	Alleles = read.table(paste("/scratch/joel/broad/ImputeDosages/HD",i,"/ibd_rel5_EU_clean_HD",i,"_HG19.core.allele",sep=""),header=F)
	fam = read.table("/scratch/ming/IBD.fam",header=F)
	temp = cbind(sex=fam[,5],temp)
	
	Alleles[,1] = gsub(":", ".", Alleles[,1], fixed = T)
	colnames(temp = gsub(":", ".", colnames(temp), fixed = T)
	inter = intersect(colnames(temp),as.character(Alleles[,1]))
	subset_SNPs =  Alleles[as.character(Alleles[,1]) %in% inter,]
	
	temp = t(temp)
	temp = temp[,order(temp[1,])]
		
	Males =  temp[2:nrow(temp),1:S_Males]
	Females =  temp[2:nrow(temp),S_Females:ncol(temp)]
	
	Males_fam = fam[as.character(fam[,1]) %in% as.character(Males[1,]),]
	Females_fam = fam[as.character(fam[,1]) %in% as.character(Females[1,]),]
	Males_fam[,5] = rep(1,nrow(Males_fam))
	Females_fam[,5] = rep(2,nrow(Females_fam))
	
	Males = Males[-1,]
	Females = Females[-1,]
	
	Males = cbind(subset_SNPs[,3:4],Males)
	Females = cbind(subset_SNPs[,3:4],Females)
	rownames(Males) = as.character(subset_SNPs[,1])
	rownames(Females) = as.character(subset_SNPs[,1])		

	write.table(as.matrix(Males_fam),col.names=F,row.names=F,quote=F,file=paste("/scratch/mahmoud/test/HD_Regions/HD",i,"_Males.fam",sep=""))
	write.table(as.matrix(Females_fam),col.names=F,row.names=F,quote=F,file=paste("/scratch/mahmoud/test/HD_Regions/HD",i,"_Females.fam",sep=""))
	write.table(as.matrix(Males),col.names=F,quote=F,file=paste("/scratch/mahmoud/test/HD_Regions/HD",i,"_Males.dosage",sep=""))
	write.table(as.matrix(Females),col.names=F,quote=F,file=paste("/scratch/mahmoud/test/HD_Regions/HD",i,"_Females.dosage",sep=""))
}














S_Males = 30396
S_Females = 30397
for(i in 51:100)
{
	print(paste("HD",i,sep=""))
	Males = read.table(paste("/scratch/mahmoud/test/HD_Regions/HD",i,"_Males.dosage",sep=""),header=F,sep=" ")
	Females = read.table(paste("/scratch/mahmoud/test/HD_Regions/HD",i,"_Females.dosage",sep=""),header=F,sep=" ")
	
	Males_fam = read.table(paste("/scratch/mahmoud/test/HD_Regions/HD",i,"_Males.fam",sep=""),header=F)
	Females_fam = read.table(paste("/scratch/mahmoud/test/HD_Regions/HD",i,"_Females.fam",sep=""),header=F)
	temp_fam = rbind(Males_fam,Females_fam)
	
	temp=cbind(Males,Females[4:ncol(Females)])
		
	mainDir = "/scratch/mahmoud/test"
	subDir = paste("HD",i,"_simulation",sep="")
	dir.create(file.path(mainDir, subDir))
	shuffle = lapply(1:100,function(obj) {sample(4:(68427+3),replace=F)})
	
	Males_shuffled = list()
	Males_shuffled_fam = list()
	Females_shuffled = list()
	Females_shuffled_fam = list()
	for(j in 1:100)
	{
		print(j)
		shuffle_M = shuffle[[j]][1:S_Males]
		Males_shuffled[[j]] = temp[,c(1:3,shuffle_M)]
		Males_shuffled_fam[[j]] = temp_fam[(shuffle_M-3),]
		
		shuffle_F = shuffle[[j]][S_Females:68427]
		Females_shuffled[[j]] = temp[,c(1:3,shuffle_F)]
		Females_shuffled_fam[[j]] = temp_fam[(shuffle_F-3),]
		
		Males_shuffled_fam[[j]][,5] = rep(1,nrow(Males_shuffled_fam[[j]]))
		Females_shuffled_fam[[j]][,5] = rep(2,nrow(Females_shuffled_fam[[j]]))
	}
	counter <<- 1
	lapply(Males_shuffled,function(obj){system.time(write.table(as.matrix(obj),col.names=F,row.names=F,quote=F,file=paste(file.path(mainDir, subDir),"/HD",i,"_s",counter,"_Males.fam",sep=""))); counter<<-counter+1})

	
}


		system(paste("awk '{print $1,$2,$3,",shuffle_M,"}' /scratch/mahmoud/test/HD_test/HDtest.dosage > /scratch/mahmoud/test/HD_test/HDtest_",j,".dosage",sep=""),intern=TRUE)



