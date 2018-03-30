Regions = read.table("Regions_details.txt",header=T,stringsAsFactor=F)
x=c("Belgium_UC_1_SNPs_v2.txt","Belgium_UC_2_SNPs_v2.txt","Belgium_UC_3_SNPs_v2.txt","Belgium_CD_1_SNPs_v2.txt","Belgium_IBD_1_SNPs_v2.txt")
temp=lapply(x,read.table,header=F,stringsAsFactor=F)
n=colnames(read.table("Belgium_CD_2_SNPs_v2.txt",header=T))
for(i in 1:length(temp))
{
	for(j in 1:nrow(temp[[i]]))
	{
		test = Regions[Regions$CHR == temp[[i]][j,1] & (Regions$start*1000000) <= temp[[i]][j,3] & (Regions$end*1000000) >= temp[[i]][j,3],]
		if(nrow(test) == 1)
		{		
			temp[[i]][j,11] = test$Genes
			temp[[i]][j,12] = test$HD_name
		} else {
			print("NOT FOUND")
		}
	}
	colnames(temp[[i]]) = n
	write.table(temp[[i]],file=x[i],quote=F,row.names=F)
}
