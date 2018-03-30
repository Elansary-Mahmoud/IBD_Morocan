args <- commandArgs(trailingOnly = TRUE)
Disease_name = as.character(args[1])
country = as.character(args[2])
Morocan_Data = read.table("Fichier_Phenotypes_final_bis.txt")
SAMPLES = read.table("AllImmunochipsMarocains-29102012-build37.fam")
rownames(Morocan_Data) = as.character(Morocan_Data[,2])
Morocan_Data=Morocan_Data[as.character(SAMPLES[,2]),]

Disease_code = as.numeric(switch(Disease_name, "IBD"="3", "CD"="4", "UC"="5") )
summary(as.factor(Morocan_Data[,Disease_code]))
Broad = read.table(paste(Disease_name,".fam",sep=""))
fam= Broad
Samples_origin = read.csv("../IBD_Morocan/Samples.csv",sep=",",header=T)
D = switch(Disease_name, "IBD"="Indeterminate", "CD"="Crohn's Disease", "UC"="Ulcerative Colitis")
if(Disease_name == "IBD") {
	Broad = Samples_origin[Samples_origin$country_origin == country & (Samples_origin$diag == "Crohn's Disease" |  Samples_origin$diag == "Ulcerative Colitis" | Samples_origin$diag == "Unaffected") ,]
} else {
	Broad = Samples_origin[Samples_origin$country_origin == country & (Samples_origin$diag == D | Samples_origin$diag == "Unaffected") ,]
}
Broad = Broad[as.character(Broad$sample_id) %in% fam[,1],]
regions_b37 = read.table("justin_UCSC.txt",stringsAsFactor=F)
positions = read.table('ibd_rel5_EU_clean.v3.HG19.updateID.bim',stringsAsFactor=F)
Negative_SNP = read.table('Negative_SNPs_list.txt',header=T)
plink = "/home/ulg/genan/mahmoud/plink-1.07-x86_64/plink"
Genofile = "/home/ulg/genan/mahmoud/Immunochip/ibd_rel5_EU_clean.v3.HG19.updateID" 
outfile = "/home/ulg/genan/mahmoud/Immunochip/Results/"
gec = "/home/ulg/genan/mahmoud/GEC/gec/gec.jar"
Target = read.table('all_SNP_region_frequency.frq',header=T,stringsAsFactor=F)
Target = Target[Target$SNP %in% setdiff(Target$SNP,Negative_SNP[,1]),]
Leila_SNP=read.table("../IBD_Morocan/Variants_IC.txt",header=T,sep="\t")
Leila_SNP = Leila_SNP[as.character(Leila_SNP$IC) %in% as.character(Target$SNP),]
risk_regions = read.table('Regions.txt',header=T) #'risk_regions_justin.txt'

for(S in 1:3)
{

	unaff = Broad[Broad[,6] == "Unaffected",]
	if(Disease_name="IBD") {
		aff = Broad[Broad[,6] == "Crohn's Disease" | Broad[,6] == "Ulcerative Colitis",]
	} else {
		aff = Broad[Broad[,6] == D,]
	}
	
	

	aff = aff[sample(x=1:nrow(aff),size=summary(as.factor(Morocan_Data[,Disease_code]))["2"]),]
	unaff = unaff[sample(x=1:nrow(unaff),size=summary(as.factor(Morocan_Data[,Disease_code]))["1"]),]
	random_samples = c(as.character(aff[,2]),as.character(unaff[,2]))
	write.table(cbind(random_samples,random_samples),file=paste("Keep_",country,"_",Disease_name,"_",S,".txt",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
	
	system(paste(plink," --noweb --bed ",Genofile,".bed --bim ",Genofile,".bim --fam ",Disease_name,".fam --keep Keep_",country,"_",Disease_name,"_",S,".txt --recode --make-bed --out ",outfile,country,"_",Disease_name,"_",S,sep=""))
	New_Genofile = paste("/home/ulg/genan/mahmoud/Immunochip/Results/",country,"_",Disease_name,"_",S,sep="")
	
	system(paste(plink," --noweb --bfile ",New_Genofile," --freq --outx ",outfile,country,"_",Disease_name,"_",S,"_Broad_GWAS_SNP_frequency",sep=""))
	
	write.table(fam[fam[,1] %in% random_samples,c(1,2,6)],file=paste("pheno_",country,"_",Disease_name,"_",S,".txt",sep=""),quote=F,col.names=F,row.names=F)
	system(paste(plink," --noweb --bfile ",New_Genofile," --covar pcs.covar.v2.txt --covar-number 1-5 --pheno ","pheno_",country,"_",Disease_name,"_",S,".txt"," --logistic --out ",outfile,"logistic_",country,"_",Disease_name,"_",S,sep=""))
	
	GWAS_frq = read.table(paste(outfile,country,"_",Disease_name,"_",S,"_Broad_GWAS_SNP_frequency.frq",sep=""),header=T)
	GWAS_frq$pos = positions[,4]

	Target_frq = GWAS_frq[as.character(GWAS_frq$SNP) %in% as.character(Leila_SNP$dbSNP),]

	################################################################################################
	####################	Remove The negative control regions	################################

	temp = intersect(as.character(GWAS_frq$SNP),as.character(Negative_SNP[,1]))
	GWAS_frq = GWAS_frq[GWAS_frq$SNP %in% setdiff(as.character(GWAS_frq$SNP),temp),]

	################################################################################################	
	####################	Remove The Positive control regions	################################

	risk_regions$start =  risk_regions$start * 1000000
	risk_regions$end = risk_regions$end * 1000000
	SNPs_region = NULL
	for(i in 1:nrow(risk_regions))
	{
		SNPs_region[[i]] = GWAS_frq[GWAS_frq$CHR == risk_regions$CHR[i] & GWAS_frq$pos >= risk_regions$start[i] & GWAS_frq$pos <= risk_regions$end[i] ,]
	}
	GWAS_frq = GWAS_frq[GWAS_frq$SNP %in% setdiff(as.character(GWAS_frq$SNP),unlist(lapply(SNPs_region,function(x) { as.character(x$SNP)}))),]

	logistic = read.table(paste(outfile,"logistic_",country,"_",Disease_name,"_",S,".assoc.logistic",sep=""),header=T)
	logistic = logistic[logistic$TEST == "ADD",]

	risk_regions$HD = 1:163
	found = NULL
	for(i in 1:nrow(Target_frq))
	{
	
		found[[i]] = risk_regions[risk_regions$CHR == Target_frq$CHR[i] & risk_regions$start <= as.numeric(Target_frq$pos[i]) & risk_regions$end >= as.numeric(Target_frq$pos[i]),]

	}
	for(i in 1:nrow(Target_frq))
	{
		if(nrow(found[[i]]) == 1)
		{
			Target_frq$HD[i] = found[[i]]$HD
		}
		else
		{
			Target_frq$HD[i] = NA
		}
	}
	#found = NULL
	#for(i in 1:nrow(risk_regions))
	#{
	#	found[[i]] = Target_frq[Target_frq$CHR == risk_regions$CHR[i] & Target_frq$pos >= (risk_regions$start[i]) & Target_frq$pos <= risk_regions$end[i] ,]
	#}
	#good_SNPs = NULL
	#count = 1
	#temp=NULL
	#for(i in 1:length(found))
	#{
	#	print(i)
	#	if(nrow(found[[i]]) >= 1)
	#	{
	##		temp=NULL
	#		temp = logistic[logistic$SNP %in% found[[i]]$SNP,]
	#		good_SNPs[count] = head(as.character(temp[order(temp$P),]$SNP),n=1)
	#		count = count + 1
	#	}
	#}
	#good_SNPs
	#Target_frq = Target_frq[Target_frq$SNP %in% good_SNPs,]
	#rownames(Target_frq ) = seq(1:nrow(Target_frq))
	################################################################################################

	GWAS_frq = GWAS_frq[GWAS_frq$MAF %in% Target_frq$MAF,]
	GWAS_frq = GWAS_frq[order(GWAS_frq$MAF),]
	Target_frq = Target_frq[order(Target_frq$MAF),]

	GWAS_hist = hist(GWAS_frq$MAF)
	Target_hist = hist(Target_frq$MAF,breaks=GWAS_hist$breaks)

	sim = list(matrix(NA,nrow=nrow(Target_frq ),ncol=6))
	all_simulation = rep(sim,1000)
	dev.off()
	################################################################################################
	################################################################################################
		system(paste("mkdir /home/ulg/genan/mahmoud/Immunochip/Results/SNP_Simulation/",country,"_",Disease_name,"_",S,sep=""))
	for(j in 1:1000)
	{
		simulation = NULL
		print(paste("Simulation no ",j,sep=""))
		for(i in 1:(length(Target_hist$breaks)-1)) {	
			if(i == 1){
				simulation = rbind(simulation,GWAS_frq[sample(x=c(1:sum(GWAS_hist$counts[1:i])),size=Target_hist$counts[i]),])
			} else {
				simulation = rbind(simulation,GWAS_frq[sample(x=c(sum(GWAS_hist$counts[1:i-1]):sum(GWAS_hist$counts[1:i])),size=Target_hist$counts[i]),])
			}	 
		}
		#write.table(simulation,file=paste('/home/ulg/genan/mahmoud/IBD_Morocan/simulations_positive/',Disease_name,'_simulation_',j,sep=""),row.names=F,quote=F)
		write.table(simulation,file=paste("/home/ulg/genan/mahmoud/Immunochip/Results/SNP_Simulation/",country,"_",Disease_name,"_",S,"/S_",S,"_",Disease_name,"_simulation_",j,sep=""),row.names=F,quote=F)
		all_simulation[[j]] = simulation
		simulation = NULL
	
	}

	################################################################################################
	################################################################################################


	f = file.path(paste("/home/ulg/genan/mahmoud/Immunochip/Results/SNP_Simulation/",country,"_",Disease_name,"_",S,sep=""),paste("S_",S,"_",Disease_name,"_simulation_",1:1000,sep=''))
	all_simulation <- lapply(f, read.table,header=T,stringsAsFactors=F)
	simulation_hist = NULL
	orignial_sum = sum(-log10(logistic[as.character(logistic$SNP) %in% as.character(Target_frq$SNP),9]))
	for(j in 1:1000)
	{
		print(j)
		simulation_hist[j] = sum(-log10(logistic[logistic$SNP %in% all_simulation[[j]]$SNP,9]))
	}
	print(orignial_sum)
	print((sum(simulation_hist >= orignial_sum) + 1)/(1000 + 1))
	temp = cbind(orignial_sum,(sum(simulation_hist >= orignial_sum) + 1)/(1000 + 1))
	write.table(simulation_hist,file=paste('/home/ulg/genan/mahmoud/Immunochip/Results/SNP_Simulation/',country,'_',Disease_name,'/S_',S,'_',Disease_name,'_result.txt',sep=""),row.names=F,quote=F)
	#png(paste('/home/ulg/genan/mahmoud/Immunochip/Results/SNP_Simulation/',country,'_',Disease_name,'/S_',S,'_',Disease_name,'_Simulation_hist.png',sep=""),width=290,height=210,units="mm",res=100)
	#hist(simulation_hist,main=paste(Disease_name,' 1000 simulation for 122 Loci',sep=""),xlab="Sum(log(1/P))",breaks=15)
	#abline(v=orignial_sum,lty=2,col="Red")
	#dev.off()
	
	save.image(paste(country,"_",Disease_name,"_",S,".RData",sep=""))
	temp = logistic[logistic$SNP %in% Target_frq$SNP,]
	temp = temp[order(temp$P),]
	res =  data.frame(matrix(NA,nrow=30,ncol=10))
	
	temp$SNP = as.character(temp$SNP)
	temp$A1 = as.character(temp$A1)
	temp$TEST = as.character(temp$TEST)
	for(i in 1:30)
	{
		
		print(i)
		simulation_hist = NULL
		res[i,1:9] = temp[1,1:9]
		removed_SNP = temp[1,2]
		temp = temp[-1,]
		orignial_sum = sum(-log10(temp[,9]))
		for(j in 1:1000)
		{
		
			remove = all_simulation[[j]][all_simulation[[j]]$MAF == Target_frq[ grep(pattern=paste("^",as.character(removed_SNP),"$",sep=""),x=as.character(Target_frq$SNP)),5],]
			count = 0
			while(nrow(remove) == 0)
			{
				#print(paste("Inside While ",j," ",count,sep=""))
				count = count + 1
				lower = (Target_frq[ grep(pattern=paste("^",as.character(removed_SNP),"$",sep=""),x=as.character(Target_frq$SNP)),5] - (count * 0.01))
				upper = (Target_frq[ grep(pattern=paste("^",as.character(removed_SNP),"$",sep=""),x=as.character(Target_frq$SNP)),5] + (count * 0.01))
				remove = all_simulation[[j]][all_simulation[[j]]$MAF >= lower & all_simulation[[j]]$MAF <= upper,]
			}
			count = 0
			all_simulation[[j]] = all_simulation[[j]][all_simulation[[j]]$SNP != remove[1,2],]
			simulation_hist[j] = sum(-log10(logistic[logistic$SNP %in% all_simulation[[j]]$SNP,9]))
		}
		res[i,10] = (sum(simulation_hist >= orignial_sum) + 1)/(1000 + 1)
		if(res[i,10] > 0.5)
		{
			break;
		}
	}
	



}





