args <- commandArgs(trailingOnly = TRUE)
Disease_name = as.character(args[1])

Disease_name = "UC"
regions_b37 = read.table("justin_UCSC.txt",stringsAsFactor=F)
GWAS_frq =  read.table('GWAS_SNP_frequency.frq',header=T)
risk_regions = read.table('Regions.txt',header=T,sep="\t") #'risk_regions_justin.txt'
positions = read.table('AllImmunochipsMarocains-29102012-build37.bim',stringsAsFactor=F)
GWAS_frq$pos = positions[,4]
Negative_SNP = read.table('Negative_SNPs_list.txt',header=T)
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
Target_frq = read.table('all_SNP_region_frequency.frq',header=T,stringsAsFactor=F)
Target_frq = Target_frq[Target_frq$SNP %in% setdiff(Target_frq$SNP,Negative_SNP[,1]),]
UC_logistic = read.table(paste(Disease_name,".assoc.logistic",sep=""),header=T)
UC_logistic = UC_logistic[UC_logistic$TEST == "ADD",]
for(i in 1:nrow(Target_frq))
{
	Target_frq$pos[i] = positions[positions[,2] == Target_frq[i,2],4]
}
#risk_regions$HD = 1:163
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
found = NULL
for(i in 1:nrow(risk_regions))
{
	found[[i]] = Target_frq[Target_frq$CHR == risk_regions$CHR[i] & Target_frq$pos >= (risk_regions$start[i]) & Target_frq$pos <= risk_regions$end[i] ,]
}
good_SNPs = NULL
count = 1
temp=NULL
for(i in 1:length(found))
{
	print(i)
	if(nrow(found[[i]]) >= 1)
	{
		temp=NULL
		temp = UC_logistic[UC_logistic$SNP %in% found[[i]]$SNP,]
		good_SNPs[count] = head(as.character(temp[order(temp$P),]$SNP),n=1)
		count = count + 1
	}

}
good_SNPs
Target_frq = Target_frq[Target_frq$SNP %in% good_SNPs,]
Target_frq[122,]  = found[[25]][2,]
rownames(Target_frq ) = seq(1:nrow(Target_frq))
################################################################################################

GWAS_frq = GWAS_frq[GWAS_frq$MAF %in% Target_frq$MAF,]
GWAS_frq = GWAS_frq[order(GWAS_frq$MAF),]
Target_frq = Target_frq[order(Target_frq$MAF),]

GWAS_hist = hist(GWAS_frq$MAF)
Target_hist = hist(Target_frq$MAF,breaks=GWAS_hist$breaks)

sim = list(matrix(NA,nrow=nrow(Target_frq ),ncol=6))
all_simulation = rep(sim,1000)

################################################################################################
################################################################################################
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
	#write.table(simulation,file=paste('/root/Leila/Vega/IBD_Morocan/simulations_positive/',Disease_name,'_simulation_',j,sep=""),row.names=F,quote=F)
	all_simulation[[j]] = simulation
	simulation = NULL
	
}

################################################################################################
################################################################################################

#Disease_name = "IBD"
f = file.path("/home/ulg/genan/mahmoud/IBD_Morocan/simulations_positive",paste(Disease_name,"_simulation_",1:1000,sep=''))
all_simulation <- lapply(f, read.table,header=T,stringsAsFactors=F)
simulation_hist = NULL
orignial_sum = sum(-log10(UC_logistic[as.character(UC_logistic$SNP) %in% as.character(Target_frq$SNP),9]),na.rm=TRUE)
for(j in 1:1000)
{
	print(j)
	simulation_hist[j] = sum(-log10(UC_logistic[UC_logistic$SNP %in% all_simulation[[j]]$SNP,9]),na.rm=TRUE)
}
print(orignial_sum)
print((sum(simulation_hist >= orignial_sum) + 1)/(1000 + 1))
png(paste(Disease_name,'_Simulation_hist.png',sep=""),width=290,height=210,units="mm",res=450)
hist(simulation_hist,main=paste(Disease_name,' 1000 simulation for 122 Loci',sep=""),xlab="Sum(log(1/P))",breaks=15)
abline(v=orignial_sum,lty=2,col="Red")
dev.off()


f = file.path("/home/ulg/genan/mahmoud/IBD_Morocan/simulations_positive",paste("simulation_",1:1000,sep=''))
f2 = file.path("/home/ulg/genan/mahmoud/IBD_Morocan/simulations_positive",paste("IBD_simulation_",1:1000,sep=''))
f3= rbind(f,f2)
all_simulation <- lapply(f3, read.table,header=T,stringsAsFactors=F)


temp = UC_logistic[UC_logistic$SNP %in% Target_frq$SNP,]
temp = temp[order(temp$P),]
results = temp
for(i in 1:120)
{
	simulation_hist = NULL
	print(as.character(temp[i,2]))
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
		simulation_hist[j] = sum(-log10(UC_logistic[UC_logistic$SNP %in% all_simulation[[j]]$SNP,9]))
	}
        results[i,10] = (sum(simulation_hist >= orignial_sum) + 1)/(1000 + 1)
	write.table(results[i,],"output_UC.txt",sep=" ",quote=F,append=TRUE,col.names=F,row.names=F)
	print((sum(simulation_hist >= orignial_sum) + 1)/(1000 + 1))
}

temp = UC_logistic[UC_logistic$SNP %in% Target_frq$SNP,]
temp = temp[order(temp$P),]
write.table(temp,file=paste("/root/Leila/Vega/IBD_Morocan/Results/",Disease_name,"_Target.txt",sep=""),quote=F,row.names=F)

Target_frq = read.table('all_SNP_region_frequency.frq',header=T,stringsAsFactor=F)
Target_frq = Target_frq[Target_frq$SNP %in% setdiff(Target_frq$SNP,Negative_SNP[,1]),]
for(i in 1:nrow(Target_frq))
{
	Target_frq$pos[i] = positions[positions[,2] == Target_frq[i,2],4]
}
#risk_regions$HD = 1:163
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

found = NULL
for(i in 1:nrow(risk_regions))
{
	found[[i]] = Target_frq[Target_frq$CHR == risk_regions$CHR[i] & Target_frq$pos >= (risk_regions$start[i]) & Target_frq$pos <= risk_regions$end[i] ,]
}

Target_frq[Target_frq$HD %in% which(unlist(lapply(found,nrow))>1),]




""""""""""""""""""""""""""""""""""""

regions_b37 = read.table("",stringsAsFactor=F)
UC_logistic = read.table('UC.assoc.logistic',header=T)
UC_logistic = UC_logistic[UC_logistic$TEST == "ADD",]

results = data.frame(matrix(NA,nrow=163,ncol=9))
for(i in 1:163)
{
	regexp = "[[:digit:]]+"
	chr = str_extract(strsplit(regions_b37[i,1],":")[[1]][1],regexp)
	pos = strsplit(strsplit(regions_b37[i,1],":")[[1]][2],"-")[[1]][1]
	print(i)
	restart = 1
	size = 0
	while(restart == 1)
	{
		if(nrow(UC_logistic[UC_logistic$CHR == as.numeric(chr) & UC_logistic$BP >= (as.numeric(pos) - size) & UC_logistic$BP <= (as.numeric(pos) + size),]) >= 1) {
			results[i,] = UC_logistic[UC_logistic$CHR == as.numeric(chr) & UC_logistic$BP >= (as.numeric(pos) - size) & UC_logistic$BP <= (as.numeric(pos) + size),]
			restart = 0
		} else {
			restart	= 1
			size = size + 30
		}
		if(size >= 150){
			restart	= 0
		}	
	}
}



#################################################################################################### 
######################################### code not for Run ##########################################
unlist(lapply(1:nrow(b37),function(i) { strsplit(strsplit(as.character(b37),":")[[1]][2],"-")[[1]][1]}))
start=NULL
for(i in 1:length(chr))
{
	print(i)
	temp = strsplit(b37[i,1],":")[[1]][2]
	start[i] = strsplit(temp,"-")[[1]][1]
}

chr= unlist(lapply(1:nrow(b37),function(i) {strsplit(b37[i,1],":")[[1]][1]})) 
SNP=map[,2]
chr = unlist(lapply(1:length(chr),function(i) {strsplit(chr[i],"chr")[[1]][2]}))

test = data.frame(,,start=unlist(lapply(1:nrow(b37),function(i) {strsplit(strsplit(as.character(b37),":")[[1]][2],"-")[[1]][1]})))
test = data.frame(chr,SNP,start,stringsAsFactors=F)


Regions = read.table('Regions.txt',header=T,stringsAsFactors=F)
b37 = read.table('AllImmunochipsMarocains-29102012-build37.bim',stringsAsFactors=F)
lookup=read.table('lookup_positive_regions.txt',header=T,stringsAsFactors=F)


for(i in 1:nrow(Regions_v2))
{
		Regions_v2$b37_start[i] = as.character(b37[grep(pattern=Regions_v2$IC[i],x=as.character(b37[,2])),4])
		print(i)
}


#################################################################################################### 
######################################### liab prediction ##########################################
####################################################################################################

Target_frq$Risk_Allele = ""
for(i in 1:nrow(risk_regions))
{
	if(nrow(Target_frq[Target_frq$CHR == risk_regions$CHR[i] & Target_frq$pos == risk_regions$start_b37[i],]) > 1)
	{
		Target_frq[Target_frq$CHR == risk_regions$CHR[i] & Target_frq$pos == risk_regions$start_b37[i],]$Risk_Allele = risk_regions$Risk_Allele[i]
	}
}


write.table(test[!(rownames(test) %in% c("25","25.1")),c(3,9)],"liab_prediction/risk_snps_reference.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)

write.table(test[!(rownames(test) %in% c("25","25.1")),c(3)],"liab_prediction/risk_snps.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)



snpgdsBED2GDS(bed.fn= "risk_snp_list.bed", fam.fn = "risk_snp_list.fam", bim.fn = "risk_snp_list.bim", "test.gds")

(genofile <- snpgdsOpen("GDR.gds"))
mat <- snpgdsGetGeno(genofile)

convert = data.frame(A1=c("A","T","C","G"),A2=c("T","A","G","C"))
convert[,1]=as.character(convert[,1])
convert[,2]=as.character(convert[,2])
risk[,2] = as.character(risk[,2])
for(i in which(bim[,5] != risk[,2]))
{
	print(i)
	if(risk[i,2] != as.character(bim[i,6]))
	{
		risk[i,2] = convert[convert$A1 == risk[i,2],2]
	}
}



liab = cbind(apply(2-mat,1,function(x) { sum(x,na.rm=TRUE) } ),fam[,4])
x11()
d1 = density(liab[liab[,2] == 1,1])
d2 = density(liab[liab[,2] == 2,1])


plot(range(d1$x, d2$x), range(d1$y, d2$y), type = "n", xlab = "x",ylab = "Density")
lines(d1, col = "blue")
lines(d2, col = "red")
