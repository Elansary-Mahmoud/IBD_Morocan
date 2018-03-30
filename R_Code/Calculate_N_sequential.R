library('doParallel')
library('foreach')
library('doMC')
################### user Defined Function  ###################
##############################################################
predict_n = function(Probe_index)
{
	temp = data.frame(X=assoc.linear[[Probe_index]][assoc.linear[[Probe_index]]$TEST == "ADD",9],Y=assoc.linear.mperm[[Probe_index]]$EMP2)
	temp = temp[temp[,2] <= 0.95 & temp[,2] >= 0.05 , ]
	temp = temp[temp[,1] < 1, ]
	temp = data.frame(X=log10(1 - temp[,1]),Y=log10(1 - temp[,2]))
	#is.na(temp) <- sapply(temp, is.infinite)
	if(sum(is.na(temp)) != nrow(temp)) {
		result = as.double(coef(lm(Y ~ X,data=temp,na.action = na.exclude))[2])
		if(is.na(result)){
			print("No enough points")
			result
		}else{
			if(result > nrow(assoc.linear[[Probe_index]])) {
				result = nrow(assoc.linear[[Probe_index]])
				result
			} else {
				result
			}
		}
	} else {
		print("All points are NA")
		result = NA
		result
	}
}

##############################################################
##############################################################
args <- commandArgs(trailingOnly = TRUE)
Disease_name = as.character(args[1])
start = as.numeric(args[2])
end = as.numeric(args[3])
print(start)
print(Disease_name)

Disease_name = "CD"
start = 1
end = 10

Disease_code = as.numeric(switch(Disease_name, "IBD"="1", "CD"="2", "UC"="3") )

f = file.path("/home/ulg/genan/mahmoud/IBD_Morocan/Regions_N",paste(Disease_name,"_HD_",1:163,".assoc.logistic.mperm",sep=""))
assoc.linear.mperm <- lapply(f[file.exists(f)], read.table,header=T,stringsAsFactors=F)

f = file.path("/home/ulg/genan/mahmoud/IBD_Morocan/Regions_N",paste(Disease_name,"_HD_",1:163,".assoc.logistic",sep=""))
assoc.linear <- lapply(f[file.exists(f)], read.table,header=T,stringsAsFactors=F)


N = unlist(lapply(1:length(f[file.exists(f)]),predict_n))
library('stringr')
regexp = "[[:digit:]]+"
results = data.frame(HD_regions=paste("HD_",as.numeric(unlist(lapply(1:length(f[file.exists(f)]),function(i){ str_extract(strsplit(f[file.exists(f)],"/")[[i]][8],regexp) }))),sep=""),predicted_N=floor(N))
write.table(results,file=paste("/home/ulg/genan/mahmoud/IBD_Morocan/",Disease_name,"_pos_regions_pred_N.txt"),quote=F,row.names=F)

###########################################################################
regions_b37 = read.table("justin_UCSC.txt",stringsAsFactor=F)
GWAS_frq =  read.table('GWAS_SNP_frequency.frq',header=T)
risk_regions = read.table('Regions.txt',header=T) #'risk_regions_justin.txt'
positions = read.table('AllImmunochipsMarocains-29102012-build37.bim')
GWAS_frq$pos = positions[,4]

for(i in 1:nrow(risk_regions))
{
	start_region = as.numeric(strsplit(as.character(risk_regions[i,2]),"-")[[1]][1])
	end_region = as.numeric(strsplit(as.character(risk_regions[i,2]),"-")[[1]][2])
	size_region = (end_region - start_region) * 1000000
	print(i)
	center = as.numeric(strsplit(strsplit(regions_b37[i,1],":")[[1]][2],"-")[[1]][1])
	GWAS_frq = GWAS_frq[! GWAS_frq$SNP %in% as.character(GWAS_frq[GWAS_frq$pos >= ((center) - (size_region/2)) & GWAS_frq$pos <= ((center) + (size_region/2))  & GWAS_frq$CHR == risk_regions[i,1],2]),]	
	print(nrow(GWAS_frq))
	print("##################")
}

Negative_SNP = read.table('Negative_SNPs_list.txt',header=T)
Target_frq = read.table('all_SNP_region_frequency.frq',header=T)
Target_frq = Target_frq[Target_frq$SNP %in% setdiff(Target_frq$SNP,Negative_SNP[,1]),]
rownames(Target_frq ) = seq(1:nrow(Target_frq))

################################################################################################
####################	Remove The negative control regions	################################
temp = intersect(as.character(GWAS_frq$SNP),as.character(Negative_SNP[,1]))
GWAS_frq = GWAS_frq[GWAS_frq$SNP %in% setdiff(as.character(GWAS_frq$SNP),temp),]
################################################################################################

logistic = read.table(paste(Disease_name,'.assoc.logistic',sep=''),header=T)
logistic = logistic[logistic$TEST == "ADD",]
logistic = logistic[logistic$SNP %in% GWAS_frq$SNP,]
	
simulation_no=1
assoc.linear = NULL
assoc.linear.mperm = NULL
simulation = NULL
sim = list(matrix(NA,nrow=nrow(results),ncol=1))
all_simulation = rep(sim,100)
sumulation_sum = NULL
restart_times = 0
plink = "/home/ulg/genan/mahmoud/plink-1.07-x86_64/plink"
Genofile = "/home/ulg/genan/mahmoud/IBD_Morocan/AllImmunochipsMarocains-29102012-build37" 
phenofile = "/home/ulg/genan/mahmoud/Leila/Fichier_Phenotypes_final_bis.txt"
outfile = "/home/ulg/genan/mahmoud/IBD_Morocan/simulations_regions_v2/"
################### user Defined Function  ###################
##############################################################
sample_Region = function(i)
{
		output = NA
		print(paste("searching for ",i," Region with N = ",results[i,2],sep=""))
		if(!is.na(results[i,2])){
			add_new_SNP = 0
			restart = 1
			restart_times = 0
			while(restart == 1){
				if(restart_times > 2){
					restart = 0
				}
				
				if(add_new_SNP > 0){
					print(paste("extending window size by ",add_new_SNP,sep=""))
				} else {
					restart_times = restart_times + 1
					print(paste(restart_times," time restart ",sep=""))
					random_start = sample(x=1:nrow(GWAS_frq),size=1)
				}
				CHR=GWAS_frq[random_start,1]
				random_region=GWAS_frq[random_start:(random_start + results[i,2]  + add_new_SNP ),]
				
				if(length(summary(as.factor(random_region[,1]))) == 1){
					system(paste(plink," --noweb --bfile ",Genofile," --pheno ",phenofile," --mpheno ",Disease_code," --chr ",as.character(CHR)," --from-bp ",as.character(random_region[1,7])," --to-bp ",as.character(random_region[nrow(random_region),7]),"  --logistic --mperm 1000 --out ",outfile,"simulation_s_",s_no,"_",i,sep=""),ignore.stdout=T)	
					assoc.linear = NULL
					assoc.linear.mperm = NULL	
					f = file.path("/home/ulg/genan/mahmoud/IBD_Morocan/simulations_regions_v2",paste("simulation_s_",s_no,"_",i,".assoc.logistic",sep=""))
					assoc.linear <<- lapply(f[file.exists(f)], read.table,header=T,stringsAsFactors=F)
		
					f = file.path("/home/ulg/genan/mahmoud/IBD_Morocan/simulations_regions_v2",paste("simulation_s_",s_no,"_",i,".assoc.logistic.mperm",sep=""))
					assoc.linear.mperm <<- lapply(f[file.exists(f)], read.table,header=T,stringsAsFactors=F)
					pred_n = predict_n(1)
					if(!is.na(pred_n)){
					 	if(pred_n >= (results[i,2] - 2) & pred_n <= (results[i,2] + 2) ) {
							restart = 0
							print(paste("HD ",i," Region matching N = ",pred_n,sep=""))
							output = -log10(min(logistic[logistic$SNP %in% random_region$SNP,9],na.rm=T))
							write.table(output,file=paste("/home/ulg/genan/mahmoud/IBD_Morocan/simulations_regions_v2/sa_",s_no,"_",i,".txt",sep=""),quote=F)	 	
					 	} else{	
					 		restart = 1	 	
					 		add_new_SNP = add_new_SNP + 50
					 		if(add_new_SNP > 200)
					 		{
					 			add_new_SNP = 0
					 			restart = 1
					 		}
					 	}
				 	} else {
				 		restart	= 1
						add_new_SNP = 0
				 	}
				 	
				} else{
					restart	= 1
					add_new_SNP = 0
				}
			}
		} else {
			print(paste("Region ",i," has NA",sep=""))
		}
	output
}
##########################################################

for(s in start:(end - 1))
{
	s_no <<- s
	print(paste("Simulation number ",s,sep=""))
	x <- foreach(i=1:163, .combine='c') %do% sample_Region(i)
	#all_simulation[[s]] = aaply(9:nrow(results),1,sample_Region, .parallel=TRUE)	
	#write.table(sum(x),file=paste("/home/ulg/genan/mahmoud/IBD_Morocan/simulations_regions/sa_",s,".txt"),quote=F)
}





