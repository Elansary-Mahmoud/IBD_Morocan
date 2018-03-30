##############################################################
################### user Defined Function  ###################
##############################################################

check_inside_risk_region = function(CHR,START,END)
{
	risk = risk_regions[as.character(risk_regions$CHR) == as.character(CHR),]
	return_value = FALSE
	y1 = START
	y2 = END
	for(i in 1:nrow(risk))
	{
		x1 = (risk$start[i]*1000000)
		x2 = (risk$end[i]*1000000)
		if((x1 >= y1 & x1 <= y2) | (x2 >= y1 & x2 <= y2) | (y1 >= x1 & y1 <= x2) | (y2 >= x1 & y2 <= x2))
		{
			return_value = TRUE
			break;
		}
	}
	return_value
}

check_inside_used_region = function(CHR,START,END)
{
	risk = used_regions[as.character(used_regions[,1]) == as.character(CHR),]
	return_value = FALSE
	if(nrow(risk) > 0)
	{
		y1 = START
		y2 = END
		for(i in 1:nrow(risk))
		{
			x1 = risk[i,2]
			x2 = risk[i,3]
			if((x1 >= y1 & x1 <= y2) | (x2 >= y1 & x2 <= y2) | (y1 >= x1 & y1 <= x2) | (y2 >= x1 & y2 <= x2))
			{
				return_value = TRUE
				break;
			}
		}
	}
	return_value
}

sample_Region = function(i)
{
        print(i)
        restart = 1
        element_found = 0
        restart_times = 0
        while(restart == 1){
               random_start = sample(x=1:nrow(CHRs),size=1)
               #print(paste("Restart",restart_times))
               #restart_times = restart_times + 1
                random_end = 0
                random_chr = 0
                summation = 0
                for(j in random_start:(nrow(CHRs)-10))
                {
                        res = try(if(CHRs[j-1,1] == CHRs[j,1]){},TRUE)
                        if(inherits(res, "try-error")) {
                                print("skipping Error")
                                break;
                        } else {

                                if(CHRs[j-1,1] == CHRs[j,1]) {
                                summation = sum(CHRs[random_start:j,5])
                                if(results[i,2] >= (summation - 2) & results[i,2] <= (summation + 2) ) {
                                        restart = 0
                                        element_found = 1;
                                        random_end = j
                                        random_chr = CHRs[j,1]
					if(check_inside_risk_region(random_chr,CHRs[random_start,2],CHRs[random_end,3]))
					{
						print("window in risk region")
						restart = 1
						element_found = 0
						break;
					} else {
						if(nrow(used_regions) > 0)
						{
							if(check_inside_used_region(random_chr,CHRs[random_start,2],CHRs[random_end,3]))
							{
								print("window already sampled")
								restart = 1
								element_found = 0
								break;
							}
						}	
                                        	break;
                                        }
                                }
                                } else {
                                        break;
                                }
                        }

                }
                if(element_found == 1)
                {
                        restart = 0
                        output_saved = logistic[logistic$CHR == random_chr & logistic$BP >= CHRs[random_start,2] & logistic$BP <= CHRs[random_end,3],]
                        if(nrow(output_saved) > 0 )
                        {
                        	value = c(random_chr,output_saved$BP[1],output_saved$BP[nrow(output_saved)])
                        	used_regions[global_counter,] <<- value
                        	global_counter <<- global_counter + 1
                        	colnames(used_regions) = c("CHR","start","end")
				write.table(output_saved,file=paste(outfile2,country,"_",Disease_name,"_",S,"_simulation_s_",s_no,"_",results[i,1],sep=""),quote=F,row.names=F,col.names=F)                        
                        } else {
	                        restart = 1
	                        element_found = 0
                        }
                }
        }
}



used_regions <<- data.frame(matrix(NA,nrow=0,ncol=3))
global_counter <<- 1
Sys.time()
args <- commandArgs(trailingOnly = TRUE)
Disease_name = as.character(args[1])
country = as.character(args[2])
S = as.numeric(args[3])
print(paste(Disease_name,country,S))

plink = "/home/ulg/genan/mahmoud/plink-1.07-x86_64/plink"
positions = read.table('/home/ulg/genan/mahmoud/Immunochip/ibd_rel5_EU_clean.v3.HG19.updateID.bim',stringsAsFactor=F)
GWAS_frq = read.table(paste("/home/ulg/genan/mahmoud/Immunochip/Results/",country,"_",Disease_name,"_",S,"_Broad_GWAS_SNP_frequency.frq",sep=""),header=T)
GWAS_frq$pos = positions[,4]
Negative_SNP = read.table('Negative_SNPs_list.txt',header=T)
temp = intersect(as.character(GWAS_frq$SNP),as.character(Negative_SNP[,1]))
GWAS_frq = GWAS_frq[GWAS_frq$SNP %in% setdiff(as.character(GWAS_frq$SNP),temp),]
risk_regions = read.table('../IBD_Morocan/Regions.txt',header=T) #'risk_regions_justin.txt'

SNPs_region = NULL
for(i in 1:nrow(risk_regions))
{
	SNPs_region[[i]] = GWAS_frq[GWAS_frq$CHR == risk_regions$CHR[i] & GWAS_frq$pos >= (risk_regions$start[i] * 1000000) & GWAS_frq$pos <= (risk_regions$end[i] *  1000000) ,]
}
GWAS_frq = GWAS_frq[GWAS_frq$SNP %in% setdiff(as.character(GWAS_frq$SNP),unlist(lapply(SNPs_region,function(x) { as.character(x$SNP)}))),]

Leila_regions =  read.table('../IBD_Morocan/Top_IBD_SNPs.txt',header=T,sep="\t")
risk_regions = risk_regions[unique(Leila_regions$HD),]

logistic = read.table(paste("/home/ulg/genan/mahmoud/Immunochip/Results/","logistic_",country,"_",Disease_name,"_",S,".assoc.logistic",sep=""),header=T)
logistic = logistic[logistic$TEST == "ADD",]

#else {
#	load(paste(country,"_",Disease_name,"_",S,".RData",sep=""))
#}
gec = "/home/ulg/genan/mahmoud/GEC/gec/gec.jar"
New_Genofile = paste("/home/ulg/genan/mahmoud/Immunochip/Results/",country,"_",Disease_name,"_",S,sep="")
outfile = paste("/home/ulg/genan/mahmoud/Immunochip/Results/SNP_Simulation/",country,"_",Disease_name,"_",S,"/",sep="")
outfile2 = paste("/home/ulg/genan/mahmoud/Immunochip/Results/Region_Simulation/",country,"_",Disease_name,"_",S,"_v2/",sep="")


options(scipen=500)
results = data.frame(matrix(NA,nrow=nrow(risk_regions),ncol=2))
for(i in 1:nrow(risk_regions))
{
	print(i)
	#system(paste("java -jar ",gec," --no-web --effect-number --plink-binary ",New_Genofile," --regions 'chr",risk_regions[i,]$CHR,":",as.character(risk_regions[i,]$start*1000000),"-",as.character(risk_regions[i,]$end*1000000),"' --out ",outfile,"Original_region_",as.numeric(rownames(risk_regions[i,])),"_",country,"_",Disease_name,"_",S,sep=""))
	results[i,1] = as.numeric(rownames(risk_regions[i,]))
	results[i,2] = read.table(paste(outfile,"Original_region_",as.numeric(rownames(risk_regions[i,])),"_",country,"_",Disease_name,"_",S,".sum",sep=""),header=T)[2]
}



y="/home/ulg/genan/mahmoud/Immunochip/Results/"
write.table(as.character(GWAS_frq$SNP),file=paste(y,"extract_SNPs_",country,"_",Disease_name,"_",S,".txt",sep=""),quote=F,col.names=F,row.names=F)
system(paste(plink," --noweb --bfile ",New_Genofile," --extract ",y,"/extract_SNPs_",country,"_",Disease_name,"_",S,".txt --recode --make-bed --out ",New_Genofile,"_regions",sep=""))
system(paste("java -jar ",gec," --no-web --effect-number --plink-binary ",New_Genofile,"_regions"," --genome --out ",outfile,"ALL_Genome_",country,"_",Disease_name,"_",S,sep="")) 

x = file.path(paste("/home/ulg/genan/mahmoud/Immunochip/Results/SNP_Simulation/",country,"_",Disease_name,"_",S,sep=""),paste("ALL_Genome_",country,"_",Disease_name,"_",S,".block.txt",sep=""))
CHRs = read.table(x,header=T,check.names=F,sep="\t",row.names=NULL)

Sys.time()

if (file.exists(paste("/home/ulg/genan/mahmoud/Immunochip/Results/Region_Simulation/",country,"_",Disease_name,"_",S,"_v2",sep=""))){
	print("Directory already exists")
} else {
	print("Creating Directory")
	system(paste("mkdir /home/ulg/genan/mahmoud/Immunochip/Results/Region_Simulation/",country,"_",Disease_name,"_",S,"_v2",sep=""))
}
for(s in 1:100)
{
	s_no <<- s
	print(paste("Simulation number ",s,sep=""))
	lapply(1:nrow(results),sample_Region)
}

Sys.time()
regions_results = matrix(NA,nrow=nrow(results),ncol=100)
counter = 1

for(i in 1:nrow(results))
{	
	print(paste("HD Region",i,sep=""))
	x = paste(outfile2,country,"_",Disease_name,"_",S,"_simulation_s_",1:100,"_",results[i,1],sep="")
	#temp = NULL
	#temp = paste("wc -l ",x," | cut -f1 -d' '")
	#cmd = which(unlist(as.numeric(lapply(temp,system,intern=TRUE))) == 0)
	#if(length(cmd >= 1)){
	#	print("Error ")
	#	for(iterator in 1:length(cmd))
	#	{
			#sp = strsplit(strsplit(x[cmd],split="/")[[iterator]][10],split="_")[[1]][6:7]
	#		s_no <<- cmd[iterator]
	#		sample_Region(i)
	#	}
	#}	

	HD_hist = unlist(lapply(x,function(K) { min(read.table(K,header=F)[,9],na.rm=T) } ))
	regions_results[counter,] = HD_hist
	counter = counter + 1

}

Sys.time()
regions_results=na.omit(regions_results)
regions_results = -log10(regions_results)
temp=NULL
for(i in 1:ncol(regions_results))
{
	temp[i] = sum(regions_results[,i])
}
original_sum = 0
original_value = NULL
counter = 1
original_names = NULL

assoc.linear = list(NULL)
for(i in 1:nrow(results))
{
	print(i)
	assoc.linear[[i]] = logistic[logistic$CHR == risk_regions$CHR[i] & logistic$BP >= (risk_regions$start[i]*1000000) & logistic$BP <= (risk_regions$end[i]*1000000),]
}


for(i in 1:nrow(regions_results))
{
	original_sum = original_sum + (-log10(min(assoc.linear[[i]]$P,na.rm=T)))
	original_value[counter] = -log10(min(assoc.linear[[i]]$P,na.rm=T))
	original_names[counter] = results[i,1]
	counter = counter + 1
}
print(original_sum)
pvalue_regions = (sum(temp >= original_sum) + 1)/(100+1)
print(pvalue_regions)

regions_results = cbind(regions_results,original_value,original_names)
write.table(regions_results,paste(outfile2,country,"_",Disease_name,"_",S,"_Simulation_Regions_Results.txt",sep=""),quote=F,col.names=F)
temp = NULL
pvalues = NULL

temp = NULL
temp = regions_results
region = matrix(NA,nrow=(nrow(regions_results)-2),ncol=7)
name = NULL
Sys.time()
for(j in 1:(nrow(regions_results)-2))
{
	print(paste("HD Region ",temp[which.max(temp[,101]),102],sep=""))
	name = temp[which.max(temp[,101]),102]
	temp = temp[-which.max(temp[,101]),]
	sim = NULL
	for(i in 1:ncol(temp))
	{
		sim[i] = sum(temp[,i]) # sim[101] contain the originl sum - removed region
	}
	print((sum(sim[1:100] >= sim[101]) + 1)/(100+1))
	region[j,1] = name
	region[j,2] = (sum(sim[1:100] >= sim[101]) + 1)/(100+1)
	region[j,3] = sim[101]
	region[j,4] = min(sim[1:100])
	region[j,5] = max(sim[1:100])
	region[j,6] = sum(sim[1:100])
	region[j,7] = sum(sim[1:100] >= sim[101])
	
}
colnames(region) = c("Region_Name","pvalue","sum_original","min_simulated","max_simulated","sum_simulated","higher_points")
write.table(region,paste(outfile2,country,"_",Disease_name,"_",S,"_sig_region.txt",sep=""))
save.image(paste(country,"_",Disease_name,"_",S,"_Regions.RData",sep=""))
Sys.time()






