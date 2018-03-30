load("New_code_workspace.RData")
args <- commandArgs(trailingOnly = TRUE)
Disease_name = as.character(args[1])
print(Disease_name)


logistic = read.table(paste(Disease_name,'.assoc.logistic',sep=''),header=T)
logistic = logistic[logistic$TEST == "ADD",]


########################################################################################

results = matrix(NA,nrow=length(unique(Leila_list$HD))+length(notused),ncol=100)
counter = 1
counter2= 1
found = NULL
notfound = NULL
Regions = read.table(paste(Disease_name,"_Simulation_Regions_Results.txt",sep=""))
Regions = Regions[,103]
for(i in Regions)
{
                x = sample(SIMULATED_DATA[[i]],100)
                x = substr(x,1,nchar(x)-6)
                sim.HD.assoc.linear = lapply(x,read.table,header=T)
                HD_hist = unlist(lapply(1:length(sim.HD.assoc.linear),function(i) { min(logistic[logistic$SNP %in% as.character(sim.HD.assoc.linear[[i]][,2]),]$P,na.rm=T) } ))
                print(paste("HD Region",i,sep=""))
                results[counter,] = HD_hist
                counter = counter + 1
}
#results=na.omit(results)
results = -log10(results)
template=NULL
for(i in 1:ncol(results))
{
        template[i] = sum(results[,i])
}
original_sum = 0
original_value = NULL
counter = 1
original_names = NULL

f = file.path("/home/ulg/genan/mahmoud/IBD_Morocan/Regions_N",paste(Disease_name,"_HD_",1:163,".assoc.logistic.mperm",sep=""))
assoc.linear.mperm <- lapply(f[file.exists(f)], read.table,header=T,stringsAsFactors=F)

f = file.path("/home/ulg/genan/mahmoud/IBD_Morocan/Regions_N",paste(Disease_name,"_HD_",1:163,".assoc.logistic",sep=""))
assoc.linear <- lapply(f[file.exists(f)], read.table,header=T,stringsAsFactors=F)

for(i in Regions)
{
        #original_sum = original_sum + (-log10(min(assoc.linear[[i]]$P,na.rm=T)))
        original_value[counter] = -log10(min(assoc.linear[[i]]$P,na.rm=T))
        original_names[counter] = i
        counter = counter + 1
}
results = cbind(results,original_value,original_names)

original_sum = sum(results[,101])

print(original_sum)
pvalue_regions = (sum(template >= original_sum) + 1)/(100+1)
print(pvalue_regions)

temp = NULL
temp = results
region = matrix(NA,nrow=nrow(results),ncol=5)
name = NULL
for(j in 1:(nrow(results)-2))
{
        print(paste("HD Region ",temp[which.max(temp[,101]),102],sep=""))
        name = temp[which.max(temp[,101]),102]
        temp = temp[-which.max(temp[,101]),]
        sim = NULL
	if(nrow(temp) == 0) {
                break;
        } else {
        for(i in 1:ncol(temp))
        {
                sim[i] = sum(temp[,i]) # sim[101] contain the originl sum - removed region
        }
        print((sum(sim[1:100] >= sim[101]) + 1)/(100+1))
        region[j,1] = name
        region[j,2] = (sum(sim[1:100] >= sim[101]) + 1)/(100+1)
        region[j,3] = sim[101]
        region[j,4] = sum(sim[1:100])
        region[j,5] = max(sim[1:100])
	}
}
write.table(region,paste("FINAL/",Disease_name,"_sig_region_141.txt",sep=""))

png(paste("FINAL/",Disease_name,'_Simulation_hist_regions_141.png',sep=""),width=290,height=210,units="mm",res=100)
hist(template,main=paste(Disease_name,' 100 simulation for 141 Loci',sep=""),xlab="Sum(log(1/P))",xlim=c(min(template),original_sum+10))
abline(v=original_sum,lty=2,col="Red")
dev.off()


##############################################################

results = matrix(NA,nrow=length(unique(Leila_list$HD)),ncol=100)
counter = 1
counter2= 1
found = NULL
notfound = NULL
for(i in unique(Leila_list$HD))
{	
		x = sample(SIMULATED_DATA[[i]],100)
		x = substr(x,1,nchar(x)-6)
		sim.HD.assoc.linear = lapply(x,read.table,header=T)
		HD_hist = unlist(lapply(1:length(sim.HD.assoc.linear),function(i) { min(logistic[logistic$SNP %in% as.character(sim.HD.assoc.linear[[i]][,2]),]$P,na.rm=T) } ))
		print(paste("HD Region",i,sep=""))
		results[counter,] = HD_hist
		counter = counter + 1
}
#results=na.omit(results)
results = -log10(results)
template=NULL
for(i in 1:ncol(results))
{
	template[i] = sum(results[,i])
}
original_sum = 0
original_value = NULL
counter = 1
original_names = NULL

f = file.path("/home/ulg/genan/mahmoud/IBD_Morocan/Regions_N",paste(Disease_name,"_HD_",1:163,".assoc.logistic.mperm",sep=""))
assoc.linear.mperm <- lapply(f[file.exists(f)], read.table,header=T,stringsAsFactors=F)

f = file.path("/home/ulg/genan/mahmoud/IBD_Morocan/Regions_N",paste(Disease_name,"_HD_",1:163,".assoc.logistic",sep=""))
assoc.linear <- lapply(f[file.exists(f)], read.table,header=T,stringsAsFactors=F)

for(i in unique(Leila_list$HD))
{
	#original_sum = original_sum + (-log10(min(assoc.linear[[i]]$P,na.rm=T)))
	original_value[counter] = -log10(min(assoc.linear[[i]]$P,na.rm=T))
	original_names[counter] = i
	counter = counter + 1
}
results = cbind(results,original_value,original_names)

original_sum = sum(results[,101])

print(original_sum)
pvalue_regions = (sum(template >= original_sum) + 1)/(100+1)
print(pvalue_regions)

temp = NULL
temp = results
region = matrix(NA,nrow=nrow(results),ncol=5)
name = NULL
for(j in 1:(nrow(results)-2))
{
	print(paste("HD Region ",temp[which.max(temp[,101]),102],sep=""))
	name = temp[which.max(temp[,101]),102]
	temp = temp[-which.max(temp[,101]),]
	sim = NULL
	if(nrow(temp) == 0) {
		break;
	} else {
	for(i in 1:ncol(temp))
	{
		sim[i] = sum(temp[,i]) # sim[101] contain the originl sum - removed region
	}
	print((sum(sim[1:100] >= sim[101]) + 1)/(100+1))
	region[j,1] = name
	region[j,2] = (sum(sim[1:100] >= sim[101]) + 1)/(100+1)
	region[j,3] = sim[101]
	region[j,4] = sum(sim[1:100])
	region[j,5] = max(sim[1:100])
	}
}
write.table(region,paste("FINAL/",Disease_name,"_sig_region_121.txt",sep=""))
save.image(paste("FINAL/",Disease_name,"simulated_region.RData",sep=""))

png(paste("FINAL/",Disease_name,'_Simulation_hist_regions_121.png',sep=""),width=290,height=210,units="mm",res=100)
hist(template,main=paste(Disease_name,' 100 simulation for 121 Loci',sep=""),xlab="Sum(log(1/P))",xlim=c(min(template),original_sum+10))
abline(v=original_sum,lty=2,col="Red")
dev.off()


