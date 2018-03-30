for(i in 65:nrow(res))
{
        if(res[i,2] >= 100) {
                x = c(system(paste('ls /home/ulg/genan/mahmoud/IBD_Morocan/simulations_regions_v2/UC_simulation_s_*_',i,'.assoc.logistic',sep=''),intern=TRUE), system(paste('ls /home/ulg/genan/mahmoud/IBD_Morocan/simulations_regions_v2/lem_simulation_s_*_',i,'.assoc.logistic',sep=''),intern=TRUE) )
                if(length(x) > 0){
                              } else {
                        print(paste("HD Region",i,sep=""))
                        print("Not Found")
                }

        }
}

Disease_name="CD"
Disease_name="CD"
load(paste(Disease_name,"_workspace.Rdata",sep=""))

res = res[-c(66,71,73,77,100,124,137,152,153,158),]
results=na.omit(results)
results = -log10(results)
temp=NULL

temp=NULL
for(i in 1:ncol(results))
{
        temp[i] = sum(results[,i])
}
original_sum = 0
original_value = NULL
counter = 1
original_names = NULL

for(i in which(res[,2] >= 100))
{
        original_sum = original_sum + (-log10(min(assoc.linear[[i]]$P,na.rm=T)))
        original_value[counter] = -log10(min(assoc.linear[[i]]$P,na.rm=T))
        original_names[counter] = i
        counter = counter + 1
}
print(original_sum)
pvalue_regions = (sum(temp >= original_sum) + 1)/(100+1)
print(pvalue_regions)
results = cbind(results,original_value,original_names)




counter = 1
counter2 = 1
found=NULL
notfound=NULL
for(i in 1:nrow(res))
{
	print(i)
        if(res[i,2] >= 100) {
                x = c(system(paste('ls /home/ulg/genan/mahmoud/IBD_Morocan/simulations_regions_v2/UC_simulation_s_*_',i,'.assoc.logistic',sep=''),intern=TRUE), system(paste('ls /home/ulg/genan/mahmoud/IBD_Morocan/simulations_regions_v2/lem_simulation_s_*_',i,'.assoc.logistic',sep=''),intern=TRUE) )
                if(length(x) > 0){
                        found[counter] = res[i,1]
                        counter = counter + 1
                } else {
			notfound[counter2] = res[i,1]
			counter2 = counter2 + 1
	                        
                }

        }
}

