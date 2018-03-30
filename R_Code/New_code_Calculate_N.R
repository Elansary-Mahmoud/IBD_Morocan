predict_n = function(Probe_index)
{
	if(length(assoc.linear[[Probe_index]][assoc.linear[[Probe_index]]$TEST == "ADD",9]) == length(assoc.linear.mperm[[Probe_index]]$EMP2)) 
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
	} else {
	print("Length problem")
	result = NA
	result
	}
}


res = matrix(0,nrow=163,ncol=2)
for(i in 1:163)
{
	print(paste("HD Region ",i,sep=""))
	res[i,1] = i
	res[i,2] = as.numeric(system(paste('ls /home/ulg/genan/mahmoud/IBD_Morocan/simulations_regions_v2/UC_simulation_s_*_',i,'.assoc.logistic | wc -l',sep=''),intern=TRUE)) + as.numeric(system(paste('ls /home/ulg/genan/mahmoud/IBD_Morocan/simulations_regions_v2/lem_sa_*_',i,'.txt | wc -l',sep=''),intern=TRUE)) 
}

pred_n = list()
Names = list()
for(i in 1:nrow(res))
{
	x = c(system(paste('ls /home/ulg/genan/mahmoud/IBD_Morocan/simulations_regions_v2/UC_simulation_s_*_',res[i,1],'.assoc.logistic.mperm',sep=''),intern=TRUE), system(paste('ls /home/ulg/genan/mahmoud/IBD_Morocan/simulations_regions_v2/lem_simulation_s_*_',res[i,1],'.assoc.logistic.mperm',sep=''),intern=TRUE) )
	y = NULL
	print(i)
	for(j in 1:length(x))
	{

		assoc.linear = NULL
		assoc.linear.mperm = NULL
		assoc.linear <<- lapply(substr(x[j],1,nchar(x[j])-6), read.table,header=T,stringsAsFactors=F)

		assoc.linear.mperm <<- lapply(x[j], read.table,header=T,stringsAsFactors=F)
		y[j] = predict_n(1)
	}
	pred_n[[i]] = y
	Names[[i]] = x
}

Regions = matrix(0,ncol=3,nrow=nrow(N))
for(i in 1:nrow(N))
{
	summation = 0
	print(i)
	if(is.na(N[i,2])) {
		Regions[i,1] =  i
		Regions[i,2] =  N[i,2]
		Regions[i,3] = NA
	} else {
		for(j in 1:length(pred_n))
		{
			summation = summation + sum((na.omit(pred_n[[j]]) > (N[i,2] - 2)) & (na.omit(pred_n[[j]]) < (N[i,2] + 2)))
		}
		Regions[i,1] =  i
		Regions[i,2] =  N[i,2]
		Regions[i,3] =  summation
	}
}
 Rep = summary(as.factor(Regions[,2]))
 names(Rep)[61] = "NA"
 template = Regions 
 template = cbind(template,rep(0,nrow(template)),rep(0,nrow(template)))
 for(i in 1:nrow(Regions))
 {
 	 template[i,4] = Rep[as.character(template[i,2])]
 	 template[i,5] = template[i,3]/template[i,4]
 }


SIMULATED_DATA = list
for(i in 1:nrow(template))
{
	index = as.numeric(row.names(na.omit(N[N[,2] == template[i,2],])))
	sum = NULL
	for(j in 1:length(index))
	{
		sum = c(sum,Names[[index[j]]])
	}
	SIMULATED_DATA[[i]] = sum	
}

save.image("New_code_workspace.RData")

