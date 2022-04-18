# File:         get_boot_data.R
# Author:       Kelly Chang
# Date:         Oct 2017
# Version:      1.0
# 
# Description:  Helper R function which computes the mean fractional current
#               from Milnes protocol data. The data are resampled according
#               to the specified cell indices.
#

get_boot_data<-function(filepath, idx,drug){
	# read original data
	#filepath="data1/BUP.csv"
	#ass
	datadf<-read.csv(filepath)
#	print(datadf$id)
#	datadf_ass=datadf0[,c("time","sweep","conc","Cmean","exp","frac1","CBax","exp2","id")]
#	#diss
#	Diss_datadf0<-read.csv("/home/mohammadreza.samieegohar/john/hERG_fitting/data/all_Diss_sub_NS_Frac1.csv")
#	Diss_datadf=Diss_datadf0[,c("time","sweep","conc","Cmean","exp","frac1","CBax","exp2","id")]
##	print(head(Diss_datadf))
##	print(head(datadf_ass))
#	datadf=rbind(datadf_ass,Diss_datadf)
#	print(head(datadf))
#	print(nrow(datadf))
#	print(nrow(Diss_datadf))
#	print(nrow(datadf_ass))
#	
	
	#head(datadf)
	#datadf<-read.csv(sprintf("%stest.csv",filepath))
	
	
	#datadf0<-VdataR5
	#    datadf <- data.frame(apply(datadf0, 2, function(x) as.numeric(as.character(x))))
	#    datadf$Drug=VdataR5$Drug
	#    datadf$file=VdataR5$file
	#    datadf$BT=VdataR5$BT
	#	datadf$conc=VdataR5$conc
	
	# print(datadf)
	#    expdf<-unique(datadf[,c("conc","exp")])
	#    expdf<-expdf[with(expdf,order(conc,exp)),]
	expdf<-unique(datadf[,c("id","exp2")])
	# expdf21<-unique(datadf[,c("id","conc","exp2")])
	#print(expdf)
	expdf<-expdf[with(expdf,order(id,exp2)),]
	# bootstrap data
	#print(unique(expdf$id))
	print(unique(expdf$id))
	print(unique(idx))
	
	cellsdf<-expdf[idx,]
	#print(unique(expdf$id))
	
	#	print(unique(expdf$conc))
	#	print("hhh")
	#	print(unique(cellsdf$conc))
	#	print(idx)
	print(unique(cellsdf$id))
	
	cellsdf=na.omit(cellsdf)
	#print(unique(expdf$id))
	#print(unique(cellsdf$id))
	
	boot_list<-lapply(1:nrow(cellsdf),
			function(i){
				datadf[datadf$id==cellsdf[i,"id"] &
								datadf$exp2==cellsdf[i,"exp2"],
						c("time","sweep","conc","Cmean","exp","frac1","CBax")]}
	)
	#print(head(boot_list))
	#print("boooo")
	
	
	
	cellsdf$id<-as.character(cellsdf$id) # convert to characters for indexing
	
	#print(head(cellsdf))
	# get mean traces with selected points
	fracdata<-list()
	
	for(id in unique(cellsdf$id)){
		#print(id)
		expidx<-which(cellsdf$id==id)
		concdf<-do.call(rbind,boot_list[expidx])
		
		# compute mean of bootstrap data
		#meandf<-aggregate(concdf[,"frac",drop=FALSE], by=list(time=concdf$time, sweep=concdf$sweep), FUN=mean)
		meandf=concdf
		fracdata[[id]]<-list()
		sweeps<-unique(meandf$sweep)
		#if(any(sweeps!=1:length(sweeps))) stop("Bad sweeps!")
		for(i in sweeps){
			sweepdf<-meandf[meandf$sweep==i,c("time","conc","Cmean","exp","frac1","CBax")]
			sweepdf<-sweepdf[order(sweepdf$time),]
			
			# select starting point from first sweep
			if(i==1) startpt<-min(c(1,which(sweepdf$frac<=20)))
			fracdata[[id]][[i]]<-sweepdf[startpt:nrow(sweepdf),]
			#fracdata[[id]]<-sweepdf[startpt:nrow(sweepdf),]
			
		}# for each sweep
	}# for each dose0
#	print(fracdata$id)
#	prin("ddddddddddddd")
	fracdata
}
