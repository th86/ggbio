#Package:	Biological Figure Rendering Package using ggPlot	
#Author:	Tai-Hsien Ou Yang
#Description:	Functions rendering customized heatmaps and aligned heatmaps


library(gridExtra)
library(reshape2)
library(ggplot2)
library(scales)

library("IlluminaHumanMethylation27k.db")

#load("panCan12v47.notintersect.rda")
#load("allattractome.rda")

#panCanHeatmap.align(  synObjList.f,  lymList, mpList,'ls', "heatmap.lymxmp.tiff",  datasetList=datasetList[,4]    )
#panCanHeatmap.align(  synObjList.f,  lymList, mmList,'ls', "heatmap.lymxmm.tiff",  datasetList=datasetList[,4]    )


cinList=c("TPX2","NCAPG","KIFC1","KIF4A","NCAPH","HJURP","CENPA","CCNA2","KIF2C","BUB1B")
lymList=c("CD53","SASH3","LAPTM5","PTPRC","CD3E","LAIR1","NCKAP1L","SPI1","HAVCR2","CD2")
mtList=c("COL3A1","COL1A2","AEBP1","COL5A1","COL1A1","THBS2","COL6A3","FBN1","COL5A2","EMILIN1")
#geneList<-allattractome$rnaseq$'AHSA2|130872'[,1]
#ashaList=getGeneSymbols(geneList)
#geneList<-allattractome$rnaseq$'ANG'[,1]
#angList=getGeneSymbols(geneList)

#mpList=getGeneSymbols(allattractome$meth$'NA_NA_29848860'[,1])
#mmList=getGeneSymbols(allattractome$meth$'BIN2_12_50003941'[c(1,3:9),1])
#mmList[2]="PCED1B_12_45896487"
#names(mmList[2])=mmList[2]


getGeneSymbols = function(innames){
	outnames = sapply(innames, function(x){
	
		if(regexpr("\\?", x) > 0){
			o = strsplit(x, "\\|")[[1]][2]
		}else{
			o = strsplit(x, "\\|")[[1]][1]
		}
		return (o)

	}
	)
}

cinList=c("TPX2","NCAPG","KIFC1","KIF4A","NCAPH","HJURP","CENPA","CCNA2","KIF2C","BUB1B")
lymList=c("CD53","SASH3","LAPTM5","PTPRC","CD3E","LAIR1","NCKAP1L","SPI1","HAVCR2","CD2")
mtList=c("COL3A1","COL1A2","AEBP1","COL5A1","COL1A1","THBS2","COL6A3","FBN1","COL5A2","EMILIN1")
geneList<-allattractome$rnaseq$'AHSA2|130872'[,1]
ashaList=getGeneSymbols(geneList)
geneList<-allattractome$rnaseq$'ANG'[,1]
angList=getGeneSymbols(geneList)



#panCanHeatmap( synObjList.f,  cinList, "heatmap.cin.tiff")

#Render one heatmap
heatmap.metagene<-function(ge, geneList, figTitle , removeOutlier=TRUE , showTitle=TRUE, showTick=TRUE ){ 

	
	metagene =colSums(ge)
	#remove outliers
	if(removeOutlier==TRUE){
	outlier  = names(boxplot(metagene)$out )
	dev.off()
	select=setdiff(colnames(ge),outlier  )
	ge=ge[,select]
	
	#Rank to metagene 
	#ge=ge[,order( metagene[select],decreasing=FALSE)]
	}else{
	#Rank to metagene 
	#ge=ge[,order( metagene,decreasing=FALSE)]
	}


	#row-wise normalization
	ge=t(apply(ge, 1, function(x) (   (x-min(x)) / (max(x)-min(x)) ) )  )
	#ge=t(apply(ge, 1, function(x) (   (x) / (max(x)) ) )  )


	#Color distribution
	qn = quantile(ge, c(0.1, 0.9), na.rm = FALSE)	
	qn01 <- rescale (c(qn, range (ge))) 	


	ge=data.frame(ge) #to use $
	
	#Conver to long format
	gec <- cbind(ge,gene=rownames(ge))
	gem <- melt(gec,id.vars=c("gene"))
	
	#Render the heatmap

	if(showTitle==FALSE)
	figTitle=NULL

	if(showTick==FALSE)
	geneList=NULL
	
	geneList=substr(geneList,1,6)


	


	heatmap.gg<-ggplot(gem, aes(x=variable,y=gene)) +
		geom_tile(aes(fill=value)) +
		#scale_fill_gradientn(colours=c("blue","white","red")) +
		scale_fill_gradientn (
      			colours=colorRampPalette (c ("blue", "white", "red")) (20),
      			values = c (0, seq (qn01 [1], qn01 [2], length.out = 18), 1)) +
		theme( axis.text.x = element_blank(),  axis.ticks= element_blank(), legend.position = "none") +
		theme(aspect.ratio=0.1)+
		scale_x_discrete(name="")+
		scale_y_discrete(name="", labels=geneList)+ # ggplot cannot change the axis to the right
	         #annotation_custom(grob = textGrob(label = figTitle , hjust = 0, gp = gpar(cex = 0.5)), ymin=5, ymax=5, xmin=-20,xmax=-20)+
	       	ggtitle(figTitle)+
		theme(plot.title=element_text(size=2))+
		theme(plot.margin= unit(c(0,0,-1,0), "lines"))+ 
		theme(axis.text.y= element_text( angle = 0,hjust = 0,  size = 2 )) #gene names
		#theme(axis.title.y = element_text(angle = 90,hjust = 1,  size = 1 ) ) 
		#http://docs.ggplot2.org/0.9.3/theme.html
		
	return(heatmap.gg)
	}



#Render PanCan Heatmaps

panCanHeatmap<-function( synObjList, geneList, opfileName,datasetList =NULL , saveTiff=TRUE ){

	if(length(datasetList)==0){
	datasetList=c('LUSC','READ','GBM','BLCA','UCEC','COAD','OV','LAML','HNSC','LUAD','BRCA','KIRC') #Cancer Types
	}

	ggheatmap=list()

	for(i in as.numeric(order(datasetList))){ #Sort by alphabets	
	ggheatmap[[i]]<-heatmap.metagene(synObjList[[i]]$e[geneList,], geneList, datasetList[i])
	cat("Rendering",datasetList[i],"...\n")
	}

	if(saveTiff==TRUE){
		
		tiff(file = opfileName, width =3.5, height = 7.3, units = "in", res = 300)
		
		#Extract ggplots from the list
		cat('Arranging heatmaps...\n')
		n <- length(ggheatmap)
		nCol <- floor(sqrt(n))
		do.call("grid.arrange", c(ggheatmap, ncol=1))

		dev.off()
	}else{

		dev.new()
		cat('Arranging heatmaps...\n')
		#n <- length(ggheatmap)
		#nCol <- floor(sqrt(n))
		do.call("grid.arrange", c(ggheatmap, ncol=1))
		#dev.off()
	}
		
	cat('Done\n')
	
}

#Align epigenomic feature heatmaps by samples
panCanHeatmap.align<-function(synObjList , geneList , methList, sortVar, opfileName,datasetList =NULL  ){

        ggheatmap=list()
	# as.numeric(order(datasetList) ) 
	for(i in c(2:4,6:12) ){ #Sort by alphabets
	load(paste(datasetList[i],".meth.rda",sep=""))
	cat(paste(datasetList[i],".meth.rda",sep=""),"\n")
        meth= synObjList.meth$e

	expNames= colnames(synObjList[[i]]$e) [  order( synObjList[[i]]$metagene[ sortVar, ])  ]
	sampleNames=   intersect( expNames, colnames(meth)   )
	
	#cat("sampleNames:",  sampleNames)

	x <- IlluminaHumanMethylation27kSYMBOL
	mapped_probes <- mappedkeys(x)
	xxx <- unlist(as.list(x[mapped_probes]))
	rn.meth=rownames(meth)
	rn.gn=xxx[rn.meth]

	chr <-IlluminaHumanMethylation27kCHR
	mapped_chr <- mappedkeys(chr)
	chrxxx <- unlist(as.list(chr[mapped_chr]))
	rn.chr=chrxxx[rownames(meth)]

	cpg <-IlluminaHumanMethylation27kCPGCOORDINATE
	mapped_cpg <- mappedkeys(cpg)
	cpgxxx <- unlist(as.list(cpg[mapped_cpg]))
	rn.cpg=cpgxxx[rownames(meth)]

	rn=paste(rn.gn,"_",rn.chr,"_",rn.cpg,sep="")
	rownames(meth)=rn

		cat("Rendering",datasetList[i],"...\n")

	ggheatmap[[i]]<-list(
		 e=heatmap.metagene(synObjList[[i]]$e[geneList,sampleNames], geneList, datasetList[i] ),
		 m=heatmap.metagene(meth[methList,sampleNames], methList, datasetList[i],removeOutlier=FALSE,showTitle=FALSE)
		 )
               

	}	
		tiff(file = opfileName, width =3.5, height = 7.3, units = "in", res = 300)
		
		#Extract ggplots from the list
		cat('Arranging heatmaps...\n')
		n <- length(ggheatmap)
		nCol <- floor(sqrt(n))
		grid.arrange( 	ggheatmap[[2]]$e,ggheatmap[[2]]$m,
				ggheatmap[[3]]$e,ggheatmap[[3]]$m,
				ggheatmap[[4]]$e,ggheatmap[[4]]$m,
				ggheatmap[[6]]$e,ggheatmap[[6]]$m,
				ggheatmap[[7]]$e,ggheatmap[[7]]$m,
				ggheatmap[[8]]$e,ggheatmap[[8]]$m,
				ggheatmap[[9]]$e,ggheatmap[[9]]$m,
				ggheatmap[[10]]$e,ggheatmap[[10]]$m,
				ggheatmap[[11]]$e,ggheatmap[[11]]$m,
				ggheatmap[[12]]$e,ggheatmap[[12]]$m, ncol=1,
				widths = unit(rep(3,20), "in")
				)

		#do.call("grid.arrange", c(ggheatmap$e, ncol=1))

		dev.off()	
	cat('Done\n')
	
}

panCanHeatmap.align3<-function(synObjList, geneList, methList1, methList2, sortVar, opfileName, datasetList =NULL  ){

        ggheatmap=list()
	# as.numeric(order(datasetList) ) 
	for(i in c(2:4,6:12) ){ #Sort by alphabets
	load(paste(datasetList[i],".meth.rda",sep=""))
	cat(paste(datasetList[i],".meth.rda",sep=""),"\n")
        meth= synObjList.meth$e
	expNames= colnames(synObjList[[i]]$e) [  order( synObjList[[i]]$metagene[ sortVar, ])  ]
	sampleNames=   intersect( expNames, colnames(meth)   )
	
	#cat("sampleNames:",  sampleNames)

	x <- IlluminaHumanMethylation27kSYMBOL
	mapped_probes <- mappedkeys(x)
	xxx <- unlist(as.list(x[mapped_probes]))
	rn.meth=rownames(meth)
	rn.gn=xxx[rn.meth]

	chr <-IlluminaHumanMethylation27kCHR
	mapped_chr <- mappedkeys(chr)
	chrxxx <- unlist(as.list(chr[mapped_chr]))
	rn.chr=chrxxx[rownames(meth)]

	cpg <-IlluminaHumanMethylation27kCPGCOORDINATE
	mapped_cpg <- mappedkeys(cpg)
	cpgxxx <- unlist(as.list(cpg[mapped_cpg]))
	rn.cpg=cpgxxx[rownames(meth)]

	rn=paste(rn.gn,"_",rn.chr,"_",rn.cpg,sep="")
	rownames(meth)=rn

		cat("Rendering",datasetList[i],"...\n")

	ggheatmap[[i]]<-list(
		 e=heatmap.metagene(synObjList[[i]]$e[geneList,sampleNames], geneList, figTitle=paste(datasetList[i], "LYM") ,removeOutlier=FALSE,showTitle=TRUE),
		 m1=heatmap.metagene(meth[methList1,sampleNames], methList1, figTitle=paste(datasetList[i], "M+"),removeOutlier=FALSE,showTitle=TRUE),
		 m2=heatmap.metagene(meth[methList2,sampleNames], methList2, figTitle=paste(datasetList[i], "M-"),removeOutlier=FALSE,showTitle=TRUE)
		 )
               

	}	
		tiff(file = opfileName, width =3.5, height = 15, units = "in", res = 300)
		
		#Extract ggplots from the list
		cat('Arranging heatmaps...\n')
		grid.arrange( 	ggheatmap[[2]]$e,ggheatmap[[2]]$m1,ggheatmap[[2]]$m2, #Stupid way 
				ggheatmap[[3]]$e,ggheatmap[[3]]$m1,ggheatmap[[3]]$m2,
				ggheatmap[[4]]$e,ggheatmap[[4]]$m1,ggheatmap[[4]]$m2,
				ggheatmap[[6]]$e,ggheatmap[[6]]$m1,ggheatmap[[6]]$m2,
				ggheatmap[[7]]$e,ggheatmap[[7]]$m1,ggheatmap[[7]]$m2,
				ggheatmap[[8]]$e,ggheatmap[[8]]$m1,ggheatmap[[8]]$m2,
				ggheatmap[[9]]$e,ggheatmap[[9]]$m1,ggheatmap[[9]]$m2,
				ggheatmap[[10]]$e,ggheatmap[[10]]$m1,ggheatmap[[10]]$m2,
				ggheatmap[[11]]$e,ggheatmap[[11]]$m1,ggheatmap[[11]]$m2,
				ggheatmap[[12]]$e,ggheatmap[[12]]$m1,ggheatmap[[12]]$m2, ncol=1,
				widths = unit(rep(3.2,30), "in")
				)
		#do.call("grid.arrange", c(ggheatmap$e, ncol=1))
		dev.off()	
	cat('Done\n')
}

#panCanHeatmap.align3(  synObjList.f,  lymList, mpList, mmList, 'ls', "heatmap.lymxmeth.tiff",  datasetList=datasetList[,4]    )


#Render 3 aligned metagene heatmap
panCanHeatmap.gene<-function(synObjList, geneList1, geneList2, geneList3, sortVar, opfileName, datasetList =NULL, geneTitle=c("CIN", "ANG", "MES")  ){

        ggheatmap=list()
	# as.numeric(order(datasetList) ) 
	for(i in 1:12 ){ #Sort by alphabets

	expNames= colnames(synObjList[[i]]$e) [  order( synObjList[[i]]$metagene[ sortVar, ])  ]
	#sampleNames=   intersect( expNames, colnames(meth)   )
	sampleNames=    expNames
	#cat("sampleNames:",  sampleNames)

	
		cat("Rendering",datasetList[i],"...\n")

	ggheatmap[[i]]<-list(
		 e1=heatmap.metagene(synObjList[[i]]$e[geneList1,sampleNames], geneList1, figTitle=paste(datasetList[i], geneTitle[1]) ,removeOutlier=FALSE,showTitle=TRUE),
		 e2=heatmap.metagene(synObjList[[i]]$e[geneList2,sampleNames], geneList2, figTitle=paste(datasetList[i], geneTitle[2]) ,removeOutlier=FALSE,showTitle=TRUE),
		 e3=heatmap.metagene(synObjList[[i]]$e[geneList3,sampleNames], geneList3, figTitle=paste(datasetList[i], geneTitle[3]) ,removeOutlier=FALSE,showTitle=TRUE)
		 )
               

	}	
		tiff(file = opfileName, width =3.5, height = 15, units = "in", res = 300)
		
		#Extract ggplots from the list
		cat('Arranging heatmaps...\n')
		grid.arrange(	ggheatmap[[1]]$e1,ggheatmap[[1]]$e2,ggheatmap[[1]]$e3, #Stupid way
			 	ggheatmap[[2]]$e1,ggheatmap[[2]]$e2,ggheatmap[[2]]$e3, 
				ggheatmap[[3]]$e1,ggheatmap[[3]]$e2,ggheatmap[[3]]$e3,
				ggheatmap[[4]]$e1,ggheatmap[[4]]$e2,ggheatmap[[4]]$e3,
				ggheatmap[[5]]$e1,ggheatmap[[5]]$e2,ggheatmap[[5]]$e3,
				ggheatmap[[6]]$e1,ggheatmap[[6]]$e2,ggheatmap[[6]]$e3,
				ggheatmap[[7]]$e1,ggheatmap[[7]]$e2,ggheatmap[[7]]$e3,
				ggheatmap[[8]]$e1,ggheatmap[[8]]$e2,ggheatmap[[8]]$e3,
				ggheatmap[[9]]$e1,ggheatmap[[9]]$e2,ggheatmap[[9]]$e3,
				ggheatmap[[10]]$e1,ggheatmap[[10]]$e2,ggheatmap[[10]]$e3,
				ggheatmap[[11]]$e1,ggheatmap[[11]]$e2,ggheatmap[[11]]$e3,
				ggheatmap[[12]]$e1,ggheatmap[[12]]$e2,ggheatmap[[12]]$e3, ncol=1,
				widths = unit(rep(3.2,30), "in")
				)
		#do.call("grid.arrange", c(ggheatmap$e, ncol=1))
		dev.off()	
	cat('Done\n')
}

#panCanHeatmap.gene(  synObjList.f, mtList, angList, cinList, 'mt', "heatmap.angxmes.tiff",  datasetList=datasetList[,4], geneTitle=c("MES", "ANG", "CIN")    )



