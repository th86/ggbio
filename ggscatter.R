#Package:	Biological Figure Rendering Package using ggPlot	
#Author:	Tai-Hsien Ou Yang
#Description:	Render scatterplots using ggplot

#dat <- data.frame(x=1:20,y=rnorm(20,0,10), v=20:1)
#qplot(x,y,colour=v,data=dat,size=4)+scale_colour_gradient(low="blue", high="red") 
#p <- ggplot(data = dat, aes(x = x, y = y, colour = v)) + geom_point(size = 2) + scale_colour_gradient(low="blue", high="red") 


#setwd("/home/tai-hsien/Desktop/")
#source("mypancan12.R")


library(gridExtra)
library(reshape2)
library(ggplot2)
library(scales)

library("IlluminaHumanMethylation27k.db")

load("panCan12v47.notintersect.rda")
load("allattractome.rda")


cinList=c("TPX2","NCAPG","KIFC1","KIF4A","NCAPH","HJURP","CENPA","CCNA2","KIF2C","BUB1B")
lymList=c("CD53","SASH3","LAPTM5","PTPRC","CD3E","LAIR1","NCKAP1L","SPI1","HAVCR2","CD2")
mtList=c("COL3A1","COL1A2","AEBP1","COL5A1","COL1A1","THBS2","COL6A3","FBN1","COL5A2","EMILIN1")
geneList<-allattractome$rnaseq$'AHSA2|130872'[,1]
ashaList=getGeneSymbols(geneList)
geneList<-allattractome$rnaseq$'ANG'[,1]
angList=getGeneSymbols(geneList)


mpList=getGeneSymbols(allattractome$meth$'NA_NA_29848860'[,1])
mmList=getGeneSymbols(allattractome$meth$'BIN2_12_50003941'[,1])
mmList[2]="PCED1B_12_45896487"
names(mmList[2])=mmList[2]



ggscatter.gene.meth<-function( synObjList, geneName, methName1, methName2, datasetList, opfileName, axislabel=NULL, ilm27k=FALSE   ){

  scatter.gg=list()

	# as.numeric(order(datasetList) ) 

	if(ilm27k==TRUE){
		DatasetAvail=c(2:4,6:12)
	}else{
		DatasetAvail=1:12
	}

	for(i in DatasetAvail  ){ #Sort by alphabets   
       
	load(paste(datasetList[i],".meth.rda",sep=""))
	cat(paste(datasetList[i],".meth.rda",sep=""),"\n")
        meth= synObjList.meth$e
	expNames= colnames(synObjList[[i]]$e) #[  order( synObjList[[i]]$metagene[ sortVar, ])  ]
	sampleNames=   intersect( expNames, colnames(meth)   )



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

	if(length(geneName)>1){
		ge=colSums(synObjList[[i]]$e[geneName,sampleNames])/length(geneName)
		m1=colSums(meth[methName1,sampleNames])/length(methName1)
		m2=colSums(meth[methName2,sampleNames])/length(methName2)
	}else{
		ge=scale(synObjList[[i]]$e[geneName,sampleNames])
		m1=scale(meth[methName1,sampleNames])
		m2=scale(meth[methName2,sampleNames])
	}

	sc.data=data.frame( 	ge = ge-median(ge),
				m1 = (m1-min(m1) )/( max(m1)-min(m1)) ,
				m2 = (m2-min(m2) )/( max(m2)-min(m2)) )
	#print(sc.data)
	medexp=median(ge)
	cat(medexp,"\n")





	 p<-ggplot(data = sc.data, aes(x = m1, y = m2, colour = ge)) +
	 			geom_point(size = 1) +
	 			#scale_colour_gradient(low="blue", high="red", midpoint= median(sc.data$ge)) +
				#scale_x_continuous(limits=c(-5,3))+ 
				#scale_y_continuous(limits=c(-3,3))+
	 			theme_bw() +
  		 		theme(
    		 		panel.grid.major = element_blank(),
    		 		panel.grid.minor = element_blank(),
    		 		panel.background = element_blank(),
				legend.position = "right",
				legend.key.size = unit(0.5, "lines"),
				legend.title = element_text( size=6, face = "plain"),
				legend.text = element_text( size=4, face = "plain")
		 		#panel.border = element_blank()
				) + xlab( axislabel[1] ) +ylab( axislabel[2] ) +
				ggtitle(datasetList[i])+
				theme(
				plot.title=element_text(size=8),
				plot.margin= unit(c(0,1,0,1), "lines"), #Top,left,down,right
				axis.title=element_text(size=6),
				axis.text=element_text(size=4)
				)
					
	p<-p+coord_fixed(ratio = 1)
		
	scatter.gg[[i]]<-p +scale_color_gradient2(axislabel[3],low="blue", high="red", mid ="gray", midpoint=0 )	#Don't use "="
	}

	cat('Arranging Scatterplots...\n')

	tiff(file = opfileName, width =7.3, height = 7.3, units = "in", res = 300)


	if(ilm27k==TRUE ){
		grid.arrange(
			scatter.gg[[2]],
			scatter.gg[[3]],
			scatter.gg[[4]],
			scatter.gg[[6]],
			scatter.gg[[7]],
			scatter.gg[[8]],
			scatter.gg[[9]],
			scatter.gg[[10]],
			scatter.gg[[11]],
			scatter.gg[[12]],
			ncol=3
			)
	}else{
		do.call("grid.arrange", c(scatter.gg, ncol=3))
	}		
		
	dev.off()
	cat('Done\n')

}


#ggscatter.gene.meth(  synObjList.f, lymList, mmList, mpList, datasetList[,4], "scatter.12.mmxmpxlym.tiff"    )



ggscatter.gene<-function( synObjList, geneNames, flabel=NULL , datasetList, opfileName    ){

  scatter.gg=list()

	# as.numeric(order(datasetList) ) 
	for(i in 1:12  ){ #Sort by alphabets   
       
	sampleNames= colnames(synObjList[[i]]$e) 

	cat("Rendering",datasetList[i],"...\n")

	if(length(geneNames[[1]])>1){
		e1=colSums(synObjList[[i]]$e[geneNames[[1]],sampleNames])/length(geneNames[[1]])
		e2=colSums(synObjList[[i]]$e[geneNames[[2]],sampleNames])/length(geneNames[[2]])
		e3=colSums(synObjList[[i]]$e[geneNames[[3]],sampleNames])/length(geneNames[[3]])
	}else{
		e1=scale(synObjList[[i]]$e[geneNames[[1]],sampleNames])
		e2=scale(synObjList[[i]]$e[geneNames[[2]],sampleNames])
		e3=scale(synObjList[[i]]$e[geneNames[[3]],sampleNames])
	}

	sc.data=data.frame( 	me1 = e1,
				me2 = e2,
				me3 = e3 )

	

	 p<-ggplot(data = sc.data, aes(x = me1, y = me2, colour = me3)) +
	 			geom_point(size = 1) +
				scale_color_gradient2("LYM",low="blue", high="red", mid ="gray", midpoint= median(sc.data$me3)  )+
	 			#scale_colour_gradient(low="blue", high="red", midpoint= median(sc.data$me3)) +
				#scale_x_continuous(limits=c(-5,3))+ 
				#scale_y_continuous(limits=c(-3,3))+
				scale_fill_discrete(name= flabel[3]   )+
	 			theme_bw() +
  		 		theme(
    		 		panel.grid.major = element_blank(),
    		 		panel.grid.minor = element_blank(),
    		 		panel.background = element_blank(),
				legend.position = "right",
				legend.key.size = unit(0.5, "lines"),
				legend.title = element_text( size=6, face = "plain"),
				legend.text = element_text( size=6, face = "plain")
		 		#panel.border = element_blank()
				) + xlab( flabel[1] ) +ylab( flabel[2] ) +
				ggtitle(datasetList[i])+
				theme(
				plot.title=element_text(size=8),
				plot.margin= unit(c(0,1,0,1), "lines"), #Top,left,down,right
				axis.title=element_text(size=6),
				axis.text=element_text(size=4)
				)
					

	scatter.gg[[i]]<-p	#Don't use "="
	p
	}

	cat('Arranging Scatterplots...\n')

	tiff(file = opfileName, width =7.3, height = 7.3, units = "in", res = 300)


	do.call("grid.arrange", c(scatter.gg, ncol=3))
				
		
	dev.off()
	cat('Done\n')

	return(scatter.gg) 
}


#ggscatter.gene(  synObjList.f, list( cinList, angList, lymList), flabel=c("CIN","ANG","LYM") , datasetList[,4], "scatter.12.3g.tiff"    )





