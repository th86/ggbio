#Package:	Biological Figure Rendering Package using ggPlot	
#Author:	Tai-Hsien Ou Yang
#Description:	Functions for data loading and preprocessing


#source('http://depot.sagebase.org/CRAN.R')
#pkgInstall("synapseClient")

#source("http://www.bioconductor.org/biocLite.R") 
#biocLite("impute") 


library('impute')
library('synapseClient')
#synapseLogin()

source('metagene.R')
#lung
#cli='syn1446125'
#mrna='syn418033'
#syn=loadsynapse(mrna,cli)

#synapseQuery("select id from entity where parentId=='syn300013' and acronym=='BRCA'")
#qu=synapseQuery("select name, id from entity where parentId=='syn395688'")


datasetList=matrix(NA,12,4)

#PanCanFreeze4.7  #LUSC 	#Read	    #GBM	   #BLCA	#UCEC		#COAD	     #OV	  #LAML		#HNSC		#LUAD	   #BRCA	 #KIRC
datasetList[,1]=c('syn1446127', 'syn1446153', 'syn1446090', 'syn1571519', 'syn1446169', 'syn1446080', 'syn1446137', 'syn1571544','syn1571433', 'syn1571481', 'syn1446067','syn1446103' ) #Clinical
datasetList[,2]=c('syn418033','syn1446276', 'syn1446214', 'syn1571504', 'syn1446289', 'syn1446197', 'syn1446264', 'syn1681084','syn1571420', 'syn1571468', 'syn417812', 'syn417925' ) #mRNA
datasetList[,3]=c('syn415758','syn416194', 'syn412284', 'syn1571509', 'syn416204', 'syn411993',	    'syn415945',  'syn1571536', 'syn1571424', 'syn1571458', 'syn411485','syn1445982') #methylation
datasetList[,4]=c('LUSC',     'READ',	    'GBM',	  'BLCA',	'UCEC',	      'COAD',	    'OV',	  'LAML',      'HNSC',	     'LUAD',	   'BRCA',	'KIRC') #Cancer Type

datasetList=datasetList[order(datasetList[,4]),]


#No Intersect Dataset
#synObjList.f=list()
#for(i in 1:12)
#synObjList.f[[i]]=loadsynapse(datasetList[i,2],datasetList[i,1],datasetList[i,4],intersectClnc=FALSE)
#save(synObjList.f, file="panCan12v47.notintersect.rda")
#ls()

#Intersected Dataset and KM curve
#dev.new()
#synObjList=list()
#par(mfrow=c(3,4)) #ver,hor
#for(i in 1:12){
#synObjList[[i]]=loadsynapse(datasetList[i,2],datasetList[i,1],datasetList[i,4])
#kmc.plot(synObjList[[i]]$sur,synObjList[[i]]$metagene["mitotic",],datasetList[i,4])
#}
#save(synObjList, file="panCan12v47.rda")
#ls()


#Generating high resolution KM curve
#tiff(file = "km.tiff", width =7.3, height = 7.3, units = "in", res = 300)
#par(
#mar = c(2,2,2,1),
#mgp = c(1, 0.4, 0),
#mfrow = c(3,4)) #ver,hor
#for(i in 1:12){
#kmc.plot(synObjList[[i]]$sur,synObjList[[i]]$metagene["mitotic",],datasetList[i,3])
#}
#dev.off()
#ls()

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

load.exp<-function (file, sep = "\t") 
{
    line = readLines(file)
    tokens = strsplit(line[1], "\t")[[1]]
    n = length(tokens) - 1
    m = length(line) - 1
    mset = matrix(NA, m, n)
    cname = tokens[2:(n + 1)]
    rname = rep(NA, m)
    b = txtProgressBar(style = 3)
    for (i in 1:m) {
        tokens = strsplit(line[i + 1], "\t")[[1]]
        tokens[tokens == "NA"] = NA
        rname[i] = tokens[1]
        mset[i, ] = as.numeric(tokens[2:(n + 1)])
        if (i%%100 == 0) {
            setTxtProgressBar(b, i/m)
        }
    }
    cat("\n")
    colnames(mset) <- cname
    rownames(mset) <- rname
    mset
}

loadsynapse <-function( mRNA.ID, cli.ID ,anno.ID , intersectClnc=TRUE) {
cat("Downloading Synapse Entities...\n")
	mrna = loadEntity(mRNA.ID) 
	cli = loadEntity(cli.ID) 

	cat("Loading mRNA Datasets...\n")
	e = load.exp(file.path(mrna$cacheDir, mrna$files[[1]][1]))


	cat("Imputing Expression data...\n")

	gs<-getGeneSymbols(rownames(e))
	#KNN impute	
	e.impute<-impute.knn(e)$data
	rownames(e.impute)=gs
	colnames(e)=substr(colnames(e),1,12)
	colnames(e.impute)=colnames(e)
	e.impute=log2(e.impute+1)

	#mean imputation
	#for(i in 1:nrow(ge)){
 	#	ge[i, is.na(ge[i,] )] <- mean(ge[i,], na.rm = TRUE) 
 	#	#ge[i, is.nan(ge[i,] )] <- mean(ge[i,], na.rm = TRUE) 
	#}

	

	cat("Loading Clinical Data...\n")
	clinical.raw =read.delim( file.path(cli$cacheDir, cli$files[[1]][1]) ,header=T)
	rownames(clinical.raw)=clinical.raw[,1]
		
	surv.days=as.numeric(clinical.raw[,"days_to_death"]) #days_to_last_known_alive
	surv.status=as.character(clinical.raw[,"vital_status"])
	surv.status=gsub('DECEASED',1,surv.status)
	surv.status=as.numeric(gsub('LIVING',0,surv.status))
	
	if(length(clinical.raw[which(surv.status==0),"days_to_last_known_alive"])!=0){
	surv.days[which(surv.status==0)]=as.numeric(clinical.raw[which(surv.status==0),"days_to_last_known_alive"])
	}

	surv=cbind(surv.days,surv.status ) 
	surv.complete=surv[complete.cases(surv),]
	rownames(surv.complete)=rownames(clinical.raw)[complete.cases(surv)]
	colnames(surv.complete)=c("days","status")
	
	cat("Creating Metagenes...\n")
	load("attractome.minimalist.rda")
	rn.ge=as.matrix(gs)
	rownames(e)=rn.ge
	dim(rn.ge)=c(length(rn.ge),1)
	colnames(rn.ge)="Gene.Symbol"
	rownames(rn.ge)=rn.ge
	metagene=CreateMetageneSpace(e,attractome.minimalist, rn.ge)$metaSpace
	
	cat("Summarizing...\n")	
	intersect.set=intersect(colnames(e),rownames(surv.complete))
	cat(length(intersect.set), "Samples\n")	
	
	if(intersectClnc==TRUE){
	synapseObj=list(e=e.impute[,intersect.set],metagene=metagene[,intersect.set],sur=surv.complete[intersect.set,], type=anno.ID    )
	}else{
	synapseObj=list(e=e.impute,metagene=metagene,sur=surv.complete[intersect.set,], type=anno.ID)
	}


	return(synapseObj)
}










