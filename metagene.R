#Create Metagene Values from Weiyi Cheng's DreamBox7


CreateMetageneSpace <- function(ge, attractome, map, chosenProbes = NULL, maxGenes = NULL){
  if(is.null(chosenProbes)) {
  nMeta = length(attractome)
  metaSpace = matrix(NA, nrow=nMeta, ncol=ncol(ge))
  dimnames(metaSpace) = list(names(attractome), colnames(ge))
  pbs = list()
  mappedGenes = rep(NA, nrow(ge))
  names(mappedGenes) = rownames(ge)
  idx = intersect(rownames(map) , rownames(ge))
  mappedGenes[idx] = map[idx,1]
  for (i in 1:nMeta){
    #cat(i, "\n")
    #flush.console()
    a = attractome[[i]]
    if(!is.null(maxGenes)){
      if(nrow(a) > maxGenes){
        genes = a[1:maxGenes,1]
      }else{
        genes = a[, 1]
      }
    }else{
	genes = a[, 1]
    }
    il = lapply(genes, function(g){which(mappedGenes == g)})
    ill = sapply(il, length)
    goodIdx = sapply(il, function(i){ if(length(i) == 1) i})
    goodIdx = goodIdx[sapply(goodIdx, function(x){!is.null(x)})]
    numGood = sum(ill == 1)
    goodMat = NULL
    if(numGood > 0){
      goodMat = ge[unlist(goodIdx),]
    }
    badIdx = il[sapply(il, function(i){length(i) > 1 })]

    numBad = length(badIdx)
    
    badMat = NULL
    chosenIdx = NULL
    if(numBad > 0) {
      geneSum = apply(ge[unlist(il), ],2,sum)
      chosenIdx = lapply(badIdx, function(idcs){
	mis = sapply(idcs, function(i){getMI(geneSum, ge[i,])} )
	idcs[which(mis > 0.5)]
      })
      chosenIdx = chosenIdx[sapply(chosenIdx, function(x){length(x)>0})]
	#badMat = ge[chosenIdx,]
	badMat = t(sapply(chosenIdx, function(idcs){
		if(length(idcs) > 1){
			apply(ge[idcs,], 2, function(x){mean(x, na.rm=TRUE)})
		}else if(length(idcs) == 1){
			ge[idcs,]
		}else{
			rep(NA, ncol(ge))
		}
	}) )
	if(length(chosenIdx) == 0) {chosenIdx = NULL; badMat = NULL}
    }    
    pbs[[i]] = sapply(c(goodIdx, chosenIdx), names)
    metaSpace[i,] = (apply(rbind(goodMat, badMat), 2, function(x){mean(x, na.rm=TRUE)}))

    #o = sapply(genes, 
    #          function(g, ge, mappedGenes){
    #           idx = which(mappedGenes == g)
    #          if (length(idx)==1) return (ge[idx,])
    #             if (length(idx)==0) return (rep(NA, ncol(ge)))
    #             return (apply(ge[idx,], 2, function(x) mean(x, na.rm=T)))
    #           },
    #           ge = ge,
    #           mappedGenes = map[rownames(ge), "Gene.Symbol"]
    #           )
    #if(length(genes)==1){metaSpace[i,] = o}
    #else {metaSpace[i,] = apply(o, 1, function(x) mean(x, na.rm=T))}
  }
  names(pbs) = names(attractome)
  o = list(metaSpace = metaSpace, pbs = pbs)
  return (o)
  }else{
	metaSpace = t(sapply(chosenProbes, function(pb){
		pb = sapply(pb, function(p){intersect(p, rownames(ge))})
		gmat = sapply(pb, function(p, ge){
			if(length(p) > 1){
				apply(ge[p,], 2, mean)
			}else{
				ge[p,]
			}
		}, ge = ge)
		if(!is.null(dim(gmat))) {apply(gmat, 1, mean)}
		else{gmat}
	}) )
	return(metaSpace)

  }
}