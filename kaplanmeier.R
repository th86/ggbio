#Package:	Biological Figure Rendering Package using ggPlot	
#Author:	Tai-Hsien Ou Yang
#Description:	Render a KaplanMeier plot



# source("http://bioconductor.org/biocLite.R")
# biocLite("survcomp")

library(survcomp)




kmc.plot<-function( surv,feature,title){

X<-cbind(surv[,1],surv[,2], feature> median(feature) )
colnames(X)=c("time","status", "x")
fit <- survfit(Surv(X[,1], X[,2]) ~ x, data = data.frame(X))
pval<-summary(coxph(Surv(X[,1], X[,2]) ~ x, data.frame(X)))$logtest[3]

plot(fit, col=c("20","18"),lwd=1:1, main=title) 
#legend(150, 1, c("GENE-", "GENE+"), lwd=1:1, col=c("20","18"),cex=1) 
#text(160,0.1,"P",font=3,cex=0.8)
text(max(surv[,1])*0.6,0.9,paste("P value=",round(pval,4)),cex=0.8)

ci=concordance.index(feature, surv[,1], surv[,2])$c.index
text(max(surv[,1])*0.6,0.8,paste("CI=",round(ci,4)),cex=0.8)

return( fit  )
}


#kmc.plot(syn$surv,syn$metagene["mitotic",],"Red>median")
