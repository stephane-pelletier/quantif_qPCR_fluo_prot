
###############  main parameter  ##########################
data=1
conc=c(0,0.5,1,2,4,6,8)



############################################
############################################
##                CODE                    ##
############################################
############################################

if (data==1){
  data<-read.csv(FILE<-file.choose(),na.strings="NA", dec=",",sep=";",header=T)
  attach(data)
  DIR<-dirname(FILE)
}else{options(warn=-1)
}
titre0<-c(paste("quantifprot",format(Sys.time(), "_%a%d%b%Y"),sep=""))
DIR0<-c(paste(DIR,titre0,sep="/"))
dir.create(DIR0)
calcule<-as.matrix(data)
abso<-c(apply(calcule,2,mean))
abs<-c(abso[1:length(conc)])
f<-lm(conc~abs)
Y <- predict(f, newdata=data.frame(abs))
coeff<-c(coefficients(lm(conc~abs)))
print(lm(abs~conc))
coeff<-c(coefficients(lm(conc~abs)))
titre1<-c(paste("quantifprot",format(Sys.time(), "_%a%d%b%Y_%H.%M.%S"),".csv",sep=""))
DIR2<-c(paste(DIR0,titre1,sep="/"))
sink(DIR2,append=T,split=T)
cat("\n","\n","\n")
cat("****************************","\n")
cat("*          REPORT          *","\n")
cat("****************************","\n","\n","\n")
print(paste0("y = ", round(coeff[2],4), "*x ", round(coeff[1],4)))
print(paste0(round(coeff[2],4), "*abs", +round(coeff[1],4)))
concsample<-c(round(coeff[2],4)*abso[8:length(abso)]+round(coeff[1],4))
cat("-----------------------------","\n")
cat("\n")
concsample=matrix(signif(concsample,3),nrow=1)
data2=data[,-(1:length(conc))]
colnames(concsample)<-colnames(data2)
rownames(concsample)<-"Concentration"
cat("Quantity of protein you need to take (in ug):")
quantiteprot <- scan(file = "", integer(0), n = 1, quiet = TRUE)
cat(quantiteprot,"\n","\n")
dixmicro<-c(quantiteprot/concsample)
dixmicro=matrix(signif(dixmicro,3),nrow=1)
rownames(dixmicro)<-"Prot to take"
cat("Finale volume (in ul):")
volumefinal <- scan(file = "", integer(0), n = 1, quiet = TRUE)
cat(volumefinal,"\n","\n")
volh2o<-c(volumefinal-dixmicro)
volh2o=matrix(signif(volh2o,4),nrow=1)
rownames(volh2o)<-"Blue to add"
print(rbind(concsample,dixmicro,volh2o))
nbsample<-c(ncol(concsample))
cat("\n")
cat("Concentration are in ug/ul","\n")
cat("Volume are in ul","\n","\n")
cat("Number of sample calculate:")
cat(nbsample,"\n","\n")
cat(format(Sys.time(), "%a %d %b %Y %H:%M:%S"),"\n","\n")
cat("-----------END------------","\n")
sink()
mat2=as.matrix(rbind(concsample,dixmicro,volh2o))
titre2<-c(paste("table-quantifprot",format(Sys.time(), "_%a%d%b%Y_%H.%M.%S"),".csv",sep=""))
DIR3<-c(paste(DIR0,titre2,sep="/"))
write.table(cbind(c("Concentration","protein to take","blue to add"),mat2),DIR3,sep=";",row.names=F)
dev.new(noRStudioGD = T)
par(mar = c(5,5,2,2))
plot(abs,conc,ylab="Concentration",xlab="Absorbance",cex.lab=2.5,cex.axis=2,lwd=1.5,xlim=c(0,1),ylim=c(0,10))
lines(abs, Y,col="blue",lwd=2)
r2<-c(paste0("R2=",round(summary(lm(abs~conc))$r.squared,4)))
eq=(paste0("y = ", round(coeff[2],4), "*x +", round(coeff[1],4)))
text(median(abs),max(conc),2,labels=eq,xpd=T,srt=0,cex=1.5,col="black",font=1,pos=1)
text(median(abs),mean(conc),labels=r2,xpd=T,srt=0,cex=1.5,col="black",font=1,pos=1)
titre3<-c(paste("plot",format(Sys.time(), "_%a%d%b%Y_%H.%M.%S"),".tiff",sep=""))
DIR4<-c(paste(DIR0,titre3,sep="/"))
dev.print(device = tiff, file = DIR4, type="cairo",res=200,units="cm",width=26,height=20)

