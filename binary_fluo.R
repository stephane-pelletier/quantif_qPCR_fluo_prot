#fluo gfp , ko =1,wt=0
###############  main parameter  ##########################
data=1
To_be_quantify = c("TPR","Nup133")
Sample = c("cl14","s.c3","s.c5","s.d3","s.f4")
z=0 ; control="wt"   #"z" is the line you use to normalize. "control" will be the name to display in the graph
chan="chan" #field of view
fluo="fluo" #binary value if there is fluorescence or not (will be used for "z")
save.matsynth=0 ; de.matsynth=0 ; save.synth=0 ; de.synth=0
calculation=mean

#################graph####################
sep1=c(0);text1=c("");sep2=c(0);text2=c("");pos.sep1=50
pos.text1=60;col.sep1=c("black");col.text1=c("black");pos.sep2=50;pos.text2=60;col.sep2=c("black");col.text2=c("black")
w=9.5;h=7
bottom=12;left=8.5;top=4;right=9
marge=c(bottom,left,top,right)
space_axis=2;space_labelY=5;space_tick=0
space=c(space_labelY,space_axis,space_tick)
labelY="Level IF"
size_lab=3;size_main=5;size_axis=2
legend_size=2;pos.legend="topright"
col1=c("black","turquoise2","red","gray23","mediumblue","yellow","springgreen1","magenta2","orange","lightgoldenrod4","slategray")
col2=c("turquoise1","thistle1","greenyellow")




############################################
############################################
##                                        ##
##                CODE                    ##
##                                        ##
############################################
############################################


prot<-To_be_quantify
cell<-Sample


value01<-c(0,1)
stopifnot(save.matsynth%in%value01)
stopifnot(de.matsynth%in%value01)
stopifnot(save.synth%in%value01)
stopifnot(de.synth%in%value01)


if (data==1){
  data<-read.csv(FILE<-file.choose(),na.strings="NA", dec=",",sep=";",header=T)
  attach(data)
  DIR<-dirname(FILE)
  }else{options(warn=-1)
}
cell2=cell[order(cell)]
fonction0 <- function(prot, chan, fluo, data){
  df <- data[, c(prot, fluo, chan)]
  sp <- split(df, df[, chan])
  difference <- lapply(
    X = sp, 
    FUN = function(sdf) { 
      sdf[, prot] - sdf[sdf[, fluo]==2, prot] 
    }
  )
  return(difference)
}
  out <- lapply(X = cell2, FUN = function(isuffix) {
    tmp <- lapply(X = prot, FUN = function(ix) {
      fonction0(
        prot = paste(ix, isuffix, sep = "_"), 
        chan = paste(chan, isuffix, sep = "_"), 
        fluo = paste(fluo, isuffix, sep = "_"), 
        data = data
      )
    })
    names(tmp) <- prot
	return(tmp)
  })
  names(out) <- cell2
mat<-matrix(ncol=length(cell2),nrow=nrow(data))
matname<-c(apply(expand.grid(chan,cell2), 1, paste, collapse = "_"))
colnames(mat)<-matname
for (i in cell2){
  for (j in matname){
    chane<-c(apply(expand.grid(chan,i), 1, paste, collapse = "_"))
    chane2vec<-c(get(j))
    chan5<-matrix(chane2vec,ncol=1)
    colnames(chan5)<-j
    mat[,j]<-chan5[,j]
  }
}
matfluo<-matrix(ncol=length(cell2),nrow=nrow(data))
matfluoname<-c(apply(expand.grid(fluo,cell2), 1, paste, collapse = "_"))
colnames(matfluo)<-matfluoname
for (i in cell2){
  for (j in matfluoname){
    chane<-c(apply(expand.grid(fluo,i), 1, paste, collapse = "_"))
    chane2vec<-c(get(j))
    chan5<-matrix(chane2vec,ncol=1)
    colnames(chan5)<-j
    matfluo[,j]<-chan5[,j]
  }
}
matchandf<-matrix(ncol=length(cell2),nrow=nrow(data))
matname2<-c(apply(expand.grid(chan,cell2), 1, paste, collapse = "_"))
colnames(matchandf)<-matname2
for (i in cell2){
  for (j in matname2){
    chane<-c(apply(expand.grid(chan,i), 1, paste, collapse = "_"))
    chane2vec<-c(get(j))
    chan5<-matrix(chane2vec,ncol=1)
    colnames(chan5)<-j
    matchandf[,j]<-chan5[,j]
  }
}
matchandf=matchandf[,rep(1:ncol(matchandf),each=length(prot))]
matfluodf<-matrix(ncol=length(cell2),nrow=nrow(data))
colnames(matfluodf)<-matfluoname
for (i in cell2){
  for (j in matfluoname){
    chane<-c(apply(expand.grid(fluo,i), 1, paste, collapse = "_"))
    chane2vec<-c(get(j))
    chan5<-matrix(chane2vec,ncol=1)
    colnames(chan5)<-j
    matfluodf[,j]<-chan5[,j]
  }
}
matfluodf=matfluodf[,rep(1:ncol(matfluodf),each=length(prot))]
vecchandf<-as.vector(matchandf[,1:ncol(matchandf)])
vecchandf<-c(na.omit(vecchandf))
print(matfluodf)
vecfluodf<-as.vector(matfluodf[,1:ncol(matfluodf)])
print(vecfluodf)
vecfluodf<-c(na.omit(vecfluodf))
vecchan<-c(as.vector(unlist(apply(mat,2,as.list))))
vecfluo<-c(as.vector(unlist(apply(matfluo,2,as.list))))
vecchan<-c(na.omit(vecchan))
vecfluo<-c(na.omit(vecfluo))
vecfluo<-c(rep(vecfluo ,length(prot)))
vecchan<-c(rep(vecchan ,length(prot)))
matchan3df=apply(matchandf,2,sort)
matchandf2<-c(as.vector(matchandf[,1:ncol(matchandf)]))
matchandf2<-c(na.omit(matchandf2))
repdf<-c(apply(expand.grid(prot,cell2), 1, paste, collapse = "_"))
matrepdf=matrix(ncol=1,nrow=length(repdf))
rownames(matrepdf)<-repdf
for (i in repdf){
  getdf<-c(na.omit(get(i)))
  mat=matrix(length(getdf),ncol=1)
  rownames(mat)<-i
  matrepdf[i,]<-mat[i,]
}
matrepdfmat=matrix(matrepdf,nrow=length(prot))
matrepdfmat2=apply(matrepdfmat,2,sum)
vecmatrepdf<-c(as.vector(matrepdfmat[1,]))
vecmatrepdf2<-c(rep(vecmatrepdf,length(prot)))
repcell2<-c(rep(cell2,each=length(prot)))
matchandf2paste<-c(rep(repcell2,matrepdf))
repprot<-c(rep(prot,length(cell2)))
matchandf2pasteab<-c(rep(repprot,matrepdf))
matchandf2paste2<-c(paste(matchandf2paste,matchandf2pasteab ,matchandf2,sep="-"))
names(vecfluodf)<-matchandf2paste2
names(vecchandf)<-matchandf2paste2
repcell22<-c(rep(cell2,length(prot)))
repprot2<-c(rep(prot,each=length(cell2)))
repvecchanfluo<-c(rep(repprot2, vecmatrepdf2))
repcell2vecchanfluo<-c(rep(repcell22, vecmatrepdf2))
vecfluopaste<-c(paste(repvecchanfluo,repcell2vecchanfluo,vecchan,sep="-"))
names(vecfluo)<- vecfluopaste
names(vecchan)<- vecfluopaste
dffluonet=data.frame(unlist(out),vecfluodf,vecchandf)
dfchanfluo=data.frame(vecfluo,vecchan)
vecmatdf<-c(as.vector(na.omit(matrepdf[,1])))
repdfvec<-c(rep(repdf,vecmatdf))
dffluonet=data.frame(repdfvec,dffluonet)
colnames(dffluonet)<-c("name","value","fluo","chan")
print(dffluonet)
sp2<-split(dffluonet,list(dffluonet$chan,dffluonet$name))
rapport <- lapply(sp2, function(value) ((value$value / calculation(value$value[value$fluo==0],na.rm=T))*z))
rapport2=data.frame(unlist(rapport))
orderrepdf=repdf[order(repdf)]
matorderrepdf=matrix(ncol=1,nrow=length(orderrepdf))
rownames(matorderrepdf)<-orderrepdf
for (i in orderrepdf){
  getdf<-c(na.omit(get(i)))
  mat=matrix(length(getdf),ncol=1)
  rownames(mat)<-i
  matorderrepdf[i,]<-mat[i,]
}
vecmatdf2<-c(as.vector(na.omit(matorderrepdf[,1])))
repdfvec2<-c(rep(orderrepdf,vecmatdf2))
rownames(rapport2)<-NULL
dffluonorma=data.frame(repdfvec2,rapport2,dfchanfluo)
colnames(dffluonorma)<-c("name","value","fluo","chan")
repdf2<-c(apply(expand.grid(prot,cell), 1, paste, collapse = "_"))
tab <- subset(dffluonorma, dffluonorma$fluo != 2)
mat1=matrix(ncol=length(repdf2),nrow=100)
colnames(mat1)<-repdf2
for (j in repdf2){
  matrepdf1=matrix(tab$value[(tab$name==j)&(tab$fluo==1)],ncol=1)
  colnames(matrepdf1)<-j
  repna=matrix(rep(NA,nrow(mat1)-nrow(matrepdf1)),ncol=1)
    if (nrow(matrepdf1)<nrow(mat1)){
      matrepdf1<-rbind(matrepdf1,repna)
      } else {
  }
  mat1[,j]<-matrepdf1[,j]
}
mat1=mat1[apply(mat1, 1, function(y) !all(is.na(y))),]
mat0=matrix(ncol=length(repdf2),nrow=100)
colnames(mat0)<-repdf2
for (j in repdf2){
  matrepdf0=matrix(tab$value[(tab$name==j)&(tab$fluo==0)],ncol=1)
  colnames(matrepdf0)<-j
  repna=matrix(rep(NA,nrow(mat0)-nrow(matrepdf0)),ncol=1)
    if (nrow(matrepdf0)<nrow(mat0)){
      matrepdf0<-rbind(matrepdf0,repna)
    } else {
  }
  mat0[,j]<-matrepdf0[,j]
}
mat0=mat0[apply(mat0, 1, function(y) !all(is.na(y))),]
matpos=matrix(rep(1,nrow(mat1)),ncol=1)
matneg=matrix(rep(0,nrow(mat0)),ncol=1)
colnames(matpos)<-"fluo"
colnames(matneg)<-"fluo"
matsynth1=cbind(matpos,mat1)
matsynth0=cbind(matneg,mat0)
if (save.matsynth==1){
  sink("matsynth.txt",append=T,split=F)
  print(matsynth1)
  print(matsynth0)
  sink()
  } else  {
} 
print(matsynth1)
print(matsynth0)
if (de.matsynth==1){
  data.entry(matsynth1)
  data.entry(matsynth0)
  } else  {
} 
for (i in prot){
  assign(paste(i,1,sep=""),matrix(matsynth1[,grep(pattern=i,colnames(matsynth1))],ncol=length(cell2)))
}
l1<-c(mget(apply(expand.grid(prot,1), 1, paste, collapse = "")))
names(l1)<-prot
#dev.new(noRStudioGD = T,width=w, height=h)
#par(mar=c(marge))
#par(mgp=c(space))
#matplot(l1[[1]],l1[[2]],type="p",pch=1,col=col1[seq_along(col1)%in%1:length(cell)],main="GFP=1",xlab=prot[seq_along(prot)==1],ylab=prot[seq_along(prot)==2],cex.lab=size_lab,cex.main=size_main,cex.axis=size_axis)
#legend("topright",legend=cell,fill=col1[seq_along(col1)%in%1:length(cell)])
for (i in prot){
  assign(paste(i,0,sep=""),matrix(matsynth0[,grep(pattern=i,colnames(matsynth0))],ncol=length(cell2)))
}
l0<-c(mget(apply(expand.grid(prot,0), 1, paste, collapse = "")))
names(l0)<-prot
#dev.new(noRStudioGD = T,width=w, height=h)
#par(mar=c(marge))
#par(mgp=c(space))
#matplot(l0[[1]],l0[[2]],type="p",pch=1,col=col1[seq_along(col1)%in%1:length(cell)],main=control,xlab=prot[seq_along(prot)==1],ylab=prot[seq_along(prot)==2],cex.lab=size_lab,cex.main=size_main,cex.axis=size_axis)
#legend("topright",legend=cell,fill=col1[seq_along(col1)%in%1:length(cell)])
calculation1=apply(matsynth1,2,calculation,na.rm=T)
calculation0=apply(matsynth0,2,calculation,na.rm=T)
sd1=apply(matsynth1,2,sd,na.rm=T)
sd0=apply(matsynth0,2,sd,na.rm=T)
#dev.new(noRStudioGD = T,width=w, height=h)
#par(mar=c(marge))
#par(mgp=c(space))
matsynth11=matsynth1[,-1]
#boxplot(rbind(matsynth11),las=2,col=col2[seq_along(col2)%in%1:length(prot)],cex.lab=size_lab,cex.main=size_main,cex.axis=size_axis)
#legend("topright",legend=prot,fill=col2[seq_along(col2)%in%1:length(prot)])
#dev.new(noRStudioGD = T,width=w, height=h)
#par(mar=c(marge))
#par(mgp=c(space))
matsynth00=matsynth0[,-1]
#boxplot(rbind(matsynth00),las=2,col=col2[seq_along(col2)%in%1:length(prot)],cex.lab=size_lab,cex.main=size_main,cex.axis=size_axis)
#legend("topright",legend=prot,fill=col2[seq_along(col2)%in%1:length(prot)])
matsynthwt000=matrix(nrow=nrow(matsynth11),ncol=length(prot))
colnames(matsynthwt000)<-prot
repna2=matrix(rep(NA,abs(nrow(matsynth11)-nrow(matsynth00)),ncol=1))
repna22<-matrix(cbind(rep(repna2,length(prot))),ncol=length(prot))
for (i in prot){
  matsynth000=apply(matsynth00[,grep(pattern=i,colnames(matsynth00))],1,calculation,na.rm=T)
  matsynth000=matrix(matsynth000,ncol=1)
  colnames(matsynth000)<-i
#repna2=matrix(rep(NA,abs(nrow(matsynth11)-nrow(matsynth000)),ncol=1))
    if (nrow(matsynth000)<nrow(matsynth11)){
      matsynth000<-rbind(matsynth000,repna2)
    } else {
  }
    if (nrow(matsynth000)>nrow(matsynth11)){
      matsynthwt0002<-rbind(matsynthwt000,repna22)
    } else {
  }
}
matsynthwtandsample=cbind(matsynthwt000,matsynth11)
repdfwt<-c(apply(expand.grid(prot,control), 1, paste, collapse = "_"))
colnames(matsynthwtandsample)<-c(repdfwt,repdf)
#dev.new(noRStudioGD = T,width=w, height=h)
#par(mar=c(marge))
#par(mgp=c(space))
#boxplot(rbind(matsynthwtandsample),las=2,cex.lab=size_lab,cex.main=size_main,cex.axis=size_axis,col=col2[seq_along(col2)%in%1:length(prot)])
#legend("topright",legend=prot,fill=col2[seq_along(col2)%in%1:length(prot)])
matfin=rbind(calculation1,calculation0,sd1,sd0)
matcalculation=rbind(calculation0,calculation1)
matcalculation<-matcalculation [,-1]
matsd=rbind(sd0,sd1)
matsd<-matsd [,-1]
matcalculationfin=matrix(ncol=length(cell2)*2,nrow=length(prot))
colnames(matcalculationfin)<-as.vector(rbind(control,cell))
rownames(matcalculationfin)<-prot
for (i in prot){
  matcalculation2=matrix(assign(paste("mat",i,"calculation",sep=""),as.vector(matcalculation[,grep(pattern=i,colnames(matcalculation))])),nrow=1)
  colnames(matcalculation2)<-as.vector(rbind(control,cell))
  rownames(matcalculation2)<-i
  matcalculationfin[i,]<-matcalculation2[i,]
}
nmatcalculation0=matrix(apply(matsynth00,2,function(x){length(x[!is.na(x)])}),nrow=1)
nmatcalculation1=matrix(apply(matsynth11,2,function(x){length(x[!is.na(x)])}),nrow=1)
rownames(nmatcalculation0)<-"n"
rownames(nmatcalculation1)<-"n"
colnames(nmatcalculation0)<-rep(cell,each=length(prot))
colnames(nmatcalculation1)<-rep(cell,each=length(prot))
bindnmatcalculation=rbind(nmatcalculation0,nmatcalculation1)
colnames(bindnmatcalculation)<-rep(cell,each=length(prot))
vecbindnmatcalculation<-as.vector(bindnmatcalculation[,unique(colnames(bindnmatcalculation))])
allnmatcalculation0=as.vector(nmatcalculation0[,unique(colnames(nmatcalculation0))])
allnmatcalculation0<-c(sum(allnmatcalculation0))
mixnmatcalculation<-as.vector(nmatcalculation1[,unique(colnames(nmatcalculation1))])
vecbindnmatcalculation2<-c(allnmatcalculation0,mixnmatcalculation)
matvecbindnmatcalculation=matrix(vecbindnmatcalculation,nrow=1)
rownames(matvecbindnmatcalculation)<-"n"
matvecbindnmatcalculation2=matrix(vecbindnmatcalculation2,nrow=1)
rownames(matvecbindnmatcalculation2)<-"n"
matcalculationfinn=rbind(matcalculationfin,vecbindnmatcalculation)
rownames(matcalculationfinn)<-c(prot,"n")
cat("calculation","\n")
print(matcalculationfinn)
matsdfin=matrix(ncol=length(cell2)*2,nrow=length(prot))
colnames(matsdfin)<-as.vector(rbind(control,cell))
rownames(matsdfin)<-prot
for (i in prot){
  matsd2=matrix(assign(paste("mat",i,"sd",sep=""),as.vector(matsd[,grep(pattern=i,colnames(matsd))])),nrow=1)
  colnames(matsd2)<-as.vector(rbind(control,cell))
  rownames(matsd2)<-i
  matsdfin[i,]<-matsd2[i,]
}
cat("\n","\n","SD","\n")
print(matsdfin)
dev.new(noRStudioGD = T,width=w, height=h)
par(mar=c(marge))
par(mgp=c(space))
bar<-c(barplot(matcalculationfin,beside=T,col=col2[seq_along(col2)%in%1:length(prot)],ylab=labelY,cex.lab=size_lab,cex.axis=size_axis,cex.main=size_main,cex.names=2,las=2,ylim=c(0,max(matcalculationfin+matsdfin))))
legend(max(bar)+0.5,1,legend=rownames(matcalculationfin),fill=col2[seq_along(col2)%in%1:length(prot)],xpd=T,bty="n",cex=legend_size)
segments (bar,matcalculationfin-matsdfin,bar,matcalculationfin+matsdfin)
segments (bar-0.1,matcalculationfin-matsdfin,bar+0.1,matcalculationfin-matsdfin)
segments (bar-0.1,matcalculationfin+matsdfin,bar+0.1,matcalculationfin+matsdfin)
bar3=matrix(bar,nrow=2)
text(apply(bar3,2,calculation), (-70*max(matcalculationfin))/100,labels=vecbindnmatcalculation,xpd=T,pos=1,cex=ifelse(ncol(matcalculationfin)>5,1,1.5))
text(-2, (-70*max(matcalculationfin))/100,labels="n = ",xpd=T,pos=1,cex=1.5)
barsep1=matrix(bar,nrow=nrow(matcalculationfin))
name<-c(rep(1:length(sep1),sep1))
colnames(barsep1)<-name
barsep1=apply(barsep1,2,calculation)
barsep1=as.vector(barsep1)
names(barsep1)<-name
matsep1=matrix(nrow=length(unique(name)),ncol=1)
for (i in name[name>1]){
  test<-c(median(c(max((barsep1)[names(barsep1)%in%unique(i-1)]),min((barsep1)[names(barsep1)%in%unique(i)]))))
  matsep1[i,]<-test
}
matvecsep1=as.vector(na.omit(matsep1))
mat0sep1=rep(0,length(matvecsep1))
pourcsep1=(-pos.sep1*max(matcalculationfin))/100
pourcrepsep1=rep(pourcsep1,length(matvecsep1))
segments(matvecsep1,mat0sep1,matvecsep1,pourcrepsep1,xpd=T,col=col.sep1)
mat2text1=matrix(nrow=length(unique(name)),ncol=1)
for (i in name){
  test<-c(calculation(barsep1[names(barsep1)==i]))
  mat2text1[i,]<-test
}
matvec2text1=as.vector(na.omit(mat2text1))
pourc2text1=(-pos.text1*max(matcalculationfin))/100
pourcrep2text1=rep(pourc2text1,length(matvecsep1))
if (length(sep1) >1){
  text(matvec2text1,pourcrep2text1,labels=text1,xpd=T,cex=1.5,col=col.text1)
  }else{
}
matcalculationfinwt=matcalculationfin[,grep(pattern=control,colnames(matcalculationfin))]
matcalculationfinsample=matcalculationfin[,-grep(pattern=control,colnames(matcalculationfin))]
matsdfinwt=matsdfin[,grep(pattern=control,colnames(matsdfin))]
matsdfinsample=matsdfin[,-grep(pattern=control,colnames(matsdfin))]
matcalculationwtcalculation=apply(matcalculationfinwt,1,calculation,na.rm=T)
print(matcalculationwtcalculation)
print(matcalculationfinsample)
matsdwtcalculation=apply(matsdfinwt,1,calculation,na.rm=T)
matcalculationfin2=cbind(matcalculationwtcalculation,matcalculationfinsample)
print(matcalculationfin2)
matsdfin2=cbind(matsdwtcalculation,matsdfinsample)
colnames(matcalculationfin2)<-c(control,cell)
colnames(matsdfin2)<-c(control,cell)
matcalculationfin2n=rbind(matcalculationfin2,vecbindnmatcalculation2)
rownames(matcalculationfin2n)<-c(prot,"n")
if (save.synth==1){
  sink("synth.txt",append=T,split=F)
  print(matcalculationfinn)
  print(matsdfin)
  print(matcalculationfin2n)
  print(matsdfin2)
  sink()
  } else  {
} 
if (de.synth==1){
  data.entry(matcalculationfin)
  data.entry(matsdfin)
  data.entry(matcalculationfin2)
  data.entry(matsdfin2)
  } else  {
} 
cat("\n","\n","calculation : control bind","\n")
print(matcalculationfin2n)
cat("\n","\n","SD : control bind","\n")
print(matsdfin2)
dev.new(noRStudioGD = T,width=w, height=h)
#par(mar=c((800*max(matcalculationfin2))/100,8.5,4,9))
par(mar=c(marge))
par(mgp=c(space))
bar2<-c(barplot(matcalculationfin2,beside=T,col=col2[seq_along(col2)%in%1:length(prot)],ylab=labelY,cex.lab=size_lab,cex.axis=size_axis,cex.main=size_main,cex.names=2,las=2,ylim=c(0,max(matcalculationfin2+matsdfin2))))
legend(max(bar2)+0.5,1,legend=rownames(matcalculationfin2),fill=col2[seq_along(col2)%in%1:length(prot)],xpd=T,bty="n",cex=legend_size)
segments (bar2,matcalculationfin2-matsdfin2,bar2,matcalculationfin2+matsdfin2)
segments (bar2-0.1,matcalculationfin2-matsdfin2,bar2+0.1,matcalculationfin2-matsdfin2)
segments (bar2-0.1,matcalculationfin2+matsdfin2,bar2+0.1,matcalculationfin2+matsdfin2)
bar4=matrix(bar2,nrow=2)
text(apply(bar4,2,calculation), ifelse(length(sep2)>1, (-70*max(matcalculationfin2))/100, (-60*max(matcalculationfin2))/100),labels=vecbindnmatcalculation2,xpd=T,pos=1,cex=ifelse(ncol(matcalculationfin2)>5,1.5,2))
text(-2,ifelse(length(sep2)>1, (-70*max(matcalculationfin2))/100, (-60*max(matcalculationfin2))/100),labels="n = ",xpd=T,pos=1,cex=2)
barsep2=matrix(bar2,nrow=nrow(matcalculationfin2))
name<-c(rep(1:length(sep2),sep2))
colnames(barsep2)<-name
barsep2=apply(barsep2,2,calculation)
barsep2=as.vector(barsep2)
names(barsep2)<-name
matsep2=matrix(nrow=length(unique(name)),ncol=1)
for (i in name[name>1]){
  test<-c(median(c(max((barsep2)[names(barsep2)%in%unique(i-1)]),min((barsep2)[names(barsep2)%in%unique(i)]))))
  matsep2[i,]<-test
}
matvecsep2=as.vector(na.omit(matsep2))
mat0sep2=rep(0,length(matvecsep2))
pourcsep2=(-pos.sep2*max(matcalculationfin))/100
pourcrepsep2=rep(pourcsep2,length(matvecsep2))
segments(matvecsep2,mat0sep2,matvecsep2,pourcrepsep2,xpd=T,col=col.sep2)
mat2text2=matrix(nrow=length(unique(name)),ncol=1)
for (i in name){
  test<-c(calculation(barsep2[names(barsep2)==i]))
  mat2text2[i,]<-test
}
matvec2text2=as.vector(na.omit(mat2text2))
pourc2text2=(-pos.text2*max(matcalculationfin))/100
pourcrep2text2=rep(pourc2text2,length(matvecsep2))
if (length(sep2) >1){
  text(matvec2text2,pourcrep2text2,labels=text2,xpd=T,cex=1.5,col=col.text2)
  }else{
}
for (i in prot){
  matsynth002<-matsynth00[,grep(i,colnames(matsynth00))]
  vecmatsynth002<-na.omit(as.vector(matsynth002))
  matsynth112<-matsynth11[,grep(i,colnames(matsynth11))]
  newna<-rep(NA,length(vecmatsynth002)-nrow(matsynth112))
  newna2<-rep(newna,ncol(matsynth112))
  matna2<-matrix(newna2,ncol=ncol(matsynth112))
  matsynth112<-rbind(matsynth112,matna2)
  matsynth112<-cbind(vecmatsynth002,matsynth112)
  dev.new(noRStudioGD = T,width=w, height=h)
  par(mar=c(marge))
  par(mgp=c(space))
  boxplot(matsynth112,main=i,names=c(control,cell),las=2)
  assign(paste("matsynth114plot",i,sep=""),matsynth112[,grep(i,colnames(matsynth112))])
  assign(paste("matsynth114plot",i,sep=""),cbind(vecmatsynth002,get(paste("matsynth114plot",i,sep=""))))
}
l14plot<-c(mget(apply(expand.grid("matsynth114plot",prot), 1, paste, collapse = "")))
print(l14plot)
dev.new(noRStudioGD = T,width=w, height=h)
par(mar=c(marge))
par(mgp=c(space))
matplot(l14plot[[1]],l14plot[[2]],type="p",pch=1,col=col1[seq_along(col1)%in%1:length(c(control,cell))],main="",xlab=prot[seq_along(prot)==1],ylab=prot[seq_along(prot)==2],cex.lab=size_lab,cex.main=size_main,cex.axis=size_axis)
legend("topright",legend=c(control,cell),fill=col1[seq_along(col1)%in%1:length(c(control,cell))])

