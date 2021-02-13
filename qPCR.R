
###############  main parameter  ##########################
data=1
primer=c("rp49","pdm3_2","pdm3_3","vvl_1")
condition=c("wt1","wt2","gd","kk1","kk2","pdm3_a","pdm3_b")   #cell line, wt, ko, mutant...
replicate=3
HousekeepingGene="rp49"
timepoint=0

####################  Optional parameters  ##################
Standard=0
dna=c(log10(1/2),log10(1/8),log10(1/16))
rev=0
Robustness=0
de=0
save=0   ############ save final result in .txt

###############  graph  #########################
w=9.5;h=7
bottom=4;left=8.5;top=4;right=9;marge=c(bottom,left,top,right)
space_axis=2;space_labelY=6;space_tick=0
title="";labelX="";labelY="";
col=c(1:100)   ####### for timeline
colbar=c("black")   ####### not for timeline
nb=c(0:5)   ####### for timeline
size_lab=3;size_main=5;size_axis=2;cex.names=2
legend_size=2;xleg=max(nb)+0.15

############## optional parameters for data ##############
CPnames = "Cp"
PRIMERnames = "gene"
LINEnames = "name"

############################################
############################################
##                CODE                    ##
############################################
############################################

cond<-condition
rq<-HousekeepingGene
repli<-replicate
space=c(space_labelY,space_axis,space_tick)

value01<-c(0,1)
stopifnot(Standard%in%value01)
stopifnot(Robustness%in%value01)
stopifnot(rev%in%value01)
stopifnot(timepoint%in%value01)
stopifnot(de%in%value01)
stopifnot(save%in%value01)

#############select data################
if (data==1){
  data<-read.csv(FILE<-file.choose(),na.strings="NA", dec=",",sep=";",header=T)
  attach(data)
  DIR<-dirname(FILE)
}else{options(warn=-1)
}

colnames(data)[grepl(CPnames,colnames(data))]<-"Cp"
colnames(data)[grepl(PRIMERnames,colnames(data))]<-"gene"
colnames(data)[grepl(LINEnames,colnames(data))]<-"name"
################ordering data##############
Cp<-data$Cp
name<-data$name
gene<-data$gene
if(Standard==1){
  type<-data$type
}  else  {
} 
a<-c(Cp[order(name,gene)])
b<-c(as.vector(gene[order(name,gene)]))
d<-c(as.vector(name[order(name,gene)]))
e=data.frame(b,d,a)

##################### check if the number of replicat is correct ################
dfs<-split(e,list(e$b,e$d))

te<-c(unique(paste(b,d,sep=".")))
for (i in te){
  if (nrow(dfs[[i]])<repli){
    dfs[[i]]<-(rbind(unname(dfs[[i]]),c(unlist(strsplit(i,".",fixed=T)),lapply(dfs[[i]]["a"],mean))))
    names(dfs[[i]])<-colnames(e)
  }else{}
}
df_new<-do.call(rbind.data.frame, dfs)
rownames(df_new)<-NULL
e<-df_new
###################################################################################


mattest=matrix(e$d,ncol=length(e$d)/repli,byrow=F)
d2<-c(as.vector(mattest[1,]))
replicat=matrix(e$a,ncol=length(d2),byrow=F)
colnames(replicat)<-d2
mattestb=matrix(b,ncol=length(e$d)/repli,byrow=F)
b2<-c(as.vector(mattestb[1,]))


#nom<-as.character(substitute(primer)) 
#nom2=nom[-1]
#nom3=nom2[order(nom2)]
nom3=primer
#################control : robustness and standard#################
if (Robustness ==1){
  if (Standard==1){
    repStandard<-c(rep("Standard",length(unique(name[type=="Standard"]))))
    repSample<-c(rep("Sample",length(unique(name[type=="Sample"]))))
    mixrep<-c(repStandard,repSample)
    f2<-c(rep(mixrep,length(unique(gene))))
     }  else  {
    f2<-c(rep("Sample",length(unique(name))))
  } 
  mean_replicat=apply(replicat,2,mean)
  mean_replicat_max=apply(replicat,2,max)
  mean_replicat_min=apply(replicat,2,min)
  sous_mat<-c(mean_replicat_max-mean_replicat_min)
  sous_mat<-abs(sous_mat)
  sd_replicat=apply(replicat,2,sd)
  e202=data.frame(mean_replicat,sd_replicat,colnames(replicat),b2,f2,sous_mat)
  par(mar = c(6,2.5,5.5,1))
  par(mfrow=c(2,3))
  names(nom3)<-seq_along(nom3)
  for (n in names(nom3)){
    n2<-as.numeric(n)
    vecnamearg=e202$colnames.replicat.[e202$colnames.replicat.%in%nom3[names(nom3)==n]]
    bar<-c(barplot(as.vector(e202)[,"mean_replicat"][e202$b%in%nom3[names(nom3)==n]],col="lightblue",las=2,names=levels(vecnamearg),cex.names=1))
    title(main=nom3[names(nom3)==n],cex.main=2,outer=F,line=0.5)
    segments (bar,as.vector(e202)[,"mean_replicat"][e202$b%in%nom3[names(nom3)==n]]-as.vector(e202)[,"sd_replicat"][e202$b%in%nom3[names(nom3)==n]],bar,as.vector(e202)[,"mean_replicat"][e202$b%in%nom3[names(nom3)==n]]+e202[,"sd_replicat"][e202$b%in%nom3[names(nom3)==n]])
    segments (bar-0.1,as.vector(e202)[,"mean_replicat"][e202$b%in%nom3[names(nom3)==n]]-as.vector(e202)[,"sd_replicat"][e202$b%in%nom3[names(nom3)==n]],bar+0.1,as.vector(e202)[,"mean_replicat"][e202$b%in%nom3[names(nom3)==n]]-as.vector(e202)[,"sd_replicat"][e202$b%in%nom3[names(nom3)==n]])
    segments (bar-0.1,as.vector(e202)[,"mean_replicat"][e202$b%in%nom3[names(nom3)==n]]+as.vector(e202)[,"sd_replicat"][e202$b%in%nom3[names(nom3)==n]],bar+0.1,as.vector(e202)[,"mean_replicat"][e202$b%in%nom3[names(nom3)==n]]+as.vector(e202)[,"sd_replicat"][e202$b%in%nom3[names(nom3)==n]])
    text(bar,(90*min(as.vector(e202)[,"mean_replicat"][e202$b%in%nom3[names(nom3)==n]],na.rm=T))/100,labels=ifelse(as.vector(e202)[,"sous_mat"][e202$b%in%nom3[names(nom3)==n]]<0.25,"*",ifelse(as.vector(e202)[,"sous_mat"][e202$b%in%nom3[names(nom3)==n]]>=0.25,"*",ifelse(as.vector(e202)[,"sous_mat"][e202$b%in%nom3[names(nom3)==n]]>0.5,"*",""))),col=ifelse(as.vector(e202)[,"sous_mat"][e202$b%in%nom3[names(nom3)==n]]<0.25,"blue",ifelse((as.vector(e202)[,"sous_mat"][e202$b%in%nom3[names(nom3)==n]]>=0.25&as.vector(e202)[,"sous_mat"][e202$b%in%nom3[names(nom3)==n]]<0.5),"green",ifelse((as.vector(e202)[,"sous_mat"][e202$b%in%nom3[names(nom3)==n]]>=0.5&as.vector(e202)[,"sous_mat"][e202$b%in%nom3[names(nom3)==n]]<1),"yellow",ifelse(as.vector(e202)[,"sous_mat"][e202$b%in%nom3[names(nom3)==n]]>=1,"red","black")))),xpd=T,cex=2)
    numberdevice<-c((n2%/%6)+1[(n2 %% 6 != 0)])
    title(main=ifelse(length(nom3)<=6 ,"Robustness",ifelse(length(nom3)>6 & n2%%6==1,paste("Robustness (",numberdevice,"/",round(length(nom3)/6),")",sep=""),"")),outer=T,line=-2,cex.main=3,font.main=2)
    if (n2 %% 6 == 0 & n2 != 1 & n2!= length(nom3)) {dev.new(noRStudioGD = T);par(mfrow=c(2,3))}
  } 
} else  {
} 



if (Standard==1){
  repStandard<-c(rep("Standard",length(unique(name[type=="Standard"]))))
  repSample<-c(rep("Sample",length(unique(name[type=="Sample"]))))
  mixrep<-c(repStandard,repSample)
  f2<-c(rep(mixrep,length(unique(gene))))
  mean_replicat=apply(replicat,2,mean)
  e2=data.frame(mean_replicat,colnames(replicat),b2,f2)
  dna<-c(sort(dna))
  if (rev==1){
    dna<-c(rev(sort(dna)))
  } else  {
  } 
  print(dna)
  valStandard=as.vector(unique(name[type=="Standard"]))
  xval<-c(e2$mean_replicat[(e2$b2%in%primer)&(d2%in%valStandard)])
  xvalmat=matrix(xval,ncol=length(xval)/length(dna),byrow=T)
  nameStandard<-c(as.vector(name[type=="Standard"]))
  uname<-c(unique(nameStandard))
  colnames(xvalmat)<-nom3
  rownames(xvalmat)<-uname
  dev.new(noRStudioGD = T)
  par(mfrow=c(2,3))
  for (j in 1:ncol(xvalmat)){
    plot(dna,xvalmat[,j],cex.lab=2.5,cex.axis=2,lwd=1.5,xlim=c(-2,0),ylim=c(min(xvalmat[,j],na.rm=T)-1,max(xvalmat[,j],na.rm=T)+1))
    #abline(lm(xvalmat[,j]~dna),col="black")
    f<-lm(xvalmat[,j]~dna)
    Y <- predict(f, newdata=data.frame(dna))
    coeff<-c(coefficients(lm(xvalmat[,j]~dna)))
    lines(dna, Y,col="blue",lwd=2)
    title(main=colnames(xvalmat)[j],cex.main=2,outer=F,line=0.5)
    print(lm(xvalmat[,j]~dna))
    coeff<-c(coefficients(lm(xvalmat[,j]~dna)))
    r2<-c(paste0("R2=",round(summary(lm(dna~xvalmat[,j]))$r.squared,4)))
    eq=(paste0("y = ", round(coeff[2],4), "*x +", round(coeff[1],4)))
    text(median(dna),max(xvalmat[,j]),2,labels=eq,xpd=T,srt=0,cex=1,col="black",font=1,pos=1)
    text(median(dna),min(xvalmat[,j]),labels=r2,xpd=T,srt=0,cex=1,col="black",font=1,pos=1)
    numberdevice<-c((j%/%6)+1[(j %% 6 != 0)])
    title(main=ifelse(length(nom3)<=6,"Standard",ifelse(length(nom3)>6 & j%%6==1,paste("Standard (",numberdevice,"/",round(length(nom3)/6),")",sep=""),"")),outer=T,line=-2,cex.main=3,font.main=2)
    if (j %% 6 == 0 & j != 1 & j!= ncol(xvalmat)) {dev.new(noRStudioGD = T);par(mfrow=c(2,3))}
  }
} else  {
} 

######################calculation for timeline####################
if (timepoint==1){
  mean_replicat=apply(replicat,2,mean)
  e2=data.frame(mean_replicat,colnames(replicat),b2[order(b2)])
  e2dec=e2[order(e2[,2],decreasing=F),]
  vec2dec=as.vector(e2dec$mean_replicat)
  names(vec2dec)<-e2dec$b2.order.b2.
  vec2dec[order(names(vec2dec))]
  matvec2dec=matrix(vec2dec,ncol=length(unique(names(vec2dec))),byrow=T)
  vecname=as.vector(e2dec$colnames.replicat)
  colnames(matvec2dec)<-unique(names(vec2dec))
  rownames(matvec2dec)<-unique(vecname)
  mate2=matrix(ncol=length(nom3),nrow=length(names(vec2dec))/length(unique(names(vec2dec))))
  rownames(mate2)<-unique(vecname)
  colnames(mate2)<-nom3
  for (p in cond){
    for (q in rownames(matvec2dec)[grep(pattern=p,rownames(matvec2dec))]){
      for (k in nom3){
        mate2[q,k]<-2^(matvec2dec[q,rq]-matvec2dec[q,k])
      }
    }
  }
  mate22=mate2[apply(mate2, 1, function(y) !all(is.na(y))),]
  reslt<-matrix(ncol=length(nb),nrow=length(nom3)*length(cond))
  colnames(reslt)<-paste("D",nb,sep="")
  rownames(reslt)<-apply(expand.grid(primer, cond, rq), 1, paste, collapse = "")
  for (i in cond){
    for (j in primer){
      assign(pos=1,paste(j,i,rq,sep=""),mate2[,j][grep(pattern=i,rownames(mate2))])
      print(paste("---",j,i,rq,"---",sep=""))
      for (n in apply(expand.grid(j, i, rq), 1, paste, collapse = "")){
        mat4=matrix(mate2[,j][grep(pattern=i,rownames(mate2))],byrow=F,ncol=length(unique(name)[grep(pattern=i,unique(name))]))
        mat4[is.na(mat4)]<-0
        print(mat4)
        if (ncol(mat4)<ncol(reslt)){
          nlength<-c((ncol(reslt)-ncol(mat4)))
          print(nlength)
          repnlength<-c(rep(NA,nlength))
          repnlength=matrix(repnlength,nrow=1)
          mat4<-cbind(mat4,repnlength)
        } else {
        }
        rownames(mat4)<- apply(expand.grid(j, i, rq), 1, paste, collapse = "")
        reslt[n,]<-mat4[n,]
      }
    }
  } 
  matfin=t(reslt)
  if (save==1){
    sink("anaqpcr_timepoint.txt",append=T,split=F)
    print(matfin)
    sink()
  } else  {
  } 
  print(reslt)
  if (de==1){
    data.entry(matfin)
  } else  {
  } 
  #nb=(nrow(matfin)-1)
  ####################plot timeline###############"
  for (j in primer){
    dev.new(noRStudioGD = T,width=w, height=h)
    par(mar=c(marge))
    par(mgp=c(space))
    matplot(nb,matfin[,grep(pattern=j,colnames(matfin))],type="o",pch=1,lty=1,main=j,xlab=labelX,ylab=labelY,cex.lab=size_lab,cex.main=size_main,cex.axis=size_axis,ylim=c(0,max(matfin[,grep(pattern=j,colnames(matfin))],na.rm=T)))
    legend(x=xleg,y=max(matfin[,grep(pattern=j,colnames(matfin))],na.rm=T),legend=cond,bty="n",pch=1,col=col,cex=legend_size,xpd=T)
  }
} else  {
##############################Calculation for barplot#####################
  mean_replicat=apply(replicat,2,mean)
  e2=data.frame(mean_replicat,colnames(replicat),b2)
  vec2dec=as.vector(e2$mean_replicat)
  names(vec2dec)<-e2$b2
  vec2dec[order(names(vec2dec))]
  matvec2dec=matrix(vec2dec,ncol=length(unique(names(vec2dec))),byrow=T)
  vecname=as.vector(e2$colnames.replicat)
  colnames(matvec2dec)<-unique(names(vec2dec))
  rownames(matvec2dec)<-unique(vecname)
  mate2=matrix(ncol=length(nom3),nrow=length(names(vec2dec))/length(unique(names(vec2dec))))
  rownames(mate2)<-unique(vecname)
  colnames(mate2)<-nom3
  for (k in nom3){
    mate2[,k]<-2^(matvec2dec[,rq]-matvec2dec[,k])
  }

  mat3=mate2[order(rownames(mate2)),]
  vecmat3<-c(mat3[rownames(mat3)%in%cond])
  mat32=matrix(vecmat3,nrow=length(cond) ,byrow=F)
  colnames(mat32)<- unique(colnames(replicat))[ unique(colnames(replicat))%in%primer]
  cond2=cond[order(cond)]
  rownames(mat32)<-cond2
  mat321=mat32[match(cond,row.names(mat32)),]
  colnames(mat321)<-nom3
  if (save==1){
    sink("anaqpcr.txt",append=T,split=T)
    print(mat321)
    sink()
  } else  {
    print(mat321)} 
  if (de==1){
    
    dfcond<-data.frame(cond)
    mat3212<-data.frame(mat321)
    mat3212<-cbind(dfcond,mat3212)
    data.entry(as.matrix(mat3212))
  } else  {
  } 
  ######################barplot###################
  for (j in primer){
    mat5=matrix(assign(pos=1,paste(j,rq,sep="_"),mat321[,j][rownames(mat321)%in%cond]),byrow=F,ncol=length(j))
    rownames(mat5)<-cond
    colnames(mat5)<-j
    num.decimals <- function(x) {
      stopifnot(class(x)=="numeric")
      x <- sub("0+$","",x)
      x <- sub("^.+[.]","",x)
      nchar(x)
    }
    dev.new(noRStudioGD = T,width=w, height=h)
    if (max(num.decimals(get(paste(j,rq,sep="_"))),na.rm=T)<3){
      par(mgp=c(space_labelY-1,space_axis,space_tick+1))
      par(mar=c(bottom+5,left-1.5,top+1,right-7))
    }else{
      par(mgp=c(space_labelY,space_axis,space_tick+1))
      par(mar=c(bottom+5,left-0.5,top+1,right-7))
    }
    get2<-c(get(paste(j,rq,sep="_")))
    bar<-c(barplot(get(paste(j,rq,sep="_")),main=j,col=colbar,ylab=labelY,cex.lab= size_lab,cex.axis= size_axis,cex.main= size_main,cex.names=cex.names,las=2))
    text(bar,0,labels=ifelse(is.na(get2),"NA",""),pos=3,col="red")
    print(paste(j,rq,sep="_"))
  } 
} 
#END


