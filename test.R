## Aoman comparison
## 2017.8.25

### step 1 : data preparation
##################################

library(qqman)
#t=read.table("./Chrall.filter.vcf.DP_out.txt.including_zero.ratio",header=F)
#head=read.table("./cmp_material_header",header=F)

cmp=read.table("P01.cmp.all.txt")
load("step1_data_preparation.RData")
colnames(t)[3:25]=as.character(head$V1)


diffr=abs(t$V3-t$V4)


### step 2 : function construction
##################################


  
  s.diffr=sort(diffr)
  r1=rank(diffr)
  
  y=-log10(1- (r1-0.01)/length(r1))
  #plot(1:length(r1), y)
  
  index.y=which(r1> length(r1)*0.99)
  
  threshold0=  -log10(1-  length(r1)*0.99/length(r1))
  print(threshold0)
  walking<-function(start,step=5, threshold= threshold0) {
    score=y[start]
    #print(score)
    left=list()
    left[[1]]=c(start,score)
    i=1
    while(score > threshold) {
      tinterval=y[(start-step*i):(start-step*(i-1)-1)]
      #print(tinterval)
      score <- score+sum(tail(sort(tinterval), 2)) -step
      i=i+1
      #print(step)
      #print(score)
      inp<-c(start-step*(i-1),score)
      left[[i]]=inp
    }
    #print(left)
    
    start_score=y[start]
    right=list()
    right[[1]]=c(start,score)
    i=1
    while(score > threshold) {
      tinterval=y[(start+step*(i-1)+1):(start+step*i)]
      #print(tinterval)
      score <- score+sum(tail(sort(tinterval), 2)) -step
      i=i+1
      #print(step)
      #print(score)
      inp<-c(start-step*(i-1),score)
      right[[i]]=inp
      
    }
    #print(right)
    return(c(left,right))
  }
  
  #start=418
  #out=walking(start)
  index.y=subset(index.y,index.y>5)
  lout=lapply(index.y,walking)
  
  outm=matrix(unlist(lout),ncol =2,byrow=T)
  colnames(outm)=c("ind","value")
  df=data.frame(outm)
  df2=data.frame(  t[df$ind,2],t[df$ind,1], t[df$ind,2],df$value )
  colnames(df2)=c("SNP","CHR","BP","P")
  write.table(df2,file=paste0("Konggan_",colnames(data),"_test.result.txt"),col.names=T,row.names=F,quote=F,sep="\t")
  cutoff=0.05/(dim(df2)[1])
  pdf(file=paste0("Konggan_",colnames(data),".pdf"),wi=20,he=5)
  manhattan(df2,logp=F,genomewideline= -log10(cutoff),suggestiveline=F)
  dev.off()
  
