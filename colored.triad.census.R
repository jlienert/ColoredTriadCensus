require(network)
library(sna)
require(ergm.count)
require(gtools)
require(Matrix)
require(gplots)
load("uniTri.RData")

kfun<-function(k,output="fourClass",directed=T){
  if(k>=3){
    ef<-choose(k,3)+2*choose(k,2)+k
    c030<-2*choose(k,3)+2*choose(k,2)+k
    ex2<-3*choose(k,3)+4*choose(k,2)+k
    uni<-6*choose(k,3)+6*choose(k,2)+k
  }
  if(k==2){
    ef<-2*choose(k,2)+k
    c030<-2*choose(k,2)+k
    ex2<-4*choose(k,2)+k
    uni<-6*choose(k,2)+k
  }
  if(k==1){
    ef<-k
    c030<-k
    ex2<-k
    uni<-k
  }
  if(output=="fourClass") return(c(ef,c030,ex2,uni))
  if(output=="total") {
    if(directed==T) return(sum(2*ef,c030,6*ex2,7*uni))
    if(directed==F) return(sum(2*ef,2*ex2))
  }
}

tr=function(mat){return(sum(diag(mat)))}

txtmats=function(){
  perms=permutations(3,3,1:3)
  ret=list()
  col=rep(seq(1:6),6,each=1)
  row=rep(seq(1:6),1,each=6)
  for(i in 1:36) ret[[i]]=list(perms[row[i],],perms[col[i],])
  return(ret)
}

samemat=function(mat1,mat2){
  txtmat=txtmats()
  testmats=list()
  check=F
  for(i in 1:length(txtmat)) testmats[[i]]=mat2[txtmat[[i]][[1]],txtmat[[i]][[2]]]
  for(i in 1:length(testmats)) if(all(mat1==testmats[[i]])) check=T
  return(check)
}

makeUniqueTri=function(triangle,colorCombs){
  out=1:nrow(colorCombs)
  uniqueTriangles=list()
  for(i in 1:nrow(colorCombs)){uniqueTriangles[[i]]= triangle ; diag(uniqueTriangles[[i]])=colorCombs[i,]+1}
  if(length(uniqueTriangles)>1){
    for(i in (length(uniqueTriangles)-1):1){
      for (j in (i+1):length(uniqueTriangles)){
        if(samemat(uniqueTriangles[[i]],uniqueTriangles[[j]])) out[j]=i
      }
    }
  }
  return(out)
}

remEx=function(triad,colorCombs,testM,testE0,testC,cmat){
  T003Col=c()
  for(i in 1:nrow(colorCombs)){
    if(triad=="003"){
      T003Col[i]=tr((t(cmat[[colorCombs[i,2]]])*cmat[[colorCombs[i,1]]]*testE0)%*%(t(cmat[[colorCombs[i,3]]])*cmat[[colorCombs[i,2]]]*testE0)%*%(t(cmat[[colorCombs[i,1]]])*cmat[[colorCombs[i,3]]]*testE0))
      names(T003Col)[i]=paste("T003",paste(colorCombs[i,],collapse="",sep=""),sep="-")
      triangle=matrix(c(0,0,0,0,0,0,0,0,0),nrow=3)
    }
    
    if(triad=="012"){
      T003Col[i]=tr((t(cmat[[colorCombs[i,2]]])*cmat[[colorCombs[i,1]]]*testE0)%*%(t(cmat[[colorCombs[i,3]]])*cmat[[colorCombs[i,2]]]*testE0)%*%(t(cmat[[colorCombs[i,1]]])*cmat[[colorCombs[i,3]]]*testC))
      names(T003Col)[i]=paste("T012",paste(colorCombs[i,],collapse="",sep=""),sep="-")
      triangle=matrix(c(0,0,0,0,0,0,1,0,0),nrow=3,byrow=T)
    }
    
    if(triad=="102"){
      T003Col[i]=tr((t(cmat[[colorCombs[i,2]]])*cmat[[colorCombs[i,1]]]*testE0)%*%(t(cmat[[colorCombs[i,3]]])*cmat[[colorCombs[i,2]]]*testM)%*%(t(cmat[[colorCombs[i,1]]])*cmat[[colorCombs[i,3]]]*testE0))
      names(T003Col)[i]=paste("T102",paste(colorCombs[i,],collapse="",sep=""),sep="-")
      triangle=matrix(c(0,0,0,0,0,1,0,1,0),nrow=3)
    }
    
    if(triad=="021D"){
      T003Col[i]=tr((t(cmat[[colorCombs[i,2]]])*cmat[[colorCombs[i,1]]]*testC)%*%(t(cmat[[colorCombs[i,3]]])*cmat[[colorCombs[i,2]]]*testE0)%*%(t(cmat[[colorCombs[i,1]]])*cmat[[colorCombs[i,3]]]*t(testC)))
      names(T003Col)[i]=paste("T021D",paste(colorCombs[i,],collapse="",sep=""),sep="-")
      triangle=matrix(c(0,1,1,0,0,0,0,0,0),nrow=3,byrow=T)
    }
    
    if(triad=="021U"){
      T003Col[i]=tr((t(cmat[[colorCombs[i,2]]])*cmat[[colorCombs[i,1]]]*t(testC))%*%(t(cmat[[colorCombs[i,3]]])*cmat[[colorCombs[i,2]]]*testE0)%*%(t(cmat[[colorCombs[i,1]]])*cmat[[colorCombs[i,3]]]*testC))
      names(T003Col)[i]=paste("T021U",paste(colorCombs[i,],collapse="",sep=""),sep="-")
      triangle=matrix(c(0,0,0,1,0,0,1,0,0),nrow=3,byrow=T)
    }
    
    if(triad=="021C"){
      T003Col[i]=tr((t(cmat[[colorCombs[i,2]]])*cmat[[colorCombs[i,1]]]*testC)%*%(t(cmat[[colorCombs[i,3]]])*cmat[[colorCombs[i,2]]]*testE0)%*%(t(cmat[[colorCombs[i,1]]])*cmat[[colorCombs[i,3]]]*testC))
      names(T003Col)[i]=paste("T021C",paste(colorCombs[i,],collapse="",sep=""),sep="-")
      triangle=matrix(c(0,1,0,0,0,0,1,0,0),nrow=3,byrow=T)
    }
    
    if(triad=="111D"){
      T003Col[i]=tr((t(cmat[[colorCombs[i,2]]])*cmat[[colorCombs[i,1]]]*testC)%*%(t(cmat[[colorCombs[i,3]]])*cmat[[colorCombs[i,2]]]*testM)%*%(t(cmat[[colorCombs[i,1]]])*cmat[[colorCombs[i,3]]]*testE0))
      names(T003Col)[i]=paste("T111D",paste(colorCombs[i,],collapse="",sep=""),sep="-")
      triangle=matrix(c(0,1,0,0,0,1,0,1,0),nrow=3,byrow=T)
    }
    
    if(triad=="111U"){
      T003Col[i]=tr((t(cmat[[colorCombs[i,2]]])*cmat[[colorCombs[i,1]]]*t(testC))%*%(t(cmat[[colorCombs[i,3]]])*cmat[[colorCombs[i,2]]]*testM)%*%(t(cmat[[colorCombs[i,1]]])*cmat[[colorCombs[i,3]]]*testE0))
      names(T003Col)[i]=paste("T111U",paste(colorCombs[i,],collapse="",sep=""),sep="-")
      triangle=matrix(c(0,0,0,1,0,1,0,1,0),nrow=3,byrow=T)
    }
    
    if(triad=="030T"){
      T003Col[i]=tr((t(cmat[[colorCombs[i,2]]])*cmat[[colorCombs[i,1]]]*t(testC))%*%(t(cmat[[colorCombs[i,3]]])*cmat[[colorCombs[i,2]]]*t(testC))%*%(t(cmat[[colorCombs[i,1]]])*cmat[[colorCombs[i,3]]]*testC))
      names(T003Col)[i]=paste("T030T",paste(colorCombs[i,],collapse="",sep=""),sep="-")
      triangle=matrix(c(0,0,0,1,0,0,1,1,0),nrow=3,byrow=T)
    }
    
    if(triad=="030C"){
      T003Col[i]=tr((t(cmat[[colorCombs[i,2]]])*cmat[[colorCombs[i,1]]]*testC)%*%(t(cmat[[colorCombs[i,3]]])*cmat[[colorCombs[i,2]]]*testC)%*%(t(cmat[[colorCombs[i,1]]])*cmat[[colorCombs[i,3]]]*testC))
      names(T003Col)[i]=paste("T030C",paste(colorCombs[i,],collapse="",sep=""),sep="-")
      triangle=matrix(c(0,1,0,0,0,1,1,0,0),nrow=3,byrow=T)
    }
    
    if(triad=="201"){
      T003Col[i]=tr((t(cmat[[colorCombs[i,2]]])*cmat[[colorCombs[i,1]]]*testM)%*%(t(cmat[[colorCombs[i,3]]])*cmat[[colorCombs[i,2]]]*testE0)%*%(t(cmat[[colorCombs[i,1]]])*cmat[[colorCombs[i,3]]]*testM))
      names(T003Col)[i]=paste("T201",paste(colorCombs[i,],collapse="",sep=""),sep="-")
      triangle=matrix(c(0,1,1,1,1,0,1,0,0),nrow=3)
    }
    
    if(triad=="120D"){
      T003Col[i]=tr((t(cmat[[colorCombs[i,2]]])*cmat[[colorCombs[i,1]]]*testC)%*%(t(cmat[[colorCombs[i,3]]])*cmat[[colorCombs[i,2]]]*testM)%*%(t(cmat[[colorCombs[i,1]]])*cmat[[colorCombs[i,3]]]*t(testC)))
      names(T003Col)[i]=paste("T120D",paste(colorCombs[i,],collapse="",sep=""),sep="-")
      triangle=matrix(c(0,1,1,0,0,1,0,1,0),nrow=3,byrow=T)
    }
    
    if(triad=="120U"){
      T003Col[i]=tr((t(cmat[[colorCombs[i,2]]])*cmat[[colorCombs[i,1]]]*t(testC))%*%(t(cmat[[colorCombs[i,3]]])*cmat[[colorCombs[i,2]]]*testM)%*%(t(cmat[[colorCombs[i,1]]])*cmat[[colorCombs[i,3]]]*testC))
      names(T003Col)[i]=paste("T120U",paste(colorCombs[i,],collapse="",sep=""),sep="-")
      triangle=matrix(c(0,0,0,1,0,1,1,1,0),nrow=3,byrow=T)
    }
    
    if(triad=="120C"){
      T003Col[i]=tr((t(cmat[[colorCombs[i,2]]])*cmat[[colorCombs[i,1]]]*testC)%*%(t(cmat[[colorCombs[i,3]]])*cmat[[colorCombs[i,2]]]*testM)%*%(t(cmat[[colorCombs[i,1]]])*cmat[[colorCombs[i,3]]]*testC))
      names(T003Col)[i]=paste("T120C",paste(colorCombs[i,],collapse="",sep=""),sep="-")
      triangle=matrix(c(0,1,0,0,0,1,1,1,0),nrow=3,byrow=T)
    }
    
    if(triad=="210"){
      T003Col[i]=tr((t(cmat[[colorCombs[i,2]]])*cmat[[colorCombs[i,1]]]*testM)%*%(t(cmat[[colorCombs[i,3]]])*cmat[[colorCombs[i,2]]]*testM)%*%(t(cmat[[colorCombs[i,1]]])*cmat[[colorCombs[i,3]]]*testC))
      names(T003Col)[i]=paste("T210",paste(colorCombs[i,],collapse="",sep=""),sep="-")
      triangle=matrix(c(0,1,0,1,0,1,1,1,0),nrow=3,byrow=T)
    }
    
    if(triad=="300"){
      T003Col[i]=tr((t(cmat[[colorCombs[i,2]]])*cmat[[colorCombs[i,1]]]*testM)%*%(t(cmat[[colorCombs[i,3]]])*cmat[[colorCombs[i,2]]]*testM)%*%(t(cmat[[colorCombs[i,1]]])*cmat[[colorCombs[i,3]]]*testM))
      names(T003Col)[i]=paste("T300",paste(colorCombs[i,],collapse="",sep=""),sep="-")
      triangle=matrix(c(0,1,1,1,0,1,1,1,0),nrow=3)
    }
  }
  imp=uniTri[names(uniTri)==paste((max(colorCombs)),triad,sep=" ")][[1]]
  #imp=makeUniqueTri(triangle=triangle,colorCombs)
  T003Col2=c()
  for(i in unique(imp)) 
  {
    T003Col2[length(T003Col2)+1]=T003Col[imp==i][1]
    names(T003Col2)[length(T003Col2)]=names(T003Col)[imp==i][1]
    
    if(triad %in% c("003","300")){
      if (colorCombs[i,1]==colorCombs[i,2] | colorCombs[i,1]==colorCombs[i,3] | colorCombs[i,2]==colorCombs[i,3]) T003Col2[length(T003Col2)]=T003Col2[length(T003Col2)]/2
      if (colorCombs[i,1]==colorCombs[i,2] & colorCombs[i,1]==colorCombs[i,3] & colorCombs[i,2]==colorCombs[i,3]) T003Col2[length(T003Col2)]=T003Col2[length(T003Col2)]/3  
    }
    
    if(triad=="030C") if (colorCombs[i,1]==colorCombs[i,2] & colorCombs[i,1]==colorCombs[i,3] & colorCombs[i,2]==colorCombs[i,3]) T003Col2[length(T003Col2)]=T003Col2[length(T003Col2)]/3 
    
    if(triad %in% c("102","021D","021U","201","120D","120U")) if (colorCombs[i,3]==colorCombs[i,2]) T003Col2[length(T003Col2)]=T003Col2[length(T003Col2)]/2
    
  }
  print(triad)
  return(T003Col2)
}

colored.triad.census=function(mat,col,directed=F){
  #Adjacency matrix
  testA=mat[,]
  #Asymmetric edges only
  testC=testA-t(testA)
  testC=matrix(pmax(0,testC),nrow=nrow(mat),byrow=F)
  #Mutual edges only
  testM=testA-testC
  #Symmetrized full matrix, and its complement
  testE=testA+t(testA)
  testE=pmin(testE,1)
  testE0=1-testE
  diag(testE0)=0
  
  testA=Matrix(testA)
  testC=Matrix(testC)
  testM=Matrix(testM)
  testE=Matrix(testE)
  testE0=Matrix(testE0)
  
  #Colors
  #number of colors
  colnum=length(unique(col))
  colors=1:length(unique(col))
  names(colors)=unique(col)
  colNums=colors[match(col,names(colors))]
  #Makes matrix for each color where for every node that is that color, that entire column is 1
  cmat=list()
  for(i in 1:colnum){
    cmat[[i]]=Matrix(rep(colNums==i,length(col)),nrow=length(col),byrow=FALSE)*1
    cmat[[i]][,]
  }
  colorCombs=permutations(colnum,3,unique(colors),repeats.allowed = TRUE)
  
  #003
  T003Col=remEx(triad="003",colorCombs=colorCombs,testM=testM,testE0=testE0,testC=testC,cmat=cmat)
  #102
  T102Col=remEx("102",colorCombs=colorCombs,testM=testM,testE0=testE0,testC=testC,cmat=cmat)
  #201
  T201Col=remEx("201",colorCombs=colorCombs,testM=testM,testE0=testE0,testC=testC,cmat=cmat)
  #300
  T300Col=remEx("300",colorCombs=colorCombs,testM=testM,testE0=testE0,testC=testC,cmat=cmat)
  x=c(T003Col,T102Col,T201Col,T300Col)
  
  if(directed==T) {
    #012
    T012Col=remEx("012",colorCombs=colorCombs,testM=testM,testE0=testE0,testC=testC,cmat=cmat)
    #021D
    T021DCol=remEx("021D",colorCombs=colorCombs,testM=testM,testE0=testE0,testC=testC,cmat=cmat)
    #021C
    T021CCol=remEx("021C",colorCombs=colorCombs,testM=testM,testE0=testE0,testC=testC,cmat=cmat)
    #021U
    T021UCol=remEx("021U",colorCombs=colorCombs,testM=testM,testE0=testE0,testC=testC,cmat=cmat)
    #030T
    T030TCol=remEx("030T",colorCombs=colorCombs,testM=testM,testE0=testE0,testC=testC,cmat=cmat)
    #030C
    T030CCol=remEx("030C",colorCombs=colorCombs,testM=testM,testE0=testE0,testC=testC,cmat=cmat)
    #120D
    T120DCol=remEx("120D",colorCombs=colorCombs,testM=testM,testE0=testE0,testC=testC,cmat=cmat)
    #120U
    T120UCol=remEx("120U",colorCombs=colorCombs,testM=testM,testE0=testE0,testC=testC,cmat=cmat)
    #120C
    T120CCol=remEx("120C",colorCombs=colorCombs,testM=testM,testE0=testE0,testC=testC,cmat=cmat)
    #210
    T210Col=remEx("210",colorCombs=colorCombs,testM=testM,testE0=testE0,testC=testC,cmat=cmat)
    #111D
    T111DCol=remEx("111D",colorCombs=colorCombs,testM=testM,testE0=testE0,testC=testC,cmat=cmat)
    #111U
    T111UCol=remEx("111U",colorCombs=colorCombs,testM=testM,testE0=testE0,testC=testC,cmat=cmat)
    x=c(T003Col,T102Col,T012Col,T201Col,T300Col,T021DCol,T021UCol,T021CCol,T030TCol,T030CCol,T120DCol,T120UCol,T120CCol,T210Col,T111DCol,T111UCol)
  }
  return(x)
}


