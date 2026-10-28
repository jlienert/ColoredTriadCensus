source("colored.triad.census.R")

#Generate a network
set.seed(9414)
bmat<-rgnm(1,10,30)

#One color colored triad census
col.1<-colored.triad.census(bmat,col =rep(1,10),directed = TRUE)

#Compare with sna
triad.census(bmat)
col.1

#Two and three colored triad censuses
col.2<-colored.triad.census(bmat,col =rbinom(10,1,.25),directed = TRUE)
col.3<-colored.triad.census(bmat,col =rbinom(10,2,.25),directed = TRUE)

#Heatmap of iso classes against coloring
iso.class<-sapply(strsplit(names(col.2),"-"),function(x) x[[1]])
col.class<-sapply(strsplit(names(col.2),"-"),function(x) x[[2]])
outmat<-matrix(0,nrow=length(unique(iso.class)),ncol=length(unique(col.class)),dimnames = list(unique(iso.class),unique(col.class)))

for(i in 1:length(col.2)){
  outmat[iso.class[i],col.class[i]]<-col.2[i]
}

heatmap(outmat,Rowv=NA,Colv=NA)

