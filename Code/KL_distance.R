library("entropy")
data_file<-na.omit(read.csv("meddata_na.csv"))
data<-data_file[!duplicated(data_file),][2:10]
#data<-data.frame(data_c)[2:10]

#calculate counts of values in each attribute
counts<-list()
for(i in 1:length(colnames(data))){
  c<-NULL
  #initialize count to 0
  for(x in 1:10){c[x]<-0.0000000000000000000000000000000000000000000000001}
  for(j in 1:10){
    for(k in 1:length(rownames(data))){
      if(data[k,i]==j){
        c[j]<-c[j]+1
      }
    }
  }
  counts[[i]]<-c
}

#KL-distances for attribute pairs
KL<-list()
for(i in 1:length(colnames(data))){
  for(j in 1:length(colnames(data))){
    KL_distance<-KL.plugin(counts[[i]],counts[[j]],unit=c("log2"))
    KL[[paste(c(i,j),collapse="")]]<-KL_distance
  }
}
print(KL)
