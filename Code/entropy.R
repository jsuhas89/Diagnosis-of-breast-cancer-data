library("entropy")
data_file<-na.omit(read.csv("meddata_na.csv"))
data_c<-data_file[!duplicated(data_file),]
data<-data.frame(data_c)[2:10]

#calculate counts of values in each attribute
counts<-list()
for(i in 1:length(colnames(data))){
  c<-NULL
  #initialize count to 0
  for(x in 1:10){c[x]<-0}
  for(j in 1:10){
    for(k in 1:length(rownames(data))){
      if(data[k,i]==j){
        c[j]<-c[j]+1
      }
    }
  }
  counts[[i]]<-c
}

#initialize the entropies
entropies<-NULL

#calculating entropies for all columns
for(i in 1:length(colnames(data))){
  entropies[i]<-entropy.empirical(counts[[i]],unit=c("log2"))
  print(entropies[i])
}

