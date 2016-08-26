data_file<-na.omit(read.csv("meddata_na.csv"))
data_c<-data_file[!duplicated(data_file),]
data<-data.frame(data_c)[2:10]
data_1<-data[,1]

variance = function(list){
  sum<-0
  for(i in 1:length(list)){
    sum<-sum+sum((list[i]-mean(list))^2)/(length(list)-1)
  }
  return(sum)
}

variance(data_1)