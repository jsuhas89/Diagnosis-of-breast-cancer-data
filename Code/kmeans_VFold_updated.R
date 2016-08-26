#import actual data
bdata_raw<-na.omit(read.csv("meddata_na.csv"))
bdata1<-bdata_raw[!duplicated(bdata_raw),][2:10]
adata_raw<-na.omit(read.csv("meddata_na.csv"))
adata1<-adata_raw[!duplicated(adata_raw),]

#initialization
sum<-0
no_of_centroids<-2
no_of_iterations<-10
threshold<-0.1
random<-NULL
centroids_initial<-list()
random_again<-NULL
centroids<-list()
centroids_list<-list()
centroids_actual<-list()
centroids_finder<-0
idata<-NULL
jdata<-NULL
centroids_before<-list()
min<-list()
mal<-0
ben<-0
TP<-0
FP<-0
PPV<-0
cc<-0
data1<-adata1[1:68,]
data1_lim<-bdata1[1:68,]
data2<-adata1[69:136,]
data2_lim<-bdata1[69:136,]
data3<-adata1[137:204,]
data3_lim<-bdata1[137:204,]
data4<-adata1[205:272,]
data4_lim<-bdata1[205:272,]
data5<-adata1[273:340,]
data5_lim<-bdata1[273:340,]
data6<-adata1[341:408,]
data6_lim<-bdata1[341:408,]
data7<-adata1[409:476,]
data7_lim<-bdata1[409:476,]
data8<-adata1[477:544,]
data8_lim<-bdata1[477:544,]
data9<-adata1[545:612,]
data9_lim<-bdata1[545:612,]
data10<-adata1[613:683,]
data10_lim<-bdata1[613:675,]
min1<-list()
centroids_finder1<-0

bdata<-rbind(data2_lim,data3_lim,data4_lim,data5_lim,data6_lim,data7_lim,data8_lim,data9_lim,data10_lim)
adata<-rbind(data2,data3,data4,data5,data6,data7,data8,data9,data10)

no_of_attributes<-length(colnames(bdata))

#Euclidean distance function
distance = function(item1,item2){
  return(sqrt(sum((item1 - item2) ^ 2)))
}


#pick initial k random centroids
random<-sample(1:length(rownames(bdata)), no_of_centroids)
for(i in 1:no_of_centroids){
  centroids_initial[[i]]<-bdata[random[i],]
}


#pick k random centroids from our data set again(these are our actual initializations)
random_again<-sample(1:length(rownames(bdata)),no_of_centroids)
for(i in 1:no_of_centroids){
  centroids[[i]]<-bdata[random_again[i],]
}


for(i in 1:no_of_iterations){
  #centroids_list contains k centroids each of which is initialized to 0
  for(m in 1:no_of_centroids){
    centroids_list[[m]]<-matrix(c(0),length(rownames(bdata)),no_of_attributes)
  }
  for(n in 1:no_of_centroids){
    centroids_actual[[n]]<-matrix(c(0),length(rownames(adata)),11)
  }
  
  #initialize the index of the centroids(In R, index starts from 1)
  centroids_index<-rep(c(1),no_of_centroids) 
  
  
  #assign our data points to the relevant centroid
  for(i in 1:length(rownames(bdata))){
    for(j in 1:no_of_centroids){
      min[j]<-distance(bdata[i,],centroids[[j]])
    }
    centroids_finder<-which.min(min)
    
    
    #get the ith data into a vector
    for(k in 1:no_of_attributes){    
      idata[k]<-bdata[i,k] 
    }
    for(l in 1:11){    
      jdata[l]<-adata[i,l] 
    }
    centroids_list[[centroids_finder]][centroids_index[centroids_finder],]<-idata
    centroids_actual[[centroids_finder]][centroids_index[centroids_finder],]<-jdata
    
    #increment the index to store the next dataset if any    
    centroids_index[centroids_finder]<-centroids_index[centroids_finder] + 1 
  }
  
  #decrement the index by 1 since the previous step would have incremented by 1 during last iteration
  for(i in 1:length(centroids_index)){
    centroids_index[i]<-centroids_index[i]-1
  }
  
  
  #calculate the new set of centroids
  centroids_before<-centroids
  for(i in 1:no_of_centroids){ 
    for(j in 1:no_of_attributes){
      sum<-0
      for(k in 1:(centroids_index[i])){
        sum<-sum + centroids_list[[i]][k,j]
      }
      sum<-sum/(centroids_index[i])
      centroids[[i]][j]<-sum
    }
  }
  
  
  #calculate f0 as per the algorithm
  f0<-0
  for(i in 1:length(centroids_before)){ 
    f0_each<-0
    for(j in 1:length(centroids_initial)){
      f0_each<-f0_each + distance(centroids_before[[i]],centroids_initial[[j]])
    }  
    f0<-f0 + f0_each
  }
  
  
  #calculate f1 as per the algorithm  
  f1<-0
  for(i in 1:length(centroids)){
    f1_each<-0
    for(j in 1:length(centroids_before)){
      f1_each<-f1_each + distance(centroids[[i]],centroids_before[[j]])
    }
    f1<-f1 + f1_each
  }
  change_per<-(f1-f0)/f0
  if(abs(change_per)<=threshold){ 
    break
  }
  centroids_initial<-centroids_before 
}

#print centroid indexes
print(centroids_index)


#print the list of centroids
#print(centroids_list)

for(k in 1:no_of_centroids){
  if(k == which.max(centroids_index))
  { 
    cc[k]<-2
    cat(sprintf("Centroid %f is: benign, value is:%f\n",k,cc[k]))
  }
  else
  {
    cc[k]<-4
    cat(sprintf("Centroid %f is: malignant, value is:%f \n",k,cc[k]))
  }
}

#add D1 to the new dataset and calculate the PPV
for(i in 1:length(rownames(data1))){
  for(j in 1: no_of_centroids){
    min1[j]<-distance(data1_lim[i,],centroids[[j]])
  }
  centroids_finder1<-which.min(min1)
  cat(sprintf("dataset %f is associated to centroid %f\n",i,centroids_finder1))
  
  if(centroids_finder1 == 1)
  {
    if((cc[centroids_finder1]==2) & (data1[i,11]==2))
    {
      TP<-TP+1
    }
    else if((cc[centroids_finder1]==4) & (data1[i,11]==4))
    {
      TP=TP+1
    }
    else
    {
      FP<-FP+1
    }
  }
  
  if(centroids_finder1 == 2)
  {
    if((cc[centroids_finder1]==2) & (data1[i,11]==2))
    {
      TP<-TP+1
    }
    else if((cc[centroids_finder1]==4) & (data1[i,11]==4))
    {
      TP=TP+1
    }
    else
    {
      FP<-FP+1
    }
  }
}
cat(sprintf("TP is: %f\n",TP))
cat(sprintf("FP is: %f\n",FP))
PPV = TP/(TP+FP)
cat(sprintf("PPV of data1 is: %f\n",PPV))