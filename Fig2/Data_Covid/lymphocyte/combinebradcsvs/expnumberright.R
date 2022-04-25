datadfall<-read.csv("allumsz.csv") #this is created by comcsv
datadfall$Sitename=datadfall$Experiment
datadfall1=datadfall
Sitenames<-read.csv("site_names.csv",header = FALSE) #this is created by comcsv
datadfall1$ss=factor(datadfall1$Sitename,Sitenames$V1,levels = Sitenames$V2)
my_data1 <- list()
numexp=list()
nams=levels(unlist(datadfall$Compound))
for (i in nams){
#i="astemizole"
  my_data= subset(datadfall, Compound == i)
  my_dataif=my_data

# m=my_data$Experiment[1]
# my_data$Experiment[1]=1
  for (j in 2:length(my_data$Experiment)){
    my_data$Sitename[j]=my_data$Experiment[j]
      #print(m)
    my_data$Experiment[1]=1
  if (my_dataif$Experiment[j]==my_dataif$Experiment[j-1]){
      my_data$Experiment[j]=my_data$Experiment[j-1]
      
      print("hiiiiiiiiii")}
  else if (my_dataif$Experiment[j]!= my_dataif$Experiment[j-1]){
        my_data$Experiment[j]=my_data$Experiment[j-1]+1
        print("nooooooooo")}
    else
      print("kkkkkkkkkkkk")
  }
# if (my_data$Experiment[j]!=m){
#      my_data$Experiment[j]=my_data$Experiment[j-1]+1
#  }    
  my_data1=rbind(my_data1,my_data)

  
  #my_data1$sitename2=factor(my_data1$Sitename,Sitenames$V1,levels = Sitenames$V3)
  
  #d1=levels(my_data$Experiment)
  #expe=unique(unlist(my_data$Experiment))
  #numexp=length(expe)
  #print(i)
  #print(numexp)
  #print(expe)
  
}
# for (i in 1:length(my_data$Experiment)){
#   if (my_data$Sitename == Sitenames$V2){
# my_data$ff=Sitenames$V1
# }}
my_data2=my_data1[-6]
    write.csv(my_data1,"allumscorrectnameandexperimentnum1.csv")
my_data1$sitenames=factor(my_data1$Sitename,Sitenames$V1,levels = Sitenames$V2)
my_data1$sitenumber=factor(my_data1$Sitename,Sitenames$V3,levels = Sitenames$V2)
write.csv(my_data1,"all_with_site_name.csv")

  