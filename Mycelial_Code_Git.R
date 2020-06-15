
#Final pydiflumetofen baseline mycelial growth code, for determining EC50s

setwd(dir="C:/Users/mbreu/OneDrive/Final_Pydiflumetofen_030520")

ipak <- function( pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[,"Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("broom","stringr", "tidyverse","ggplot2", "drc","lattice","car", "lme4", "lsmeans", "plyr", "plotrix", "knitr", "ggplot2", "lmtest", "lmerTest", "Rmisc", "gridExtra", "plotly", "webshot", "ggpmisc", "ggsci","scales")
ipak(packages)
options(scipen = 999) 

ls() #lists all
rm(list = ls(all=TRUE))


AdepidynMycelialRG<-read_delim(file = "Combined_073119.csv",delim=",",na=".")
AdepidynInfo<-read_delim(file = "IsolateMetaData_073119.csv",delim=",",na=".")
AdepidynAll<-left_join(AdepidynMycelialRG,AdepidynInfo, by="Isolate")
View(AdepidynAll)


View(AdepidynMycelialOthers)
##process data to go from raw diameters to a relative growth for each isolate
Adepidynspread<-AdepidynAll %>%
  spread(key="measurement_rep", value=length)


AdepidynFull<-Adepidynspread %>% 
  mutate(average_length=(((Adepidynspread$one+Adepidynspread$two)/2)-5.35))  %>%
  group_by(Isolate, set) %>% 
  mutate(avgC=mean(average_length[ppm==0])) %>% 
  mutate(relative=average_length/avgC) %>% 
  filter(Species==" F.graminearum")

AdepidynFull$Isolate <- factor(AdepidynFull$Isolate) #sets as factor

View(AdepidynFull)

write.csv(x=AdepidynFull, file="Mycelial_processed.csv")

AdepidynFull
#split into each individual set,so that EC50s are determined separately for each run individually
 set1<-AdepidynFull %>% 
   filter(set==1)
 set2<-AdepidynFull %>%   
   filter(set==2)
 set3<-AdepidynFull %>% 
   filter(set==3)
 set4<-AdepidynFull %>% 
   filter(set==4)
     
 ###EC50 loop SET1 
 
 isolates1<-unique(set1$Isolate)
 LL3set1abs<-NULL
 
 for (i in seq_along(isolates1))
 {
   curve.LL3<- drm(100 * relative[set1$Isolate==isolates1[[i]]] ~ ppm[set1$Isolate==isolates1[[i]]], data = set1, fct = LL.3(), na.action = na.omit)
   ED50LL3_i<-data.frame(ED(curve.LL3,respLev=c(50), type=c("absolute"), interval="delta", level=0.95))
   names_i<-cbind.data.frame(ED50LL3_i,isolates1[[i]])
   LL3set1abs<-rbind.data.frame(LL3set1abs,names_i)
 }
 
 ###set2
 
 isolates2<-unique(set2$Isolate)
 
 LL3set2abs<-NULL
 
 
 for (i in seq_along(isolates2))
 {
   curve.LL3<- drm(100 * relative[set2$Isolate==isolates2[[i]]] ~ ppm[set2$Isolate==isolates2[[i]]], data = set2, fct = LL.3(), na.action = na.omit)
   ED50LL3_i<-data.frame(ED(curve.LL3,respLev=c(50), type=c("absolute"), interval="delta", level=0.95))
   names_i<-cbind.data.frame(ED50LL3_i,isolates2[[i]])
   LL3set2abs<-rbind.data.frame(LL3set2abs,names_i)
 }
 
 
 
 ##set3
 
 isolates3<-unique(set3$Isolate)
 
 LL3set3abs<-NULL
 
 for (i in seq_along(isolates3))
 {
   curve.LL3<- drm(100 * relative[set3$Isolate==isolates3[[i]]] ~ ppm[set3$Isolate==isolates3[[i]]], data = set3, fct = LL.3(), na.action = na.omit)
   ED50LL3_i<-data.frame(ED(curve.LL3,respLev=c(50), type=c("absolute"), interval="delta", level=0.95))
   names_i<-cbind.data.frame(ED50LL3_i,isolates3[[i]])
   LL3set3abs<-rbind.data.frame(LL3set3abs,names_i)
 }
 
 
 
 ##set4
 
 isolates4<-unique(set4$Isolate)
 
 LL3set4abs<-NULL
 
 for (i in seq_along(isolates4))
 {
   curve.LL3<- drm(100 * relative[set4$Isolate==isolates4[[i]]] ~ ppm[set4$Isolate==isolates4[[i]]], data = set4, fct = LL.3(), na.action = na.omit)
   ED50LL3_i<-data.frame(ED(curve.LL3,respLev=c(50), type=c("absolute"), interval="delta", level=0.95))
   names_i<-cbind.data.frame(ED50LL3_i,isolates4[[i]])
   LL3set4abs<-rbind.data.frame(LL3set4abs,names_i)
 }
 
 ##renaming 
 colnames(LL3set1abs)[5]<-"isolates[[i]]"
 colnames(LL3set2abs)[5]<-"isolates[[i]]"
 colnames(LL3set3abs)[5]<-"isolates[[i]]"
 colnames(LL3set4abs)[5]<-"isolates[[i]]"

 LL3set1abs$set=1
 LL3set2abs$set=2
 LL3set3abs$set=3
 LL3set4abs$set=4

 
 
 #combined 5 sets into one master data frame     

  Combined_Myc_EC50<-rbind.data.frame(LL3set1abs,LL3set2abs,LL3set3abs,LL3set4abs)
 colnames(Combined_Myc_EC50)<-c("Estimate","StdError","LowerBound","UpperBound","Isolate","Set")
 
write.csv(x=Combined_Myc_EC50, file="AllEc50MYC.csv")
 
 
View(Combined_Myc_EC50)
 #get the average for each isolate
 SummaryEC50<-Combined_Myc_EC50 %>% 
   group_by(Isolate) %>%  
   summarise(Mean=mean(Estimate),Median=median(Estimate),StdDev=sd(Estimate),n=n(),se=StdDev/sqrt(n),qt=qt(0.975,(n-1)),CI.min=Mean-(se*qt),CI.max=Mean+(se*qt))
 
 
View(SummaryEC50)
 #95% CI, alpha=0.05, qt calc number from distribution using (0.975,1)=12.7062
 
 write.csv(x=SummaryEC50, file="SummaryEC50Myc_LL3.csv")
 

 
     
     