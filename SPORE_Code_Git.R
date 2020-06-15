

##Final code for published pydiflumetofen baseline sensitivity paper- for spore germ EC50 determination 

ipak <- function( pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[,"Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("broom","dr4pl","stringr", "tidyverse","ggplot2", "drc","lattice","car", "lme4", "lsmeans", "plyr", "plotrix", "knitr", "ggplot2", "lmtest", "lmerTest", "Rmisc", "gridExtra", "plotly", "webshot", "ggpmisc", "ggsci","scales")
ipak(packages)

rm(list = ls(all=TRUE)) # removes all variables in the global environment 
setwd(dir="C:/Users/mbreu/Desktop/pydi_git/")

CombinedRaw<-read_delim(file = "CombinedSporeFinal_GIT.csv",delim=",",na=".")

Metadata<-read_delim(file = "IsolateMetaData_073119.csv",delim=",",na=".")

AdepidynSpore<-left_join(CombinedRaw,Metadata, by="Isolate")
View(Data)

sporepercent<-AdepidynSpore %>% 
  mutate(percentgerm = AdepidynSpore$YES/(AdepidynSpore$YES+ AdepidynSpore$NO)) %>% 
  group_by(Isolate, set) %>% 
  mutate(avgcontrol=mean(percentgerm[ppm==0])) %>% 
  mutate(relative=percentgerm/avgcontrol) 





write.csv(x=sporepercent, file="Spore_processed.csv")

View(sporepercent)

 #split into each individual set,so that EC50s are determined separately for each run individually
     set1<-sporepercent %>% 
       filter(set==1)
     set2<-sporepercent%>%   
       filter(set==2)
     set3<-sporepercent %>% 
       filter(set==3)
     set4<-sporepercent %>% 
       filter(set==4)
     set5<-sporepercent %>% 
       filter(set==5)

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

##set5
     
isolates5<-unique(set5$Isolate)
     
LL3set5abs<-NULL
     
for (i in seq_along(isolates5))
     {
       curve.LL3<- drm(100 * relative[set5$Isolate==isolates5[[i]]] ~ ppm[set5$Isolate==isolates5[[i]]], data = set5, fct = LL.3(), na.action = na.omit)
       ED50LL3_i<-data.frame(ED(curve.LL3,respLev=c(50), type=c("absolute"), interval="delta", level=0.95))
       names_i<-cbind.data.frame(ED50LL3_i,isolates5[[i]])
       LL3set5abs<-rbind.data.frame(LL3set5abs,names_i)
     }
     

#rename columns in each data frame and add set as a variable     
     
     
colnames(LL3set1abs)[5]<-"isolates[[i]]"
colnames(LL3set2abs)[5]<-"isolates[[i]]"
colnames(LL3set3abs)[5]<-"isolates[[i]]"
colnames(LL3set4abs)[5]<-"isolates[[i]]"
colnames(LL3set5abs)[5]<-"isolates[[i]]"
     
LL3set1abs$set=1
LL3set2abs$set=2
LL3set3abs$set=3
LL3set4abs$set=4
LL3set5abs$set=5
     
   
#combined 5 sets into one master data frame     
Combined_SPORE_EC50<-rbind.data.frame(LL3set1abs,LL3set2abs,LL3set3abs,LL3set4abs,LL3set5abs)
colnames(Combined_SPORE_EC50)<-c("Estimate","StdError","LowerBound","UpperBound","Isolate","Set")
View(Combined_SPORE_EC50)
write.csv(x=Combined_SPORE_EC50, file="AllEC50SPORE.csv")
#get the average EC50 for each isolate and summary statistics of two EC50 estimations

SummaryEC50SPORE<-Combined_SPORE_EC50 %>% 
   group_by(Isolate) %>%  
   summarise(Mean=mean(Estimate),Median=median(Estimate),StdDev=sd(Estimate),n=n(),se=StdDev/sqrt(n),qt=qt(0.975,(n-1)),CI.min=Mean-(se*qt),CI.max=Mean+(se*qt)) 
View(SummaryEC50SPORE)
write.csv(x=SummaryEC50SPORE, file="SummaryEC50SPORE_LL3.csv")

     