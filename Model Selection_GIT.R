##Model Selection Script for Baseline Sensitivity of Pydiflumetofen for MYCELIAL GROWTH data##

rm(list = ls(all=TRUE)) 
options(scipen = 999) 

##function to check, then load any packages that we desire
ipak <- function( pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[,"Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("broom","stringr", "tidyverse","ggplot2", "drc","lattice","car", "lme4", "lsmeans", "plyr", "plotrix", "knitr", "ggplot2", "lmtest", "lmerTest", "Rmisc", "gridExtra", "plotly", "webshot", "ggpmisc", "ggsci","scales")
ipak(packages)



##Load data

AdepidynMycelialRG<-read_delim(file = "Combined_073119_GIT.csv",delim=",",na=".")
#AdepidynInfo<-read_delim(file = "IsolateMetaData_073119.csv",delim=",",na=".")
ID<- read_delim(file = "JoinedIsolateList.csv",delim=",",na=".")
IsolateData<-unique(AdepidynMycelialRG$Isolate)
AdepidynAll<-left_join(AdepidynMycelialRG,ID, by="Isolate")
 
View(AdepidynFinal)
#what species do I have data on? 
Species<-unique(AdepidynAll$species)
View(used)


used<-unique(AdepidynFinal$Isolate)

##remove isolates that are different species that I won't use 
AdepidynFinal<-AdepidynMycelialRG %>% 
  filter(Isolate!="WMISO_2L-3-22") %>% 
  filter(Isolate!="NTR1_102") %>% 
  filter(Isolate!="MI-SF-1-17") %>% 
  filter(Isolate!="61A") %>% 
  filter(Isolate!="C1E3C") %>% 
  filter(Isolate!="C2A") %>% 
  filter(Isolate!="C9E4C") %>% 
  filter(Isolate!="P3-R5-1")  


View(AdepidynMycelialOthers)
##process data to go from raw diameters to a relative growth for each isolate
Adepidynspread<- AdepidynFinal %>%
  spread(key="measurement_rep", value=length)
AdepidynFull <- Adepidynspread %>% 
  mutate(average_length=(((Adepidynspread$one+Adepidynspread$two)/2)-5.35))  %>%
  group_by(Isolate, set) %>%
  mutate(avgC=mean(average_length[ppm==0])) %>%
  mutate(relative=average_length/avgC) %>%
  filter(set==3)
#filtered for a single set (run) at a time

AdepidynFull$Isolate <- factor(AdepidynFull$Isolate) #sets as factor

##creates function that will determine the AIC and Loglikhood for each isolate individually(group_by isolate), input is the model
model.stats <- function(func.number){
  nested_df <- AdepidynFull %>% 
    group_by(Isolate) %>% 
    nest() %>% 
    mutate(drc = map(data, ~ drm(100*relative ~ ppm, fct = func.list[[func.number]], data = .))) 
  stats <- map(nested_df$drc, ~data.frame(loglik = logLik(.x), 
                                          aic = AIC(.x))) %>%
    do.call(rbind.data.frame, .)
  stats$model <- func.list[[func.number]]$name
  stats$isolate <- levels(AdepidynFull$Isolate)
  return(stats)
}

###### func.list , would add any additional models here, there position becomes their number in the next step 
func.list <- list(LL.3(), LL.4())

##runs function according to that model number and outputs it in df 
LL3 <- model.stats(1)
LL4 <- model.stats(2)


all.together <- rbind.data.frame(LL3,LL4)


min.aic <- all.together %>%
  group_by(isolate) %>% 
  nest() %>% 
  mutate(best.model = map(data, ~.$model[which.min(.$aic)])) 
View(min.aic)

max.LL <- all.together %>%
  group_by(isolate) %>% 
  nest() %>% 
  mutate(best.model = map(data, ~.$model[which.max(.$loglik)])) 

#results numbered by set number
best.model.aic1<- data.frame(unlist(min.aic$best.model), levels(AdepidynFull$Isolate))
best.model.Log1<- data.frame(unlist(max.LL$best.model), levels(AdepidynFull$Isolate))

best.model.aic2<- data.frame(unlist(min.aic$best.model), levels(AdepidynFull$Isolate))
best.model.Log2<- data.frame(unlist(max.LL$best.model), levels(AdepidynFull$Isolate))

best.model.aic3<- data.frame(unlist(min.aic$best.model), levels(AdepidynFull$Isolate))
best.model.Log3<- data.frame(unlist(max.LL$best.model), levels(AdepidynFull$Isolate))

best.model.aic4<- data.frame(unlist(min.aic$best.model), levels(AdepidynFull$Isolate))
best.model.Log4<- data.frame(unlist(max.LL$best.model), levels(AdepidynFull$Isolate))

View(best.model.aic1)
View(best.model.aic2)




write.csv(x = besteachset,file = "besteachsetMycTrial073119.csv")

colnames(best.model.aic1) <- c("model", "isolate")
colnames(best.model.aic2) <- c("model", "isolate")
colnames(best.model.aic3) <- c("model", "isolate")
colnames(best.model.aic4) <- c("model", "isolate")

colnames(best.model.Log1) <- c("model", "isolate")
colnames(best.model.Log2) <- c("model", "isolate")
colnames(best.model.Log3) <- c("model", "isolate")
colnames(best.model.Log4) <- c("model", "isolate")

##4
set4proportions<-best.model.aic4 %>%
  group_by(model) %>%
  summarise(counts = n()) %>%
  mutate(prop = round(counts*100/sum(counts), 1))

set4proportionsLOG<-best.model.Log4 %>%
  group_by(model) %>%
  summarise(counts = n()) %>%
  mutate(prop = round(counts*100/sum(counts), 1))

##3

set3proportions<-best.model.aic3 %>%
  group_by(model) %>%
  summarise(counts = n()) %>%
  mutate(prop = round(counts*100/sum(counts), 1))

set3proportionsLOG<-best.model.Log3 %>%
  group_by(model) %>%
  summarise(counts = n()) %>%
  mutate(prop = round(counts*100/sum(counts), 1))

##2

set2proportions<-best.model.aic2 %>%
  group_by(model) %>%
  summarise(counts = n()) %>%
  mutate(prop = round(counts*100/sum(counts), 1))

set2proportionsLOG<-best.model.Log2 %>%
  group_by(model) %>%
  summarise(counts = n()) %>%
  mutate(prop = round(counts*100/sum(counts), 1))

##1

set1proportions<-best.model.aic1 %>%
  group_by(model) %>%
  summarise(counts = n()) %>%
  mutate(prop = round(counts*100/sum(counts), 1))

set1proportionsLOG<-best.model.Log1 %>%
  group_by(model) %>%
  summarise(counts = n()) %>%
  mutate(prop = round(counts*100/sum(counts), 1))


besteachset<-rbind.data.frame(best.model.aic1,best.model.aic2,best.model.aic3,best.model.aic4)
besteachsetlog<-rbind.data.frame(best.model.Log1,best.model.Log2,best.model.Log3,best.model.Log4)

View(besteachset)
#besteachset<-rbind.data.frame(best.model.aic2,best.model.aic4)
allsetsAICproportions<-rbind.data.frame( set1proportions,set2proportions, set3proportions,set4proportions)
allsetsLOGproportions<-rbind.data.frame( set1proportionsLOG,set2proportionsLOG, set3proportionsLOG,set4proportionsLOG)

allsets
View(allsetsAICproportions)

write.csv(x = besteachset,file = "allset_model.csv")

