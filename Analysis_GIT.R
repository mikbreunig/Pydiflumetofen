

##Figures, tables, analysis, and statistical tests for pydiflumetofen baseline sensitivity paper 

ipak <- function( pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[,"Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("broom","ggmap","sf","stats","ggrepel","pwr","stringr", "tidyverse","ggplot2", "drc","lattice","car", "lme4", "lsmeans", "plyr", "plotrix", "knitr", "ggplot2", "lmtest", "lmerTest", "Rmisc", "gridExtra", "plotly", "webshot", "ggpmisc", "ggsci","scales")
ipak(packages)

ls()
rm(list = ls(all=TRUE))

EC50SPORE<-read_delim(file = "SummaryEC50SPORE_LL3.csv",delim=",",na=".")
EC50Mycelial<-read_delim(file = "SummaryEC50Myc_LL3.csv",delim=",",na=".")
View(EC50SPORE)
SampleInfo<-read_delim(file = "IsolateMetaData_073119.csv",delim=",",na=".")

GeoData<-read_delim(file="Collection_DataCSV_final.csv",delim=",",na=".")

View(EC50SPORE)

SPORE50<-as.tibble(EC50SPORE) %>% 
  mutate_at(vars(-Isolate), funs(round(., 3))) %>% 
  select(-X1) %>%
  unite("C.I",CI.min, CI.max, sep=" - ")

#View(SPORE50)
Myc50<-as.tibble(EC50Mycelial) %>% 
  select(-X1) %>% 
 mutate_at(vars(-Isolate), funs(round(., 3))) %>% 
  unite("C.I",CI.min, CI.max, sep=" - ")
View(Myc50)

Joined<-SPORE50 %>% 
  left_join(Myc50,by="Isolate") %>% 
  mutate(difference=Mean.y-Mean.x)

### to get meta data with myc 
MycTable<-Myc50 %>% 
  left_join(Meta,by="Isolate")


#Summary statistics 
summary(SPORE50$Mean)
summary(Myc50$Mean)




###Visualize trends (or lack of) by factor
MycTable<-Myc50 %>% 
  left_join(Meta,by="Isolate")
View(MycTable)
host<-ggplot(MycTable, aes(x=Mean,fill=factor(MycTable$host)))+
  geom_histogram()+ 
  labs(x=expression("Mean EC"[50]*" Estimation (µg/ml)"))+
  guides(fill=guide_legend(title=NULL))
year<-ggplot(MycTable, aes(x=Mean,fill=factor(MycTable$year)))+
  geom_histogram()+
  labs(x=expression("Mean EC"[50]*" Estimation (µg/ml)"))+
  guides(fill=guide_legend(title=NULL))
  
county<-ggplot(MycTable, aes(x=Mean,fill=factor(MycTable$county)))+
  geom_histogram()+
  labs(x=expression("Mean EC"[50]*" Estimation (µg/ml)"))+
  guides(fill=guide_legend(title=NULL))

ggsave(host, filename = "host.png",  bg = "transparent", dpi=600,units="in", height=5, width=7)
ggsave(year, filename = "year.png",  bg = "transparent", dpi=600,units="in", height=5, width=7)
ggsave(county, filename = "county.png",  bg = "transparent", dpi=600,units="in", height=5, width=7)


###Shapiro-Wilk Test for normality 
shapiro.test(Myc50$Mean)
shapiro.test(SPORE50$Mean)


##correlation between 

cor(log(Joined$Mean.x),log(Joined$Mean.y),method="spearman")
cor.test(Joined$Mean.x, Joined$Mean.y,method="spearman")
t.test(Joined$Mean.x, Joined$Mean.y, paired = TRUE, alternative = "two.sided")

###Figure correlation between spore germ and myc
cor<-ggplot(data=Joined, aes(x=Mean.x, y=Mean.y))+
  geom_point() +
  #geom_smooth(method="lm",se=FALSE,color="dodgerblue2")+
  ylab("Mycelial Growth EC50")+
  xlab("Spore Germination EC50")+
  #xlim(0,0.2)+
  #geom_smooth(method="lm")
  theme(plot.margin = unit(c(6,6,6,6),"mm"),
        plot.title = element_text(size=20, face="bold"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

ggsave(cor, filename = "corrnoline.png",  bg = "transparent", dpi=600,units="in", height=5, width=7)


######FIGURES


plot1<-hist(Myc50$Mean,
     breaks=c(0,seq(0.0125,0.30,0.0125)), 
     xlim=c(0,0.3), ylim=c(0,25),
     xaxs="i", yaxs="i", las=1, 
     main="", xlab=expression("Mean EC"[50]*" Estimation (µg/ml)"), ylab="Number of isolates")
     
plot1<-abline(v=0.060319 ,lty=2,col="Black")
plot1<-text(x=0.235,y=20,"Mycelial Growth Assay\n Population Mean  = 0.063 \n N = 94", cex=1)



plot2<-hist(SPORE50$Mean,
            breaks=c(0,seq(0.02,1.2,0.02)), 
            xlim=c(0,0.7), ylim=c(0,10),
            xaxs="i", yaxs="i", las=1, 
            main="", xlab=expression("Mean EC"[50]*" Estimation (µg/ml)"), ylab="Number of isolates")

plot2<-abline(v=0.321 ,lty=2,col="Black")
plot2<-text(x=0.5,y=7,"Spore Germination Assay\n Population Mean  = 0.321 \n N = 21", cex=1)


ggsave(plot1, filename = "MycHisto_030620.png",  bg = "transparent", dpi=600,units="mm", height=60, width=90)
ggsave(plot2, filename = "SPOREHisto_030620.png",  bg = "transparent", dpi=600,units="mm", height=60, width=90)





#### example of growth curve 
MycelialFull<-read_delim(file = "Mycelial_processed.csv",delim=",",na=".")
SporeFull<-read_delim(file = "Spore_processed.csv",delim=",",na=".")
MycelialFull$relative<-as.numeric(MycelialFull$relative)
SporeFull$relative<-as.numeric(SporeFull$relative)

spore.ex<-SporeFull %>% 
  filter(Isolate=="C5-2D") 
Myc.ex<-MycelialFull %>% 
  filter(Isolate=="C5-2D")


curve.exMyc<-drm(100*relative~ ppm, data = Myc.ex, fct = LL.3())
curve.exSPore<-drm(100*relative~ ppm, data = spore.ex, fct = LL.3())

plot(curve.exSPore,xlab="Pydiflumetofen Concentration (µg/ml)",ylab="Relative Growth (%)", lty=1, xlim=c(0,5))
plot(curve.exMyc,add=T,lty=2)
legend("topright", legend = c("Spore Germination", "Mycelial Growth"), col = c("black", "black"), lty = 1:2, cex = 0.8)
     
    

####Power Analysis
p.out <-pwr.p.test(h = ES.h(p1 = 0.01, p2 = 0.0),
                   sig.level = 0.05, 
                   n = 94,
                   alternative = "greater")

p.out2 <-pwr.p.test(h = ES.h(p1 = 0.05, p2 = 0.0),
                   sig.level = 0.05, 
                   n = 94,
                   alternative = "greater")


##META-Data and maps 

##Merge all the metadata
View(metadata)
metadata<-left_join(SampleInfo,GeoData,by="Field") 
metadata$long<-as.double(metadata$long)
metadata$lat<-as.double(metadata$lat)
Counties<-read_delim(file = "CountyGPS.csv",delim=",",na=".")
AllMeta<-left_join(metadata,Counties,by=c("long","lat"))
Meta<-distinct(AllMeta,.keep_all = TRUE)

##st.as.sf won't take empty fields so have to get ride of few isolates I don't have coordinates for
metadata2<-left_join(Myc50,metadata,by="Isolate") %>% 
  filter(Isolate!="C17A") %>% 
  filter(Isolate!="ph-1") %>% 
  filter(Isolate!="F_15_182")

#Find corresponding counties to the GPS coordinates 
library(maps) #must do this so uses right maps
US <- st_as_sf(map("county", plot = FALSE, fill = TRUE))##turn states map into object
neededpoints<-data.frame(x=metadata2$long, y=metadata2$lat)
neededpoints <- st_as_sf(neededpoints, coords = c("x", "y"), crs = st_crs(US))
finished<-st_join(neededpoints, US)
write.csv(x=finished, file="CountyGPS.csv")
#Edited filed by hand to separate long/lat etc 

###MAKE MAP
as.numeric(MycTable$long) %>% 
  as.numeric(MycTable$lat)
edited$long<-as.numeric(as.character(edited$long))

class(edited$lat)
as.factor(edited$year)
Mi<-map_data("county","michigan")

map<-ggplot(data = Mi) + 
  geom_polygon(aes(x = long, y = lat,group = group),fill="grey76", color = "white") + 
  coord_fixed(1.3) +
  guides(fill=FALSE)+
  geom_point(data=MycTable, aes(x=MycTable$long, y=MycTable$lat))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggsave(map, filename = "map.png",  bg = "transparent", dpi=600,units="in", height=5, width=7)

