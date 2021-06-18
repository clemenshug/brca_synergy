setwd(dirname(rstudioapi::getSourceEditorContext()$path))

file.list <- list.files(pattern='*.csv')

library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(scales)
library(tidyverse)

AllData<- lapply(file.list, read.csv)%>%
  enframe("experiments", "data")%>%
  unnest(data)


#AllData<-list(data0301_v2,data0305, data10,data01, data09)%>%
# enframe("experiments", "data")%>%
# unnest(data)

DrugConcentrations<- tribble(
  ~agent, ~concentration, 
 "Alpelisib; ", "10.0; 0.0",
 "Trametinib; ", "0.316; 0.0", 
 "Trametinib; ",  "0.3162; 0.0",
 "Dasatinib; ", "0.1; 0.0", 
 "Dasatinib; ",  "0.0999; 0.0",
 "Alpelisib; Trametinib","10.0; 0.316", 
 "Alpelisib; Trametinib","10.0; 0.316" , 
 "Alpelisib; Trametinib", "10.0; 0.3162",
 "Torin2; ", "0.0316; 0.0", 
 "Dasatinib; Torin2", "0.1; 0.0316", 
 "Dasatinib; Torin2", "0.0999; 0.0316",
 "Dasatinib; ZZ1-33B", "0.1; 0.316",
 "THZ1; ", "0.1; 0.0", 
 "Paclitaxel; THZ1",  "0.01; 0.1", 
 "Bleomycin; THZ1", "1.0; 0.1", 
 "Bleomycin; Paclitaxel", "1.0; 0.01",
 "Paclitaxel; ", "0.01; 0.0",
 "BMS-265246; ", "1.0; 0.0",
 "ZZ1-33B; ", "0.0316; 0.0",
 "BMS-265246; Cabozantinib", "1.0; 10.0", 
 "Abemaciclib; AZD7762", "1.0; 1.0",
 "AZD7762; ", "1.0; 0.0", 
 "AZD7762; LY3023414", "1.0; 0.1",
 "Cabozantinib; Paclitaxel", "10.0; 0.01", 
)


AllDataFiltered<-AllData%>%
  filter(!agent%in%DrugConcentrations$agent)%>%
  bind_rows(AllData%>%semi_join(DrugConcentrations)) %>% 
  mutate(across(c(concentration1, concentration2), round, digits=3))


AllDataFilteredMean<-AllDataFiltered%>%
  group_by(experiments, agent, cell_line, concentration, agent1, agent2, concentration1, concentration2)%>%
  summarize(across(c(cell_count, cell_count__ctrl, increase_fraction_dead, GRvalue) , mean), .groups = "drop")%>%
  filter(agent2!="YKL-5-124")%>%
  filter(agent2!="ZZ1-33B")%>%
  filter(agent1!="YKL-5-124")%>%
  filter(agent1!="ZZ1-33B")


Labels1=c("AZD7762 \n LY3023414"  ,     "Abemaciclib \n AZD7762"   ,   "Abemaciclib \n JQ1" , "BMS-265246 \n Cabozantinib",
          "Bleomycin \n Cabozantinib" , "Bleomycin \n Paclitaxel"  ,  "Bleomycin \n THZ1"  ,       
          "Cabozantinib \n Paclitaxel",  "Dasatinib \n Torin2"  ,   
          "Paclitaxel \n THZ1"    ,     "Ribociclib \n Bleomycin"  ,  "Alpelisib \n Trametinib" )



AllDataFilteredMeanCombo<-AllDataFilteredMean%>%
  filter(agent2!="")%>%mutate(agent=factor(str_replace(agent,";", " \n"), levels=Labels1))

LogKillsv1<- ggplot(AllDataFilteredMeanCombo, aes(x=agent, 
                                                  y=log10(cell_count__ctrl/cell_count), color=cell_line))+
  stat_summary(size=1.2,position=position_dodge(width = 0.2))+
  theme_bw()+
  labs(y="Log Kills", color= "Cell Line", x= "Drug Combination")+ 
  theme(axis.title.y = element_text(angle=0, vjust=0.5))
print(LogKillsv1)

ggsave("LogKills_0608.pdf", LogKillsv1, width = 14, height =7
)

AllDataFilteredMeanHSA= AllDataFilteredMeanCombo%>%
  left_join(
    AllDataFilteredMean%>%filter(agent2==""), 
    by= c("experiments","cell_line", "agent1", "concentration1" ), 
    suffix= c("", "_agent1")
  ) %>%
  left_join(
    AllDataFilteredMean%>%filter(agent2==""), 
    by= c("experiments","cell_line", "agent2"="agent1", "concentration2"="concentration1" ), 
    suffix= c("", "_agent2")
  )

AllDataFilteredMeanAdditivity= AllDataFilteredMeanHSA%>%
  mutate(agent1LogKills= log10(cell_count__ctrl_agent1/cell_count_agent1))%>%
  mutate(agent2LogKills= log10(cell_count__ctrl_agent2/cell_count_agent2))

AllDataFilteredMeanAdditivity$agent1LogKills[AllDataFilteredMeanAdditivity$agent1LogKills<0] <- 0
AllDataFilteredMeanAdditivity$agent2LogKills[AllDataFilteredMeanAdditivity$agent2LogKills<0] <- 0


AllDataFilteredMeanAdditivity=  AllDataFilteredMeanAdditivity%>%
  mutate(ComboLogKills=log10(cell_count__ctrl/cell_count))%>%
  mutate(HSA= pmax(agent1LogKills ,agent2LogKills))%>%
  mutate(Additvity= agent1LogKills + agent2LogKills)%>%
  mutate(DeltaHSA= ComboLogKills - HSA)%>%
  mutate(DeltaAdd= ComboLogKills - Additvity)


HSA<-ggplot(AllDataFilteredMeanAdditivity, aes(x=agent, y=DeltaHSA, color=cell_line))+
  stat_summary(size=1.2,position=position_dodge(width = 0.2))+
  theme_bw()+
  labs(y="Log Kills greater than HSA", color= "Cell Line", x= "Drug Combination")+ 
  theme(axis.title.y = element_text(angle=0, vjust=0.5))
print(HSA)

Bliss<-ggplot(AllDataFilteredMeanAdditivity, aes(x=agent, y=DeltaAdd, color=cell_line))+
  stat_summary(size=1.2,position=position_dodge(width = 0.2))+
  theme_bw()+
  labs(y= "Log Kills \n greater than additivity", color= "Cell Line", x= "Drug Combination")+ 
  theme(axis.title.y = element_text(angle=0, vjust=0.5))
print(Bliss)

ggsave("HSA_06082021.pdf", HSA, width = 14, height =7
)
ggsave("Bliss_06082021.pdf", Bliss, width = 14, height =7
)


AllDataFilteredMeanSummary<-AllDataFilteredMeanAdditivity%>%
  mutate(across(c(concentration1, concentration2), round, digits=3))%>%
  group_by(agent, cell_line, agent1, agent2, concentration1, concentration2)%>%
  summarize(across(c(DeltaAdd, ComboLogKills) , mean), .groups = "drop")

AllDataFilteredMeanCellLines<-AllDataFilteredMeanAdditivity%>%
  mutate(across(c(concentration1, concentration2), round, digits=3))%>%
  group_by(agent, agent1, agent2, concentration1, concentration2)%>%
  summarize(across(c(DeltaAdd, ComboLogKills) , mean), .groups = "drop")

AllDataFilteredMedianCellLines<-AllDataFilteredMeanAdditivity%>%
  mutate(across(c(concentration1, concentration2), round, digits=3))%>%
  group_by(agent, agent1, agent2, concentration1, concentration2)%>%
  summarize(across(c(DeltaAdd, ComboLogKills) , median), .groups = "drop")

AllDataFilteredMeanSummaryAT<-AllDataFilteredMeanSummary%>%
  filter(agent1=="Alpelisib")

mean(AllDataFilteredMeanSummaryAT$DeltaAdd)
1-1/10^mean(AllDataFilteredMeanSummaryAT$DeltaAdd)

#FORMULA FOR LOG-KILL CONVERSON
1-1/10^0.8 #replace 0.78 with cell fraction 
1-1/10^1.2






#IN CARBO PACLITAXEL, TAXEL HIGHEST. TRY FOR THIS COMBO, THEN DO FOR EVERY COMBO
AllDataFilteredMeanAdditivityCP<-AllDataFilteredMeanAdditivity%>%
  mutate(across(c(concentration1, concentration2), round, digits=3))%>%
  group_by(agent, cell_line, agent1, agent2, concentration1, concentration2)%>%
  summarize(across(c(DeltaAdd, ComboLogKills, agent1LogKills, agent2LogKills) , mean), .groups = "drop")%>%
  filter(agent=="Cabozantinib \n Paclitaxel")%>%
   mutate(HSA=ComboLogKills-agent2LogKills)



CP<-ggplot(data=AllDataFilteredMeanAdditivityCP, aes(x=agent1LogKills, y=HSA, group=1, color=cell_line)) +
    stat_summary(size=1.2)+
    theme_bw()+
    labs(y="Log kills greater \n than Paclitaxel", color= "Cell Line", x= "Log kills Cabozantinib")+ 
    theme(axis.title.y = element_text(angle=0, vjust=0.5))
print(CP)
ggsave("CP.pdf", CP, width = 14, height =7
)

CP2<-ggplot(data=AllDataFilteredMeanAdditivityCP, aes(x=agent2LogKills, y=HSA, group=1, color=cell_line)) +
  stat_summary(size=1.2)+
  theme_bw()+
  labs(y="Log kills greater \n than Paclitaxel", color= "Cell Line", x= "Log kills Paclitaxel")+ 
  theme(axis.title.y = element_text(angle=0, vjust=0.5))
print(CP2)

ggsave("CP2.pdf", CP2, width = 14, height =7
)

AllDataFilteredMeanAdditivityCP=AllDataFilteredMeanAdditivity%>%
  filter(agent=="Cabozantinib \n Paclitaxel")%>%
  mutate(ComboLogKills=log10(cell_count__ctrl/cell_count))%>%
  mutate(HSA=ComboLogKills-agent1LogKills)

CP3<-ggplot(data=AllDataFilteredMeanAdditivityCP, aes(x=agent1LogKills, y=HSA, group=1, color=cell_line)) +
  stat_summary(size=1.2)+
  theme_bw()+
  labs(y="Log kills greater \n than Cabozantinib", color= "Cell Line", x= "Log kills Cabozantinib")+ 
  theme(axis.title.y = element_text(angle=0, vjust=0.5))
print(CP3)
ggsave("CP3.pdf", CP3, width = 14, height =7
)

CP4<-ggplot(data=AllDataFilteredMeanAdditivityCP, aes(x=agent2LogKills, y=HSA, group=1, color=cell_line)) +
  stat_summary(size=1.2)+
  theme_bw()+
  labs(y="Log kills greater \n than Cabozantinib", color= "Cell Line", x= "Log kills Paclitaxel")+ 
  theme(axis.title.y = element_text(angle=0, vjust=0.5))
print(CP4)

ggsave("CP4.pdf", CP4, width = 14, height =7
)


AllDataFilteredMeanAdditivityAT=AllDataFilteredMeanAdditivity%>%
  filter(agent=="Alpelisib \n Trametinib")%>%
  mutate(ComboLogKills=log10(cell_count__ctrl/cell_count))%>%
  mutate(HSA=ComboLogKills-agent1LogKills)

AT<-ggplot(data=AllDataFilteredMeanAdditivityAT, aes(x=agent1LogKills, y=HSA, group=1, color=cell_line)) +
  stat_summary(size=1.2)+
  theme_bw()+
  labs(y="Log kills greater \n than Alpelisib", color= "Cell Line", x= "Log kills Alpelisib")+ 
  theme(axis.title.y = element_text(angle=0, vjust=0.5))
print(AT)
ggsave("AT.pdf", AT, width = 14, height =7
)

AT2<-ggplot(data=AllDataFilteredMeanAdditivityAT, aes(x=agent2LogKills, y=HSA, group=1, color=cell_line)) +
  stat_summary(size=1.2)+
  theme_bw()+
  labs(y="Log kills greater \n than Alpelisib", color= "Cell Line", x= "Log kills Trametinib")+ 
  theme(axis.title.y = element_text(angle=0, vjust=0.5))
print(AT2)

ggsave("AT2.pdf", AT2, width = 14, height =7
)

AllDataFilteredMeanAdditivityAT=AllDataFilteredMeanAdditivity%>%
  filter(agent=="Alpelisib \n Trametinib")%>%
  mutate(ComboLogKills=log10(cell_count__ctrl/cell_count))%>%
  mutate(HSA=ComboLogKills-agent2LogKills)

AT3<-ggplot(data=AllDataFilteredMeanAdditivityAT, aes(x=agent2LogKills, y=HSA, group=1, color=cell_line)) +
  stat_summary(size=1.2)+
  theme_bw()+
  labs(y="Log kills greater \n than Trametinib", color= "Cell Line", x= "Log kills Trametinib")+ 
  theme(axis.title.y = element_text(angle=0, vjust=0.5))
print(AT3)

ggsave("AT3.pdf", AT3, width = 14, height =7
)


AT4<-ggplot(data=AllDataFilteredMeanAdditivityAT, aes(x=agent2LogKills, y=HSA, group=1, color=cell_line)) +
  stat_summary(size=1.2)+
  theme_bw()+
  labs(y="Log kills greater \n than Trametinib", color= "Cell Line", x= "Log kills Trametinib")+ 
  theme(axis.title.y = element_text(angle=0, vjust=0.5))
print(AT4)

ggsave("AT4.pdf", AT4, width = 14, height =7
)




         

AllDataFilteredMeanHSA= AllDataFilteredMeanCombo%>%
  left_join(
    AllDataFilteredMean%>%filter(agent2==""), 
    by= c("experiments","cell_line", "agent1", "concentration1" ), 
    suffix= c("", "_agent1")
  ) %>%
  left_join(
    AllDataFilteredMean%>%filter(agent2==""), 
    by= c("experiments","cell_line", "agent2"="agent1", "concentration2"="concentration1" ), 
    suffix= c("", "_agent2")
  )

AllDataFilteredMeanAdditivity$agent1LogKills[AllDataFilteredMeanAdditivity$agent1LogKills<0] <- 0
AllDataFilteredMeanAdditivity$agent2LogKills[AllDataFilteredMeanAdditivity$agent2LogKills<0] <- 0

  
#HighestComboSASelect<-HighestComboSA%>%
#  if (which.pmax(c(LogKills_agent1, LogKills_agent2)==1)){
#    mutate(BestSA=agent1)
#  }else {
#    mutate(BestSA=agent2)
# }




#HighestComboSASelect <-HighestComboSA%>%
 # mutate(BestSA =pmax(LogKills_agent1, LogKills_agent2))%>%
  




LogKillsv1<- ggplot(AllDataFilteredMeanCombo, aes(x=agent, 
    y=log10(cell_count__ctrl/cell_count), color=cell_line))+
  stat_summary(size=1.2,position=position_dodge(width = 0.2))+
  theme_bw()+
  labs(y="Log Kills", color= "Cell Line", x= "Drug Combination")+ 
  theme(axis.title.y = element_text(angle=0, vjust=0.5))
print(LogKillsv1)

ggsave("LogKillsv1.pdf", LogKillsv1, width = 14, height =7
       )

AllDataFilteredMeanHSA= AllDataFilteredMeanCombo%>%
  left_join(
    AllDataFilteredMean%>%filter(agent2==""), 
    by= c("experiments","cell_line", "agent1", "concentration1" ), 
    suffix= c("", "_agent1")
  ) %>%
  left_join(
    AllDataFilteredMean%>%filter(agent2==""), 
    by= c("experiments","cell_line", "agent2"="agent1", "concentration2"="concentration1" ), 
    suffix= c("", "_agent2")
  )

AllDataFilteredMeanAdditivity= AllDataFilteredMeanHSA%>%
  mutate(agent1LogKills= log10(cell_count__ctrl_agent1/cell_count_agent1))%>%
  mutate(agent2LogKills= log10(cell_count__ctrl_agent2/cell_count_agent2))

AllDataFilteredMeanAdditivity$agent1LogKills[AllDataFilteredMeanAdditivity$agent1LogKills<0] <- 0
AllDataFilteredMeanAdditivity$agent2LogKills[AllDataFilteredMeanAdditivity$agent2LogKills<0] <- 0


AllDataFilteredMeanAdditivity=  AllDataFilteredMeanAdditivity%>%
  mutate(ComboLogKills=log10(cell_count__ctrl/cell_count))%>%
  mutate(HSA= pmax(agent1LogKills ,agent2LogKills))%>%
  mutate(Additvity= agent1LogKills + agent2LogKills)%>%
  mutate(DeltaHSA= ComboLogKills - HSA)%>%
  mutate(DeltaAdd= ComboLogKills - Additvity)


HSA<-ggplot(AllDataFilteredMeanAdditivity, aes(x=agent, y=DeltaHSA, color=cell_line))+
  stat_summary(size=1.2,position=position_dodge(width = 0.2))+
  theme_bw()+
  labs(y="Log Kills greater than HSA", color= "Cell Line", x= "Drug Combination")+ 
  theme(axis.title.y = element_text(angle=0, vjust=0.5))
print(HSA)

Bliss<-ggplot(AllDataFilteredMeanAdditivity, aes(x=agent, y=DeltaAdd, color=cell_line))+
  stat_summary(size=1.2,position=position_dodge(width = 0.2))+
  theme_bw()+
  labs(y= "Log Kills \n greater than additivity", color= "Cell Line", x= "Drug Combination")+ 
  theme(axis.title.y = element_text(angle=0, vjust=0.5))
print(Bliss)


AllDataFilteredMeanSummary<-AllDataFilteredMeanAdditivity%>%
  mutate(across(c(concentration1, concentration2), round, digits=3))%>%
  group_by(agent, cell_line, agent1, agent2, concentration1, concentration2)%>%
  summarize(across(c(DeltaAdd, ComboLogKills) , mean), .groups = "drop")

AllDataFilteredMeanCellLines<-AllDataFilteredMeanAdditivity%>%
  mutate(across(c(concentration1, concentration2), round, digits=3))%>%
  group_by(agent, agent1, agent2, concentration1, concentration2)%>%
  summarize(across(c(DeltaAdd, ComboLogKills) , mean), .groups = "drop")

AllDataFilteredMedianCellLines<-AllDataFilteredMeanAdditivity%>%
  mutate(across(c(concentration1, concentration2), round, digits=3))%>%
  group_by(agent, agent1, agent2, concentration1, concentration2)%>%
  summarize(across(c(DeltaAdd, ComboLogKills) , median), .groups = "drop")

AllDataFilteredMeanSummaryAT<-AllDataFilteredMeanSummary%>%
  filter(agent1=="Alpelisib")

mean(AllDataFilteredMeanSummaryAT$DeltaAdd)
1-1/10^mean(AllDataFilteredMeanSummaryAT$DeltaAdd)

#FORMULA FOR LOG-KILL CONVERSON
1-1/10^0.78 #replace 0.78 with cell fraction 





AllDataFilteredMean<-AllDataFiltered%>%
  group_by(experiments, agent, cell_line, concentration, agent1, agent2, concentration1, concentration2)%>%
  summarize(across(c(cell_count, cell_count__ctrl, increase_fraction_dead, GRvalue) , mean), .groups = "drop")

ggsave("HSAv1.pdf", HSA, width = 14, height =7
)
ggsave("Blissv1.pdf", Bliss, width = 14, height =7
)


GR<-ggplot(AllDataFilteredMeanAdditivity, aes(x=agent, y =GRvalue, color=cell_line))+
  stat_summary(position=position_dodge(width = 0.1))+
  theme_bw()+
  labs(y= "GR value", color= "Cell Line", x= "Drug Combination")+ 
  theme(axis.title.y = element_text(angle=0, vjust=0.5))+
  scale_y_reverse()


ggsave("GRv1.pdf", GR, width = 14, height =7
)

AllDataFilteredMeanComboGR<-AllDataFilteredMeanAdditivity %>%
  mutate(GRHSA=pmax(GRvalue_agent1,GRvalue_agent2)) %>%
  mutate(deltaGRHSA= GRvalue- GRHSA)


GRHSA<-ggplot(AllDataFilteredMeanComboGR, aes(x=agent, y =deltaGRHSA, color=cell_line))+
  stat_summary(position=position_dodge(width = 0.1))+
  theme_bw()+
  labs(y= "GR value \n difference from HSA", color= "Cell Line", x= "Drug Combination")+ 
  theme(axis.title.y = element_text(angle=0, vjust=0.5))+
  scale_y_reverse()

ggsave("GRDeltav1.pdf", GRHSA, width = 14, height =7
)

FractionDead<-ggplot(AllDataFilteredMeanComboGR, aes(x=agent, y =increase_fraction_dead, color=cell_line))+
  stat_summary(position=position_dodge(width = 0.1))+
  theme_bw()+
  labs(y= "Increase Fraction Dead", color= "Cell Line", x= "Drug Combination")+ 
  theme(axis.title.y = element_text(angle=0, vjust=0.5))
  
FractionDead

ggsave("IncreaseFractionDeadv1.pdf", FractionDead, width = 14, height =7
)



Drugs<- tribble(
  ~agent,
  "Alpelisib; ",
  "Trametinib; ",
  "Alpelisib; Trametinib",
)


AllDataFilteredAT<-AllData%>%
  filter(agent%in%Drugs$agent) %>% 
  bind_rows(AllData%>%semi_join(Drugs)) %>% 
  mutate(across(c(concentration1, concentration2), round, digits=3))

AllDataFilteredMeanAT<-AllDataFilteredAT%>%
  group_by(agent, cell_line, concentration, agent1, agent2, concentration1, concentration2)%>%
  summarize(across(c(cell_count, cell_count__ctrl, increase_fraction_dead, GRvalue) , mean), .groups = "drop")

Labels1=c("Alpelisib \n Trametinib")  

AllDataFilteredMeanComboAT<-AllDataFilteredMeanAT%>%
  filter(agent2!="")%>%mutate(agent=factor(str_replace(agent,";", " \n"), levels=Labels1))

AllDataFilteredMeanHSAAT= AllDataFilteredMeanComboAT%>%
  left_join(
    AllDataFilteredMeanAT%>%filter(agent2==""), 
    by= c("cell_line", "agent1", "concentration1" ), 
    suffix= c("", "_agent1")
  ) %>%
  left_join(
    AllDataFilteredMeanAT%>%filter(agent2==""), 
    by= c("cell_line", "agent2"="agent1", "concentration2"="concentration1" ), 
    suffix= c("", "_agent2")
  )

AllDataFilteredMeanAdditivityAT= AllDataFilteredMeanHSAAT%>%
  mutate(agent1LogKills= log10(cell_count__ctrl_agent1/cell_count_agent1))%>%
  mutate(agent2LogKills= log10(cell_count__ctrl_agent2/cell_count_agent2))

AllDataFilteredMeanAdditivityAT$agent1LogKills[AllDataFilteredMeanAdditivityAT$agent1LogKills<0] <- 0
AllDataFilteredMeanAdditivityAT$agent2LogKills[AllDataFilteredMeanAdditivityAT$agent2LogKills<0] <- 0


AllDataFilteredMeanAdditivityAT=  AllDataFilteredMeanAdditivityAT%>%
  mutate(ComboLogKills=log10(cell_count__ctrl/cell_count))%>%
  mutate(HSA= pmax(log10(cell_count__ctrl_agent1/cell_count_agent1),log10(cell_count__ctrl_agent2/cell_count_agent2)))%>%
  mutate(Additvity= agent1LogKills + agent2LogKills)%>%
  mutate(DeltaHSA= ComboLogKills - HSA)%>%
  mutate(DeltaAdd= ComboLogKills - Additvity)


for (i in 1:length(unique(AllDataFilteredMeanAdditivityAT$cell_line))) {

AllDataFilteredMeanAdditivityAT1<-AllDataFilteredMeanAdditivityAT%>%
  filter(cell_line==unique(AllDataFilteredMeanAdditivityAT$cell_line)[i])

AllDataFilteredMeanAdditivityAT2<-AllDataFilteredMeanAdditivityAT1%>%
  group_by(agent, cell_line, concentration, agent1, agent2, concentration1, concentration2)%>%
  summarize(across(c(DeltaAdd) , mean), .groups = "drop")

v <- ggplot(AllDataFilteredMeanAdditivityAT2, aes(x=concentration1, y=concentration2, z = DeltaAdd))+
  labs(title=unique(AllDataFilteredMeanAdditivityAT$cell_line)[i], y= "Trametinib \n Concentration", fill = "Log Kills DIfference \n from Bliss", x= "Alpelisib Concentration")+ 
  theme_bw()+
  theme(axis.title.y = element_text(angle=0, vjust=0.5))
w<-v + geom_contour_filled()
print(w)
ggsave(str_c(unique(AllDataFilteredMeanAdditivityAT$cell_line)[i], "_PI3K+MEKi.pdf"),w, width = 14, height =7)
}


Drugs<- tribble(
  ~agent,
  "Dasatinib; ",
  "Torin2; ",
  "Dasatinib; Torin2",
)



AllDataFilteredAT<-AllData%>%
  filter(agent%in%Drugs$agent) %>% 
  bind_rows(AllData%>%semi_join(Drugs)) %>% 
  mutate(across(c(concentration1, concentration2), round, digits=3))

AllDataFilteredMeanAT<-AllDataFilteredAT%>%
  group_by(agent, cell_line, concentration, agent1, agent2, concentration1, concentration2)%>%
  summarize(across(c(cell_count, cell_count__ctrl, increase_fraction_dead, GRvalue) , mean), .groups = "drop")

Labels1=c("Dasatinib \n Torin2")  

AllDataFilteredMeanComboAT<-AllDataFilteredMeanAT%>%
  filter(agent2!="")%>%mutate(agent=factor(str_replace(agent,";", " \n"), levels=Labels1))

AllDataFilteredMeanHSAAT= AllDataFilteredMeanComboAT%>%
  left_join(
    AllDataFilteredMeanAT%>%filter(agent2==""), 
    by= c("cell_line", "agent1", "concentration1" ), 
    suffix= c("", "_agent1")
  ) %>%
  left_join(
    AllDataFilteredMeanAT%>%filter(agent2==""), 
    by= c("cell_line", "agent2"="agent1", "concentration2"="concentration1" ), 
    suffix= c("", "_agent2")
  )

AllDataFilteredMeanAdditivityAT= AllDataFilteredMeanHSAAT%>%
  mutate(agent1LogKills= log10(cell_count__ctrl_agent1/cell_count_agent1))%>%
  mutate(agent2LogKills= log10(cell_count__ctrl_agent2/cell_count_agent2))

AllDataFilteredMeanAdditivityAT$agent1LogKills[AllDataFilteredMeanAdditivityAT$agent1LogKills<0] <- 0
AllDataFilteredMeanAdditivityAT$agent2LogKills[AllDataFilteredMeanAdditivityAT$agent2LogKills<0] <- 0


AllDataFilteredMeanAdditivityAT=  AllDataFilteredMeanAdditivityAT%>%
  mutate(ComboLogKills=log10(cell_count__ctrl/cell_count))%>%
  mutate(HSA= pmax(log10(cell_count__ctrl_agent1/cell_count_agent1),log10(cell_count__ctrl_agent2/cell_count_agent2)))%>%
  mutate(Additvity= agent1LogKills + agent2LogKills)%>%
  mutate(DeltaHSA= ComboLogKills - HSA)%>%
  mutate(DeltaAdd= ComboLogKills - Additvity)

##BT20, CAL851, CAL120, CAL51, BT549, HC1143 did not have full grid
for (i in 1:length(unique(AllDataFilteredMeanAdditivityAT$cell_line))) {
  
  AllDataFilteredMeanAdditivityAT1<-AllDataFilteredMeanAdditivityAT%>%
    filter(cell_line==unique(AllDataFilteredMeanAdditivityAT$cell_line)[i])
  
  AllDataFilteredMeanAdditivityAT2<-AllDataFilteredMeanAdditivityAT1%>%
    group_by(agent, cell_line, concentration, agent1, agent2, concentration1, concentration2)%>%
    summarize(across(c(DeltaAdd) , mean), .groups = "drop")
  
  v <- ggplot(AllDataFilteredMeanAdditivityAT2, aes(x=concentration1, y=concentration2, z = DeltaAdd))+
    labs(title=unique(AllDataFilteredMeanAdditivityAT$cell_line)[i], y= "Torin2 \n Concentration", fill = "Log Kills DIfference \n from Bliss", x= "Dasatinib Concentration")+ 
    theme_bw()+
    theme(axis.title.y = element_text(angle=0, vjust=0.5))
  w<-v + geom_contour_filled()
  print(w)
  ggsave(str_c(unique(AllDataFilteredMeanAdditivityAT$cell_line)[i], "_IA.pdf"),w, width = 14, height =7)
}



##Contour plot for 1 cell line, then for loop for many

##CONTOUR PLOT AND DOING DASHED LINE
##X axis, y axis, z value. stat_contour, geom_contour 
#hline/vline geom_hline(yintercept=)
