## GWAS like plot of ROH density for Kintyre pop 

library(plyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(forcats)

setwd("H:/PHD_2ndYR")
CM_kint<-read.table("Recombination_ROH/Data_files/KINTYRE_Pall_CM_filt_moreSNPs.txt", header = TRUE,stringsAsFactors = FALSE)
CM_kint$CHR<-as.factor(CM_kint$CHR)
#CM_kint$CHR<-ordered(CM_kint$CHR, levels = c("5", "18", "20", "9","11","12","19","15",
 #                                          "30","21","23","1","14","33","25","13","17","29","28", "4",
  #                                         "27","22", "24","8","3","31","6","7","2","16","32","10","26"))


CM_kint<-CM_kint%>%arrange(CHR)%>%mutate(perc=pall*100)
CM_kint$Order <- 1:nrow(CM_kint)

mean(CM_kint$perc)


axis.set <- CM_kint %>% 
  dplyr::group_by(CHR) %>% 
  dplyr::summarize(center = (max(Order) + min(Order)) / 2)

CM_kint$CHR<-as.numeric(CM_kint$CHR)

CM_plot_kint<-ggplot(CM_kint, aes(Order, perc, col = as.factor(CHR%% 2), group=CHR)) +
  scale_colour_manual(values = c("coral", "cornflowerblue")) +
  geom_point() +
  theme_classic()+
  labs(y="ROH density \n(% of ids with a ROH)", title="Using genetic map positions",x="Chromosome")+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 8, vjust = 0.5),
        axis.text.y = element_text(size = 12, vjust = 0.5)) +
  geom_hline(yintercept = 8.28, linetype="dashed", color = "red") +
  geom_hline(yintercept = 3.4, colour = "red")+
  scale_x_continuous(limits = c(0,25194), expand = c(0, 0), label = axis.set$CHR, breaks = axis.set$center)+
  scale_y_continuous(limits = c(0,20))+
  geom_line(size=1)


quantile(CM_kint$perc, c(.01, .99)) ## 
length(which(CM_kint$perc > 8.280255)) ##310 hotspot snps

cm_hotspots_kint<-CM_kint%>%filter(perc> 8.280255)



###### BP genomwide plot using same script ################

BP_kint<-read.table("Recombination_ROH/Data_files/Kintyre_BP_UpdatedMbSANGER_filtSNPS.txt", header = TRUE,stringsAsFactors = FALSE)
BP_kint$CHR<-as.factor(BP_kint$CHR)
#BP_kint$CHR<-ordered(BP_kint$CHR, levels = c("5", "18", "20", "9","11","12","19","15",
 #                                          "30","21","23","1","14","33","25","13","17","29","28", "4",
  #                                         "27","22", "24","8","3","31","6","7","2","16","32","10","26"))


BP_kint<-BP_kint%>%arrange(CHR)%>%mutate(perc=pall*100)
BP_kint$Order <- 1:nrow(BP_kint)


axis.set_BP <- BP_kint %>% 
  dplyr::group_by(CHR) %>% 
  dplyr::summarize(center = (max(Order) + min(Order)) / 2)

BP_kint$CHR<-as.numeric(BP_kint$CHR)

BP_plot_kint<-ggplot(BP_kint, aes(Order, perc, col = as.factor(CHR%% 2), group=CHR)) +
  scale_colour_manual(values = c("coral", "cornflowerblue")) +
  geom_point() +
  theme_classic()+
  labs(tag = "B", y="ROH density \n(% of ids with a ROH)", title = "Using physical map positions")+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 8, vjust = 0.5),
        axis.text.y = element_text(size = 12, vjust = 0.5),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 8.92, linetype="dashed", color = "red") +
  geom_hline(yintercept = 3.2, colour = "red")+
  scale_x_continuous(limits = c(0,28055), expand = c(0, 0), label = axis.set_BP$CHR, breaks = axis.set_BP$center)+
  scale_y_continuous(limits = c(0,20))+
  geom_line(size=1)


quantile(BP_kint$perc, c(.01, .99)) ## 
length(which(BP_kint$perc > 8.917197 )) ##310 hotspot snps

bp_hotspots_kint<-BP_kint%>%filter(perc> 8.917197 )


B<-BP_plot_kint+CM_plot_kint+ plot_layout(ncol=1)
B



ggsave(file="PHD_2ndYR/Manuscript_draft/images/GWAS_Kintyre_31.png",plot=B)





rum<-cm_hotspots%>%select(CHR, SNP, perc, BP)%>%rename(perc_cm=perc, bp_cm=BP)%>%join(bp_hotspots)%>%na.omit

BP_kint$BP[BP_kint$SNP=="cela1_red_8_26504836"]

