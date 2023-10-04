## includes filtering out first and last 40 sNPs and density control ### 

library(plyr)
library(dplyr)
library(ggplot2)
library(windowscanr)
library(purrr)

setwd("M:/PHD_3RDYR/Data_files")
unaff <- read.table("ROH_output/ROH_search_UpdatedMb_01_2022.hom.summary", header = T, stringsAsFactors = F)
###############################################################################
#### removing the first and last 40 SNPs (i.e. min no of SNPs for a run) ####
##############################################################################
UNAFF_perchr<-list() 
for (i in 1:33){
  UNAFF_perchr[[i]]<-subset(unaff, CHR ==i)
}

remove_func<-function(x){
  n<-nrow(x)
  removed<-x[41:(n-40),]
}

removed_unaff_list<-map(UNAFF_perchr, remove_func)#using function to remove first and last 40 snps
removed<-bind_rows(removed_unaff_list, .id = "CHR")

##########################################################################
##### using window scan to check density in relation to number of ROH ####
#########################################################################
unaff_rem <- removed %>%
  mutate(KB = BP / 1000)#adding KB column to use for density control

#window scan to see how many snps there are in 1500kb windows 
pos_win <- winScan(x = unaff_rem, 
                   groups = "CHR", 
                   position = "KB", 
                   values = c("UNAFF"), 
                   win_size = 1500,
                   win_step = 100, #step every 100 kb
                   funs = c("mean"))

pos_win2<-pos_win
pos_win2$UNAFF_n <- factor(pos_win$UNAFF_n) ##only for boxplot 
ggplot(pos_win2, aes(UNAFF_n, UNAFF_mean)) +  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x="Number of SNPs in 1500Kb window", y="Mean number of ROH called in 1500Kb window", title = "Rum BP QC")+
  geom_vline(xintercept=23.5, colour="red", size=1)

cor.test(pos_win$UNAFF_n, pos_win$UNAFF_mean) ## num of ROH found was correlated with number of snps in a window 

win_filtered <- filter(pos_win, (UNAFF_n >= 23)) ##based on boxplots 20 snps was where density leveled off and was less correlated
cor.test(win_filtered$UNAFF_n, win_filtered$UNAFF_mean)
win_filtered$UNAFF_n <- factor(win_filtered$UNAFF_n) ##only for boxplot 
ggplot(win_filtered, aes(UNAFF_n, UNAFF_mean)) +  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x="Number of SNPs in 1500Kb window", y="Mean number of ROH called in 1500Kb window", title = "Rum BP Post-QC")


#############################################################################
##### getting list of all snps that pass density control ##################
###########################################################################

head(win_filtered) ## list of all windows with enough denisty 
snp_list <- list()#
counter <-0
for (i in 1:nrow(win_filtered)){
  
  counter <- counter + 1  
  
  x<-win_filtered$win_start[i] #start of window
  y<-win_filtered$win_end[i] #end of window
  z<-win_filtered$CHR[i] #chr  to prevent doubling
  
  
  KB_start <- first(which(unaff_rem$KB > x & unaff_rem$CHR == z)) #getting all snps between start and end of the window 
  KB_end <- last(which(unaff_rem$KB < y& unaff_rem$CHR == z))
  
  
  snps_between <- unaff_rem[c(KB_start:KB_end), ]
  
  snp_list[[counter]] <- snps_between #list of all snps that pass density control
  
}

snps_keep<-dplyr::bind_rows(snp_list, .id = "window_number")#unlisting
snps_filt<-snps_keep[!duplicated(snps_keep[,("SNP")]),]#filtering repeated snps


#################################################################################
#################### working out prop of ids with ROH at each snp ###############
#################################################################################

snps_filt$pall <- NA

for(i in 1:nrow(snps_filt)){
  snps_filt$pall[i] <- snps_filt$UNAFF[i]/3067
}

mean(snps_filt$pall) ##mean pall cm= 0.0577 (11/01/2021) lower than before when had less snps
sem<-sd(snps_filt$pall)/sqrt(length(snps_filt$pall)) # standard error pall
sem
quantile(snps_filt$pall, c(.01, .99)) ## 
length(which(snps_filt$pall > 0.1436224)) ##310 hotspot snps
length(which(snps_filt$pall <= 0)) ## 441 


write.table(snps_filt,
            file = "Rum_Pall_BP_filt_UPDATED_SANGER_POS.txt",
            row.names = F, quote = F, sep = "\t")
