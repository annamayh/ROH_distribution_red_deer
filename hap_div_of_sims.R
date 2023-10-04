#library(tidyverse)
library(tidyverse)
library(janitor)

setwd("H:/")

###################################
### first needed to thin out the snps ####
##########################################


# for (i in 1:25){
# 
#   system(paste0("plink --vcf PhD_4th_yr/Heredity_ROH_density_manuscript/hap_div_sims/Rum_neutral_output/SLiM_output_sim_",i,
#                 " --thin-count 1634 --out PhD_4th_yr/Heredity_ROH_density_manuscript/hap_div_sims/Rum_neutral_thinned/ROH_scan_sim_thinned",i,
#                 " --autosome --autosome-num 33 ",
#                 "--homozyg --homozyg-window-snp 35 --homozyg-snp 40 --homozyg-kb 2500 ",
#                 "--homozyg-density 70 --homozyg-window-missing 4 ",
#                 "--homozyg-het 0 ",
#                 "--freq --missing"))
# }
# 
# 
# for (i in 1:25){
# 
#   system(paste0("plink --vcf PhD_4th_yr/Heredity_ROH_density_manuscript/hap_div_sims/Rum_strong_sel/SLiM_output_sim_",i,
#                 " --thin-count 1634 --out PhD_4th_yr/Heredity_ROH_density_manuscript/hap_div_sims/Rum_strong_sel_thinned/ROH_scan_sim_thinned",i,
#                 " --autosome --autosome-num 33 ",
#                 "--homozyg --homozyg-window-snp 35 --homozyg-snp 40 --homozyg-kb 2500 ",
#                 "--homozyg-density 70 --homozyg-window-missing 4 ",
#                 "--homozyg-het 0 ",
#                 "--freq --missing"))
# }
# 




all_sims_hap_div_strong_sel=list()


for (sim_num in 1:23){
  

    thinned_snps=read.table(paste0("PhD_4th_yr/Heredity_ROH_density_manuscript/hap_div_sims/Rum_strong_sel_thinned/ROH_scan_sim_thinned",sim_num,".hom.summary"), header=T, stringsAsFactors = F)%>%
      select(BP)
  
  
    sim=read.table(paste0("PhD_4th_yr/Heredity_ROH_density_manuscript/hap_div_sims/Rum_strong_sel/SLiM_output_sim_",sim_num))%>%
      select(-V1,-V3,-V4, -V5, -V6, -V7, -V8, -V9)##removing other random columns you get from vcf files
    
    names(sim)[1]="BP"
    
    sims_sub=inner_join(thinned_snps,sim)%>%na.omit()%>%mutate(BP = ifelse(duplicated(BP), paste0(BP,"_2"),BP))
    
    sims_snps=sims_sub%>%
      t()%>%as_tibble()  ##trnaspose so SNPs are columns and IDs are rows
    
    haps=sims_snps%>%row_to_names(row_number = 1)
    phased_haps=haps%>%separate_rows(1:nrow(haps), sep="|") #seperating into phased haps
    phased_sims<-phased_haps[!(phased_haps[,1]=="" | phased_haps[,1]=="|"),] #removing the rows with the spaces cos didnt work properly 
    phased_lab=phased_sims
    
    ## hap div function taken from own github and modified slightly
    Haplotype_diversity <- list()
    
    for (i in seq(from = 1, to = ncol(phased_lab), by = 10)) {
          until <- ifelse((i+19) > ncol(phased_lab),  ncol(phased_lab), i+19) #10 SNP windows 
          if(until==i){break} #if window is only 1 snp long 
          window_temp <- phased_lab[, c(i:(until))]
          if (ncol(window_temp)<20){break}
          window_unite <- window_temp %>% 
            unite("string_haplo", na.rm = TRUE, remove = FALSE) %>% 
            group_by(string_haplo)%>% #grouping by unique haplotype
            tally() #counting number of unique haplotypes
          window_unite$ns<-NA
          for (j in 1:nrow(window_unite)){
            window_unite$ns[j]<-window_unite$n[j]*((window_unite$n[j])-1)
            
          }
      
      tops<-sum(window_unite$ns)
      bottom<-sum(window_unite$n)*(sum(window_unite$n)-1)
      D<-1-(tops/bottom)
      all_snps=colnames(phased_lab[i:until])## all snps in the window
      focal=matrix(c(D, all_snps), nrow=1)#combining into matrix to rbind into columns
      Haplotype_diversity <-(rbind(Haplotype_diversity,focal))
    }
    
    #######
    ## read in ROH to fund hotspots 
    ROH=read.table(paste0("PhD_4th_yr/Heredity_ROH_density_manuscript/hap_div_sims/Rum_strong_sel_thinned/ROH_scan_sim_thinned",sim_num,".hom.summary"), header=T, stringsAsFactors = F)
    hot_threshold=as.numeric(quantile(ROH$UNAFF, .99))
    ##set uo if else for yes no over threshold
    hotspot_snp_bp=ROH%>%filter(UNAFF>=hot_threshold)%>%pull(BP)%>%as.character()##filtering for only snps at hotspots
    ## need to search each row for the wheher it includes a hotspot snp or not
    hotspots_to_search=paste(hotspot_snp_bp, collapse = "|")
    
    view=as.data.frame(hotspot_snp_bp)
    
    Hap_div_df=as.data.frame(Haplotype_diversity) #duplicating to use in loop later
    Hap_div_df$hotspot=NA
    
    
    for (j in 1:nrow(Hap_div_df)){
    
        row=Hap_div_df[j,] ## pull out one rofrom haplotype div table
        match=grepl(hotspots_to_search, row)#searching whether any hotspot snps match the first row of the table
        Hap_div_df$hotspot[[j]]=ifelse(length(match[match==TRUE]), "yes","no") # adding in column whih states whether window has a hotspot snp in it 
    }
    
    
    ##making it pretty df
    pretty_hap_div=Hap_div_df%>%select(V1, hotspot)
    pretty_hap_div$window_number<-(1:nrow(pretty_hap_div))
    names(pretty_hap_div)[1]<-"D"
    pretty_hap_div$D=as.numeric(pretty_hap_div$D)
    
    all_sims_hap_div_strong_sel[[sim_num]]=pretty_hap_div

}







### scan for hap div in neutral sims now ###



all_sims_hap_div_neutral=list()


for (sim_num in 1:25){
  
   if (sim_num==7)next #'for failed sims on dud machine'
   if (sim_num==15)next
  
  thinned_snps=read.table(paste0("PhD_4th_yr/Heredity_ROH_density_manuscript/hap_div_sims/Rum_neutral_thinned/ROH_scan_sim_thinned",sim_num,".hom.summary"), header=T, stringsAsFactors = F)%>%
    select(BP)
  
  
  sim=read.table(paste0("PhD_4th_yr/Heredity_ROH_density_manuscript/hap_div_sims/Rum_neutral_output/SLiM_output_sim_",sim_num))%>%
    select(-V1,-V3,-V4, -V5, -V6, -V7, -V8, -V9)##removing other random columns you get from vcf files
  
  names(sim)[1]="BP"
  #filtering snps
  sims_sub=inner_join(thinned_snps,sim)%>%na.omit()%>%mutate(BP = ifelse(duplicated(BP), paste0(BP,"_2"),BP))
  
  
  sims_snps=sims_sub%>%
    t()%>%as_tibble()  ##trnaspose so SNPs are columns and IDs are rows
  
  haps=sims_snps%>%row_to_names(row_number = 1)
  phased_haps=haps%>%separate_rows(1:nrow(haps), sep="|") #seperating into phased haps
  phased_sims<-phased_haps[!(phased_haps[,1]=="" | phased_haps[,1]=="|"),] #removing the rows with the spaces cos didnt work properly 
  phased_lab=phased_sims
  
  ## hap div function taken from own github and modified slightly
  Haplotype_diversity <- list()
  
  for (i in seq(from = 1, to = ncol(phased_lab), by = 10)) {
    until <- ifelse((i+19) > ncol(phased_lab),  ncol(phased_lab), i+19) #10 SNP windows 
    if(until==i){break} #if window is only 1 snp long 
    window_temp <- phased_lab[, c(i:(until))]
    if (ncol(window_temp)<20){break}
    window_unite <- window_temp %>% 
      unite("string_haplo", na.rm = TRUE, remove = FALSE) %>% 
      group_by(string_haplo)%>% #grouping by unique haplotype
      tally() #counting number of unique haplotypes
    window_unite$ns<-NA
    for (j in 1:nrow(window_unite)){
      window_unite$ns[j]<-window_unite$n[j]*((window_unite$n[j])-1)
      
    }
    
    tops<-sum(window_unite$ns)
    bottom<-sum(window_unite$n)*(sum(window_unite$n)-1)
    D<-1-(tops/bottom)
    all_snps=colnames(phased_lab[i:until])## all snps in the window
    focal=matrix(c(D, all_snps), nrow=1)#combining into matrix to rbind into columns
    Haplotype_diversity <-(rbind(Haplotype_diversity,focal))
  }
  
  #######
  ## read in ROH to fund hotspots 
  ROH=read.table(paste0("PhD_4th_yr/Heredity_ROH_density_manuscript/hap_div_sims/Rum_neutral_thinned/ROH_scan_sim_thinned",sim_num,".hom.summary"), header=T, stringsAsFactors = F)
  hot_threshold=as.numeric(quantile(ROH$UNAFF, .99))
  ##set uo if else for yes no over threshold
  hotspot_snp_bp=ROH%>%filter(UNAFF>=hot_threshold)%>%pull(BP)%>%as.character()##filtering for only snps at hotspots
  ## need to search each row for the wheher it includes a hotspot snp or not
  hotspots_to_search=paste(hotspot_snp_bp, collapse = "|")
  
  view=as.data.frame(hotspot_snp_bp)
  
  Hap_div_df=as.data.frame(Haplotype_diversity) #duplicating to use in loop later
  Hap_div_df$hotspot=NA
  
  
  for (j in 1:nrow(Hap_div_df)){
    
    row=Hap_div_df[j,] ## pull out one rofrom haplotype div table
    match=grepl(hotspots_to_search, row)#searching whether any hotspot snps match the first row of the table
    Hap_div_df$hotspot[[j]]=ifelse(length(match[match==TRUE]), "yes","no") # adding in column whih states whether window has a hotspot snp in it 
  }
  
  
  ##making it pretty df
  pretty_hap_div=Hap_div_df%>%select(V1, hotspot)
  pretty_hap_div$window_number<-(1:nrow(pretty_hap_div))
  names(pretty_hap_div)[1]<-"D"
  pretty_hap_div$D=as.numeric(pretty_hap_div$D)
  
  all_sims_hap_div_neutral[[sim_num]]=pretty_hap_div
  
}









#load(file="PhD_4th_yr/Heredity_ROH_density_manuscript/hap_div_sims/hap_div.RData")
############################################################################################################################
########## plot ############################################################################################################
########################################################################################

library(wesanderson)


all_sims_unlisted_sel=bind_rows(all_sims_hap_div_strong_sel, .id = "simulation_num")

all_sims_unlisted_sel$simulation_num=as.numeric(all_sims_unlisted_sel$simulation_num)


# B=ggplot(all_sims_unlisted_sel, aes(hotspot, D, fill=hotspot)) + geom_boxplot()+
#   theme_classic()+
#   scale_x_discrete(labels = c("Remainder of genome","ROH Hotspots" ))+
#   theme(legend.position = "none", 
#         axis.title.x=element_blank())+
#   ylab("Window haplotype diversity")+
#   labs(title = "Strong selection simulation")+
#   scale_fill_brewer(palette="Dark2")+
#   ylim(0.2,1)+
#     geom_point(aes(x = hotspot, y = D, fill=factor(hotspot)), position = "jitter", alpha = 0.05)


ggplot(all_sims_unlisted_sel, aes(window_number, D)) +
  geom_line()+
  theme_bw()+
  facet_wrap(~simulation_num)+
  ylim(0.6,1)+
  labs(x = "Window", y="Haplotype diversity")+
  geom_rect(data = data.frame(simulation_num = 1), aes(xmin = 19, xmax = 25, ymin = 0.6, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  geom_rect(data = data.frame(simulation_num = 2), aes(xmin = 127, xmax = 131, ymin = 0.6, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  geom_rect(data = data.frame(simulation_num = 3), aes(xmin = 85, xmax = 88, ymin = 0.6, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  geom_rect(data = data.frame(simulation_num = 4), aes(xmin = 42, xmax = 46, ymin = 0.6, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  geom_rect(data = data.frame(simulation_num = 5), aes(xmin = 71, xmax = 73, ymin = 0.6, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  geom_rect(data = data.frame(simulation_num = 6), aes(xmin = 5, xmax = 9, ymin = 0.6, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  geom_rect(data = data.frame(simulation_num = 6), aes(xmin = 38, xmax = 40, ymin = 0.6, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  
  geom_rect(data = data.frame(simulation_num = 7), aes(xmin = 22, xmax = 26, ymin = 0.6, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  geom_rect(data = data.frame(simulation_num = 8), aes(xmin = 112, xmax = 114, ymin = 0.6, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  geom_rect(data = data.frame(simulation_num = 9), aes(xmin = 123, xmax = 127, ymin = 0.6, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  geom_rect(data = data.frame(simulation_num = 10), aes(xmin = 134, xmax = 137, ymin = 0.6, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  geom_rect(data = data.frame(simulation_num = 11), aes(xmin = 42, xmax = 43, ymin = 0.6, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  geom_rect(data = data.frame(simulation_num = 12), aes(xmin = 118, xmax = 120, ymin = 0.6, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  geom_rect(data = data.frame(simulation_num = 12), aes(xmin = 87, xmax = 88, ymin = 0.6, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  
  geom_rect(data = data.frame(simulation_num = 13), aes(xmin = 142, xmax = 144, ymin = 0.6, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  geom_rect(data = data.frame(simulation_num = 14), aes(xmin = 46, xmax = 49, ymin = 0.6, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  geom_rect(data = data.frame(simulation_num = 15), aes(xmin = 125, xmax = 129, ymin = 0.6, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  geom_rect(data = data.frame(simulation_num = 16), aes(xmin = 110, xmax = 115, ymin = 0.6, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  geom_rect(data = data.frame(simulation_num = 17), aes(xmin = 94, xmax = 99, ymin = 0.6, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  geom_rect(data = data.frame(simulation_num = 18), aes(xmin = 73, xmax = 78, ymin = 0.6, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  geom_rect(data = data.frame(simulation_num = 19), aes(xmin = 153, xmax = 158, ymin = 0.6, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  geom_rect(data = data.frame(simulation_num = 20), aes(xmin = 89, xmax = 93, ymin = 0.6, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  geom_rect(data = data.frame(simulation_num = 21), aes(xmin = 158, xmax = 161, ymin = 0.6, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  geom_rect(data = data.frame(simulation_num = 22), aes(xmin = 110, xmax = 113, ymin = 0.6, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  geom_rect(data = data.frame(simulation_num = 23), aes(xmin = 106, xmax = 109, ymin = 0.6, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)
  #geom_rect(data = data.frame(simulation_num = 24), aes(xmin = 82, xmax = 88, ymin = 0.6, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  #geom_rect(data = data.frame(simulation_num = 25), aes(xmin = 9, xmax = 11, ymin = 0.6, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  #geom_rect(data = data.frame(simulation_num = 25), aes(xmin = 136, xmax = 137, ymin = 0.6, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)
  


#####


all_sims_unlisted_neu=bind_rows(all_sims_hap_div_neutral, .id = "simulation_num")

all_sims_unlisted_neu$simulation_num=as.numeric(all_sims_unlisted_neu$simulation_num)


# A=ggplot(all_sims_unlisted_neu, aes(hotspot, D, fill=hotspot)) + geom_boxplot()+
#   theme_classic()+
#   scale_x_discrete(labels = c("Remainder of genome","ROH Hotspots" ))+
#    theme(legend.position = "none", 
#         axis.title.x=element_blank())+
#   ylab("Window haplotype diversity")+
#   labs(title = "Neutral simulation")+
#     scale_fill_brewer(palette="Dark2")+
#   ylim(0.2,1)+
#   geom_point(aes(x = hotspot, y = D, fill=factor(hotspot)), position = "jitter", alpha = 0.05)


ggplot(all_sims_unlisted_neu, aes(window_number, D)) +
  geom_line()+
  theme_bw()+
  facet_wrap(~simulation_num)+
     ylim(0.6,1)+
  labs(x = "Window", y="Haplotype diversity")+
geom_rect(data = data.frame(simulation_num = 1), aes(xmin = 32, xmax = 35, ymin = 0.7, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  geom_rect(data = data.frame(simulation_num = 2), aes(xmin = 137, xmax = 139, ymin = 0.7, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  geom_rect(data = data.frame(simulation_num = 3), aes(xmin = 81, xmax = 86, ymin = 0.7, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  geom_rect(data = data.frame(simulation_num = 4), aes(xmin = 60, xmax = 69, ymin = 0.7, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  geom_rect(data = data.frame(simulation_num = 5), aes(xmin = 30, xmax = 40, ymin = 0.7, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
  geom_rect(data = data.frame(simulation_num = 5), aes(xmin = 156, xmax = 161, ymin = 0.7, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
geom_rect(data = data.frame(simulation_num = 6), aes(xmin = 131, xmax = 133, ymin = 0.7, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
geom_rect(data = data.frame(simulation_num = 7), aes(xmin = 10, xmax = 12, ymin = 0.7, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
geom_rect(data = data.frame(simulation_num = 8), aes(xmin = 102, xmax = 105, ymin = 0.7, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
geom_rect(data = data.frame(simulation_num = 8), aes(xmin = 135, xmax = 137, ymin = 0.7, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+

geom_rect(data = data.frame(simulation_num = 9), aes(xmin = 120, xmax = 124, ymin = 0.7, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
geom_rect(data = data.frame(simulation_num = 10), aes(xmin = 81, xmax = 84, ymin = 0.7, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
geom_rect(data = data.frame(simulation_num = 10), aes(xmin = 105, xmax = 112, ymin = 0.7, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+

geom_rect(data = data.frame(simulation_num = 11), aes(xmin = 41, xmax = 43, ymin = 0.7, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
geom_rect(data = data.frame(simulation_num = 12), aes(xmin = 20, xmax = 25, ymin = 0.7, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
geom_rect(data = data.frame(simulation_num = 12), aes(xmin = 45, xmax = 49, ymin = 0.7, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+

geom_rect(data = data.frame(simulation_num = 13), aes(xmin = 79, xmax = 82, ymin = 0.7, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
geom_rect(data = data.frame(simulation_num = 14), aes(xmin = 5, xmax = 7, ymin = 0.7, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
geom_rect(data = data.frame(simulation_num = 15), aes(xmin = 28, xmax = 32, ymin = 0.7, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
geom_rect(data = data.frame(simulation_num = 15), aes(xmin = 78, xmax = 80, ymin = 0.7, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+

geom_rect(data = data.frame(simulation_num = 16), aes(xmin = 141, xmax = 142, ymin = 0.7, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
geom_rect(data = data.frame(simulation_num = 17), aes(xmin = 8, xmax = 12, ymin = 0.7, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
geom_rect(data = data.frame(simulation_num = 18), aes(xmin = 125, xmax = 132, ymin = 0.7, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
geom_rect(data = data.frame(simulation_num = 19), aes(xmin = 42, xmax = 44, ymin = 0.7, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
geom_rect(data = data.frame(simulation_num = 20), aes(xmin = 131, xmax = 134, ymin = 0.7, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
geom_rect(data = data.frame(simulation_num = 21), aes(xmin = 13, xmax = 19, ymin = 0.7, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
geom_rect(data = data.frame(simulation_num = 22), aes(xmin = 20, xmax = 22, ymin = 0.7, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)+
geom_rect(data = data.frame(simulation_num = 23), aes(xmin = 101, xmax = 104, ymin = 0.7, ymax = Inf), alpha = 0.2, fill="blue", inherit.aes = FALSE)




library(patchwork)

A+B



###############################################################################

all_sims_unlisted_neu$order=1:nrow(all_sims_unlisted_neu)


all_sims_unlisted_neu_plot=all_sims_unlisted_neu%>%filter(window_number!=1)%>%filter(window_number!=2)%>%
  filter(window_number!=3)%>%filter(window_number!=4)%>%filter(window_number!=5)%>%filter(window_number!=6)%>%#
  filter(window_number!=7)%>%filter(window_number!=8)%>%filter(window_number!=9)%>%filter(window_number!=10)

ggplot(all_sims_unlisted_neu_plot, aes(order, D, colour=hotspot, group=simulation_num)) + geom_point()+
  theme_classic()+
  coord_flip() +  # flip coordinates (puts labels on y axis)
  
  









all_sims_unlisted_sel_plot=all_sims_unlisted_sel%>%filter(window_number!=1)%>%filter(window_number!=2)%>%
  filter(window_number!=3)%>%filter(window_number!=4)%>%filter(window_number!=5)%>%filter(window_number!=6)%>%#
  filter(window_number!=7)%>%filter(window_number!=8)%>%filter(window_number!=9)%>%filter(window_number!=10)

all_sims_unlisted_sel_plot$order=1:nrow(all_sims_unlisted_sel_plot)


ggplot(all_sims_unlisted_sel_plot, aes(order, D, colour=hotspot)) + geom_point()+
  theme_classic()





hap_div_1=as.data.frame(all_sims_hap_div_neutral[1])


ggplot(hap_div_1, aes(window_number, D, colour=hotspot)) + geom_point()+
  theme_classic()










save(all_sims_unlisted_neu, all_sims_unlisted_sel, file="PhD_4th_yr/Heredity_ROH_density_manuscript/hap_div_sims/hap_div.RData" )
