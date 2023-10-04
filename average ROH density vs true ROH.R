library(tidyverse)

setwd("H:/") # make sure directory is the same location as the plink application #

#Average ROH density in ROH found 
roh_all <- read.table("PHD_3RDYR/Data_files/ROH_output/ROH_search_UpdatedMb_01_2022.hom", header = T, stringsAsFactors = F)
mean(roh_all$DENSITY)
## mean ROH denisty is 61.20 i.e. 1 SNP every 61.2Kb 
# = 1/0.0612 = 16.339 SNPs per 1Mb region 
# = 1634 SNPs in a 100Mb region 

#For true number of ROH in 1 simulation  
#take raw output from 1 simulation and search for ROH with default values 
system(paste0("plink --vcf PhD_4th_yr/Heredity_ROH_density_manuscript/sims_less_SNPs/Rum_neutral_SLiM_output_sim_1 ",
              "--autosome --autosome-num 33 --out PhD_4th_yr/Heredity_ROH_density_manuscript/sims_less_SNPs/ROH_scan_sim_1_full_no_mins ",
              "--homozyg  ",
              "--freq --missing"))

full_no_mins<-read.table("PhD_4th_yr/Heredity_ROH_density_manuscript/sims_less_SNPs/ROH_scan_sim_1_full_no_mins.hom.summary", header=T)%>%
  select(BP,UNAFF)%>% add_column(iteration="true")



### now to thin number of SNPs in simulations to those in eperical dataset

for (i in 1:30){
system(paste0("plink --vcf PhD_4th_yr/Heredity_ROH_density_manuscript/sims_less_SNPs/Rum_neutral_SLiM_output_sim_1 ",
              "--make-bed --thin-count 1634 --out PhD_4th_yr/Heredity_ROH_density_manuscript/sims_less_SNPs/sims_thinned_1634_its/iter",i))}

#--thin-count: 92748 variants removed (1000 remaining).


for (i in 1:30){
system(paste0("plink --bfile PhD_4th_yr/Heredity_ROH_density_manuscript/sims_less_SNPs/sims_thinned_1634_its/iter",i," ",
              "--autosome --autosome-num 33 --out PhD_4th_yr/Heredity_ROH_density_manuscript/sims_less_SNPs/sims_thinned_1634_ROH_search/ROH_out",i," ",
              "--homozyg --homozyg-window-snp 35 --homozyg-snp 40 --homozyg-kb 2500 ",
              "--homozyg-density 70 --homozyg-window-missing 4 ",
              "--homozyg-het 0 ",
              "--freq --missing"))}



thinned_ROH <- lapply(Sys.glob("PhD_4th_yr/Heredity_ROH_density_manuscript/sims_less_SNPs/sims_thinned_1634_ROH_search/ROH_out*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)

all_thinned=bind_rows(thinned_ROH, .id = "iteration")%>%select(BP,UNAFF,iteration)


thinned_plus_true=all_thinned%>%rbind(full_no_mins)



##plotting 

thinned_plus_true%>%
  mutate(highlight=ifelse(iteration=="true", "true_ROH", "thinned"))%>%
  ggplot(aes(x=BP, y=UNAFF, group=iteration,colour=highlight, size=highlight))+
  geom_line() +
  scale_color_manual(values = c("lightgrey","lightseagreen")) +
  scale_size_manual(values=c(1.3,1.5)) +
  theme(legend.position="none") +
  labs(y="ROH denisty", x="Basepair position (bp)")+
  ggtitle("Comparison of true ROH denisty (coloured) with ROH found when SNPs are thinned to the average SNP density (thinned for 30 iterations)") +
  theme_classic()+
  theme(
    legend.position="none")
  

