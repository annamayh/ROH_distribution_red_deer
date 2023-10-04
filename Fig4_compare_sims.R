library("purrr")
library(dplyr)
library(plyr)
library(ggplot2)
library(forcats)
library(wesanderson)
library(patchwork)
library(bobfunctions2)

setwd("H:/")

## function to add pall values to a list containining all iterations 
add_pall<-function(iteration_list){
  iteration_list%>%mutate(pall=UNAFF/100)
}

add_pall_7500<-function(iteration_list){
  iteration_list%>%mutate(pall=UNAFF/7500)
}


## 
add_simulation_num<-function(list){
  for (i in 1:length(list)){
    list[[i]]$simulation_num<-i} 
  
}

## function to unlist a list containing all iterations 
unlist<-function(iteration_list){bind_rows(iteration_list, .id = "simulation_num")}

## function to get the mean pall per simulation for the unlisted dataframe
unlist_then_get_mean_pall_per_sim<-function(listed_sims){
  bind_rows(listed_sims, .id = "simulation_num")%>%
    select(simulation_num, pall)%>%
    group_by(simulation_num)%>% #grouping by simulation iteration
    dplyr::summarise(Mean_pall=mean(pall))#getting mean pall of simulation
}


#func to get top 1% in a datafram
top_1_func<-function (x){quantile(x$pall, c(.01,.99))}

#function using function above to apply to all dataframes in a list of iterations
get_top1_pall<-function(iteration_list){
  
  map(iteration_list,top_1_func)%>%bind_rows(.id = "Value")%>%
    setNames(.,c("simulation_num","Bottom_1","Top_1"))
  
  
}



unlist_max_pall<-function(listed_sims){
  bind_rows(listed_sims, .id = "simulation_num")%>%
    select(simulation_num, pall)%>%
    group_by(simulation_num)%>% #grouping by simulation iteration
    dplyr::summarise(Max_pall=max(pall))#getting mean pall of simulation
}


######################################################################################################################################
######## load all simulation files in #####

## neutral model
neutral_list <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Dec_21_current_out/Rum/Neutral_re_an_noMAF.2/*/*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(neutral_list)){
  neutral_list[[i]]$simulation_num<-i #adding sim number for grouping later
} 

#selection and recombination
sel_reomb_list <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Dec_21_current_out/Rum/NewH_sel00_r/sel_r_*/ROH_out*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(sel_reomb_list)){
  sel_reomb_list[[i]]$simulation_num<-i #adding sim number for grouping later
} 

## stronger selection coefficient 
sel_s0.05_list <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Dec_21_current_out/Rum/NewH_strongsel00_r/*/ROH_out*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(sel_s0.05_list)){
  sel_s0.05_list[[i]]$simulation_num<-i #adding sim number for grouping later
}

## stronger selection coefficient 
recomb <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Dec_21_current_out/Rum/Recomb_noMAF/*/ROH_out*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(sel_s0.05_list)){
  sel_s0.05_list[[i]]$simulation_num<-i #adding sim number for grouping later
}


#############################
## servere bottleneck neutral 
##############################
bottleneck <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Dec_21_current_out/btl/Neutral_re_an_noMAF.2/*/*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(bottleneck)){
  bottleneck[[i]]$simulation_num<-i #adding sim number for grouping later
}

##bottlneck s+r
bottleneck_s_r <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Dec_21_current_out/btl/NewH_sel00_r/sel_r_*/ROH_out*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(bottleneck_s_r)){
  bottleneck_s_r[[i]]$simulation_num<-i #adding sim number for grouping later
}

## bottleneck s = 0.05
bottleneck_s0.05_r <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Dec_21_current_out/btl/NewH_strongsel00_r/*/ROH_out*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(bottleneck_s0.05_r)){
  bottleneck_s0.05_r[[i]]$simulation_num<-i #adding sim number for grouping later
}

bottleneck_r <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Dec_21_current_out/btl/recomb/*/ROH_out*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(bottleneck_s0.05_r)){
  bottleneck_s0.05_r[[i]]$simulation_num<-i #adding sim number for grouping later
}

##########################
## no bottleneck neutral##
##########################
Nobottleneck <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Dec_21_current_out/NObtl/Neutral_noMAF/*/ROH_out*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(Nobottleneck)){
  Nobottleneck[[i]]$simulation_num<-i #adding sim number for grouping later
}

##NObottlneck s+r
NObottleneck_s_r <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Dec_21_current_out/NObtl/NewH_sel00_r*/*/ROH_out*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(NObottleneck_s_r)){
  NObottleneck_s_r[[i]]$simulation_num<-i #adding sim number for grouping later
}

## NObottleneck s = 0.05

NObottleneck_s0.05_r <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Dec_21_current_out/NObtl/NewH_strongsel00_r/*/ROH_out*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(NObottleneck_s0.05_r)){
  NObottleneck_s0.05_r[[i]]$simulation_num<-i #adding sim number for grouping later
}

NObottleneck_r <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Dec_21_current_out/NObtl/Recomb_noMAF/*/ROH_out*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(NObottleneck_s0.05_r)){
  NObottleneck_s0.05_r[[i]]$simulation_num<-i #adding sim number for grouping later
}
###########################################################################################

#save.image(file="PHD_2ndYR/Slim/List_SLiM_output_new_sims.RData")

#load("PHD_2ndYR/Slim/List_all_SLiM_output.RData")

## rum pop history simulations
pall_neutral<-neutral_list%>%map(add_pall)%>%unlist_then_get_mean_pall_per_sim
pall_sel_recomb<-sel_reomb_list%>%map(add_pall)%>%unlist_then_get_mean_pall_per_sim
pall_sel_s0.05<-sel_s0.05_list%>%map(add_pall)%>%unlist_then_get_mean_pall_per_sim
pall_recomb<-recomb%>%map(add_pall)%>%unlist_then_get_mean_pall_per_sim

top1_neutral<-neutral_list%>%map(add_pall)%>%get_top1_pall()
top1_sel_recomb<-sel_reomb_list%>%map(add_pall)%>%get_top1_pall()
top1_sel_s0.05<-sel_s0.05_list%>%map(add_pall)%>%get_top1_pall()
top1_recomb<-recomb%>%map(add_pall)%>%get_top1_pall()

max_neutral<-neutral_list%>%map(add_pall)%>%unlist_max_pall()
max_sel_recomb<-sel_reomb_list%>%map(add_pall)%>%unlist_max_pall()
max_sel_s0.05<-sel_s0.05_list%>%map(add_pall)%>%unlist_max_pall()
max_recomb<-recomb%>%map(add_pall)%>%unlist_max_pall()

pall_neutral<-join(pall_neutral,top1_neutral)%>%join(max_neutral)%>%mutate(diff=Top_1-Mean_pall)
pall_sel_recomb<-join(pall_sel_recomb,top1_sel_recomb)%>%join(max_sel_recomb)%>%mutate(diff=Top_1-Mean_pall)
pall_sel_s0.05<-join(pall_sel_s0.05,top1_sel_s0.05)%>%join(max_sel_s0.05)%>%mutate(diff=Top_1-Mean_pall)
pall_recomb<-join(pall_recomb,top1_recomb)%>%join(max_recomb)%>%mutate(diff=Top_1-Mean_pall)

pall_neutral$Model<-"Neutral"
pall_sel_recomb$Model<-"Varied\nrecombination\n+ selection"
pall_sel_s0.05$Model<-"Varied\nrecombination\n+strong\nselection"
pall_recomb$Model<-"Varied\nrecombination"

pall_neutral$Pop_History<-"Rum"
pall_sel_recomb$Pop_History<-"Rum"
pall_sel_s0.05$Pop_History<-"Rum"
pall_recomb$Pop_History<-"Rum"

### bottleneck simulations
pall_bottleneck<-bottleneck%>%map(add_pall)%>%unlist_then_get_mean_pall_per_sim
pall_bottleneck_s_r<-bottleneck_s_r%>%map(add_pall)%>%unlist_then_get_mean_pall_per_sim
pall_bottleneck_s0.05_r<-bottleneck_s0.05_r%>%map(add_pall)%>%unlist_then_get_mean_pall_per_sim
pall_bottleneck_r<-bottleneck_r%>%map(add_pall)%>%unlist_then_get_mean_pall_per_sim

Btl_top1_neutral<-bottleneck%>%map(add_pall)%>%get_top1_pall()
Btl_top1_sel_recomb<-bottleneck_s_r%>%map(add_pall)%>%get_top1_pall()
Btl_top1_sel_s0.05<-bottleneck_s0.05_r%>%map(add_pall)%>%get_top1_pall()
Btl_top1_r<-bottleneck_r%>%map(add_pall)%>%get_top1_pall()

Btlmax_neutral<-bottleneck%>%map(add_pall)%>%unlist_max_pall()
Btlmax_sel_recomb<-bottleneck_s_r%>%map(add_pall)%>%unlist_max_pall()
Btlmax_sel_s0.05<-bottleneck_s0.05_r%>%map(add_pall)%>%unlist_max_pall()
Btlmax_recomb<-bottleneck_r%>%map(add_pall)%>%unlist_max_pall()

pall_bottleneck<-join(pall_bottleneck,Btl_top1_neutral)%>%join(Btlmax_neutral)%>%mutate(diff=Top_1-Mean_pall)
pall_bottleneck_s_r<-join(pall_bottleneck_s_r,Btl_top1_sel_recomb)%>%join(Btlmax_sel_recomb)%>%mutate(diff=Top_1-Mean_pall)
pall_bottleneck_s0.05_r<-join(pall_bottleneck_s0.05_r,Btl_top1_sel_s0.05)%>%join(Btlmax_sel_s0.05)%>%mutate(diff=Top_1-Mean_pall)
pall_bottleneck_r<-join(pall_bottleneck_r,Btl_top1_r)%>%join(Btlmax_recomb)%>%mutate(diff=Top_1-Mean_pall)


pall_bottleneck$Model<-"Neutral"
pall_bottleneck_s_r$Model<-"Varied\nrecombination\n+ selection"
pall_bottleneck_s0.05_r$Model<-"Varied\nrecombination\n+strong\nselection"
pall_bottleneck_r$Model<-"Varied\nrecombination"

pall_bottleneck$Pop_History<-"Severe bottleneck"
pall_bottleneck_s_r$Pop_History<-"Severe bottleneck"
pall_bottleneck_s0.05_r$Pop_History<-"Severe bottleneck"
pall_bottleneck_r$Pop_History<-"Severe bottleneck"

##### no bottleneck simulations
pall_nobtl<-Nobottleneck%>%map(add_pall_7500)%>%unlist_then_get_mean_pall_per_sim
pall_nobtl_s_r<-NObottleneck_s_r%>%map(add_pall_7500)%>%unlist_then_get_mean_pall_per_sim
pall_nobtl_s0.05_r<-NObottleneck_s0.05_r%>%map(add_pall_7500)%>%unlist_then_get_mean_pall_per_sim
pall_nobtl_r<-NObottleneck_r%>%map(add_pall_7500)%>%unlist_then_get_mean_pall_per_sim

NOBtl_top1_neutral<-Nobottleneck%>%map(add_pall_7500)%>%get_top1_pall()
NOBtl_top1_sel_recomb<-NObottleneck_s_r%>%map(add_pall_7500)%>%get_top1_pall()
NOBtl_top1_sel_s0.05<-NObottleneck_s0.05_r%>%map(add_pall_7500)%>%get_top1_pall()
NOBtl_top1_r<-NObottleneck_r%>%map(add_pall_7500)%>%get_top1_pall()

NOBtlmax_neutral<-Nobottleneck%>%map(add_pall_7500)%>%unlist_max_pall()
NOBtlmax_sel_recomb<-NObottleneck_s_r%>%map(add_pall_7500)%>%unlist_max_pall()
NOBtlmax_sel_s0.05<-NObottleneck_s0.05_r%>%map(add_pall_7500)%>%unlist_max_pall()
NOBtlmax_recomb<-NObottleneck_r%>%map(add_pall_7500)%>%unlist_max_pall()

pall_nobtl<-join(pall_nobtl,NOBtl_top1_neutral)%>%join(NOBtlmax_neutral)%>%mutate(diff=Top_1-Mean_pall)
pall_nobtl_s_r<-join(pall_nobtl_s_r,NOBtl_top1_sel_recomb)%>%join(NOBtlmax_sel_recomb)%>%mutate(diff=Top_1-Mean_pall)
pall_nobtl_s0.05_r<-join(pall_nobtl_s0.05_r,NOBtl_top1_sel_s0.05)%>%join(NOBtlmax_sel_s0.05)%>%mutate(diff=Top_1-Mean_pall)
pall_nobtl_r<-join(pall_nobtl_r,NOBtl_top1_r)%>%join(NOBtlmax_recomb)%>%mutate(diff=Top_1-Mean_pall)

pall_nobtl$Model<-"Neutral"
pall_nobtl_s_r$Model<-"Varied\nrecombination\n+ selection"
pall_nobtl_s0.05_r$Model<-"Varied\nrecombination\n+strong\nselection"
pall_nobtl_r$Model<-"Varied\nrecombination"

pall_nobtl$Pop_History<-"No bottleneck"
pall_nobtl_s_r$Pop_History<-"No bottleneck"
pall_nobtl_s0.05_r$Pop_History<-"No bottleneck"
pall_nobtl_r$Pop_History<-"No bottleneck"





plot_pall<-rbind(pall_neutral,pall_sel_recomb,pall_sel_s0.05,pall_recomb,
                 pall_bottleneck,pall_bottleneck_s_r,pall_bottleneck_s0.05_r,pall_bottleneck_r,
                 pall_nobtl,pall_nobtl_s_r,pall_nobtl_s0.05_r,pall_nobtl_r)%>%
  filter(simulation_num!=23)%>%
  filter(simulation_num!=24)%>%
  filter(simulation_num!=25)%>%
  filter(simulation_num!=26)%>%
  filter(simulation_num!=27)%>%
  filter(simulation_num!=28)
  

plot_pall<-plot_pall%>%mutate(Top1_perc=Top_1*100)

no_bottle<-rbind(pall_nobtl,pall_nobtl_s_r,pall_nobtl_s0.05_r,pall_nobtl_r)%>%mutate(Top1_perc=Top_1*100)

## This function allows us to specify which facet to annotate
## This function allows us to specify which facet to annotate
annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) 
{
  layer(data = data, stat = StatIdentity, position = PositionIdentity, 
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = TRUE, params = list(grob = grob, 
                                          xmin = xmin, xmax = xmax, 
                                          ymin = ymin, ymax = ymax))
}

###################################################################################################################################
notbl_inset<-no_bottle%>%#mutate(Top1_perc=Top_1*100)%>%
  mutate(Model = fct_relevel(Model, "Neutral","Varied\nrecombination","Varied\nrecombination\n+ selection","Varied\nrecombination\n+strong\nselection"))%>%
  ggplot(aes(x=Model, y=Top1_perc, fill=Model))+
  geom_violin(position="dodge", alpha=0.5)+
  theme_classic()+
  scale_fill_manual(values = wes_palette("Darjeeling1", n=5))+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    axis.text.x = element_blank(), 
    axis.text = element_text(size=10)
    
  )+
  theme(legend.position = "none")+
  stat_summary(fun=mean, geom="point",size=2)




top_ann_text <- data.frame(Model = "Varied\nrecombination",Top1_perc = 17,lab = "Rum actual value",
                           Pop_History = factor("No bottleneck",levels = c("Rum","Severe bottleneck","No bottleneck")))


## plotting top 1%
Top_1<-plot_pall%>%#mutate(Top1_perc=Top_1*100)%>%
  mutate(Model = fct_relevel(Model, "Neutral","Varied\nrecombination","Varied\nrecombination\n+ selection","Varied\nrecombination\n+strong\nselection"))%>%
  mutate(Pop_History = fct_relevel(Pop_History,"No bottleneck","Rum","Severe bottleneck"))%>%
  ggplot(aes(x=Model, y=Top1_perc, fill=Model))+
  geom_violin(position="dodge", alpha=0.5)+
  theme_bw()+
  scale_fill_manual(values = wes_palette("Darjeeling1", n=4))+
  xlab("Simulation Model")+
  ylab("ROH hotspot threshold per simulation\n(99th percentile ROH density)")+
  labs(tag = "A")+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=12),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"))+
  facet_wrap(~Pop_History)+
  theme(legend.position = "none")+
  stat_summary(fun=mean, geom="point",size=2)+ # adds the mean dot 
  geom_hline(yintercept = 14.9, linetype=4, color = "black", size = 1)+
  geom_text(data = top_ann_text,label = "Rum actual value")+
  annotation_custom2(ggplotGrob(notbl_inset),
                    data=data.frame(Pop_History="No bottleneck", Model="Neutral"),
                    ymin=30, ymax=80, xmin=-Inf, xmax=Inf)
 # gg_inset(ggplot2::ggplotGrob(notbl_inset), data = data.frame(Pop_History = "No bottleneck"),Model = "Varied\nrecombination",
  #         xmin = -Inf, xmax=Inf,ymin = -Inf, ymax = Inf)





####################################################################################################
max_ann_text <- data.frame(Model = "Varied\nrecombination",Max_perc = 26,lab = "Rum actual value",
                           Pop_History = factor("No bottleneck",levels = c("Rum","Severe bottleneck","No bottleneck")))


max_notbl_inset<-no_bottle%>%mutate(Max_perc=Max_pall*100)%>%
  mutate(Model = fct_relevel(Model, "Neutral","Varied\nrecombination","Varied\nrecombination\n+ selection","Varied\nrecombination\n+strong\nselection"))%>%
  ggplot(aes(x=Model, y=Max_perc, fill=Model))+
  geom_violin(position="dodge", alpha=0.5)+
  theme_classic()+
  scale_fill_manual(values = wes_palette("Darjeeling1", n=5))+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    axis.text.x = element_blank(), 
    axis.text = element_text(size=10)
    
  )+
  theme(legend.position = "none")+
  stat_summary(fun=mean, geom="point",size=2)




Max_plot<-plot_pall%>%mutate(Max_perc=Max_pall*100)%>%
  mutate(Model = fct_relevel(Model, "Neutral","Varied\nrecombination","Varied\nrecombination\n+ selection","Varied\nrecombination\n+strong\nselection"))%>%
  mutate(Pop_History = fct_relevel(Pop_History,"No bottleneck","Rum","Severe bottleneck"))%>%
  ggplot(aes(x=Model, y=Max_perc, fill=Model))+
  geom_violin(position="dodge", alpha=0.5)+
  theme_bw()+
  scale_fill_manual(values = wes_palette("Darjeeling1", n=5))+
  xlab("Simulation Model")+
  ylab("Maximum ROH density per simulation")+
  labs(tag = "B")+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=12),
    axis.text.y = element_text(face="bold"),
    axis.text.x = element_text(face="bold")
    
  )+
  facet_wrap(~Pop_History)+
  theme(legend.position = "none")+
  stat_summary(fun=mean, geom="point",size=2)+
  geom_hline(yintercept = 23.4, linetype=4, color = "black", size = 1)+
  geom_text(data = max_ann_text,label = "Rum actual value")+
  annotation_custom2(grob=ggplot2::ggplotGrob(max_notbl_inset),
                     data=data.frame(Pop_History="No bottleneck", Model="Neutral"),
                     ymin=30, ymax=80, xmin=-Inf, xmax=Inf)


Fig4=Top_1/Max_plot

ggsave("PhD_4th_yr/Heredity_ROH_density_manuscript/Main_images/Violin_ggsave.png",
       plot=Fig4, width = 15, height = 10)




#####################################################################################################

mean_ann_text <- data.frame(Model = "Varied\nrecombination\n+ selection",Mean_perc = 7.5,lab = "Rum actual value",
                            Pop_History = factor("No bottleneck",levels = c("Rum","Severe bottleneck","No bottleneck")))


mean_notbl_inset<-no_bottle%>%mutate(Mean_perc=Mean_pall*100)%>%
  mutate(Model = fct_relevel(Model, "Neutral","Varied\nrecombination","Varied\nrecombination\n+ selection","Varied\nrecombination\n+strong\nselection"))%>%
  ggplot(aes(x=Model, y=Mean_perc, fill=Model))+
  geom_violin(position="dodge", alpha=0.5)+
  theme_classic()+
  scale_fill_manual(values = wes_palette("Darjeeling1", n=5))+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    axis.text.x = element_blank(), 
    axis.text = element_text(size=10)
    
  )+
  theme(legend.position = "none")+
  stat_summary(fun=mean, geom="point",size=2)




Mean_plot<-plot_pall%>%mutate(Mean_perc=Mean_pall*100)%>%
  mutate(Model = fct_relevel(Model, "Neutral","Varied\nrecombination","Varied\nrecombination\n+ selection","Varied\nrecombination\n+strong\nselection"))%>%
  mutate(Pop_History = fct_relevel(Pop_History,"No bottleneck","Rum","Severe bottleneck"))%>%
  ggplot(aes(x=Model, y=Mean_perc, fill=Model))+
  geom_violin(position="dodge", alpha=0.5)+
  theme_bw()+
  scale_fill_manual(values = wes_palette("Darjeeling1", n=5))+
  xlab("Simulation Model")+
  ylab("Mean ROH density per simulation")+
  #labs(tag = "A")+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=10),
    axis.text.y = element_text(face="bold"),
    axis.text.x = element_text(face="bold")
    
  )+
  facet_wrap(~Pop_History)+
  theme(legend.position = "none")+
  stat_summary(fun=mean, geom="point",size=2)+
  geom_hline(yintercept = 6.4, linetype=4, color = "black", size = 1)+
  geom_text(data = mean_ann_text,label = "Rum actual value")+
  annotation_custom2(grob=ggplot2::ggplotGrob(mean_notbl_inset),
                     data=data.frame(Pop_History="No bottleneck", Model="Neutral"),
                     ymin=10, ymax=25, xmin=-Inf, xmax=Inf)


## inset plot 
no_bottle%>%
  mutate(Model = fct_relevel(Model, "Neutral","Varied\nrecombination\n+ selection","Varied\nrecombination\n+stronger\nselection"))%>%
  ggplot(aes(x=Model, y=Mean_pall, fill=Model))+
  geom_violin(position="dodge", alpha=0.5)+
  theme_classic()+
  scale_fill_manual(values = wes_palette("Darjeeling1", n=5))+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    axis.text.x = element_blank(), 
    axis.text = element_text(size=20)
    
  )+
  theme(legend.position = "none")+
  stat_summary(fun=mean, geom="point",size=2)


#####################################################################################################

#### plotting difference between top 1% and mean
plot_pall%>%
  mutate(Model = fct_relevel(Model, "Neutral","Selection","Varied\nRecomb","Varied\nrecomb\n+ selection","Higher\nselection\ncoefficients"))%>%
  mutate(Pop_History = fct_relevel(Pop_History,"Rum","Severe bottleneck","No bottleneck"))%>%
  ggplot(aes(x=Model, y=diff, fill=Model))+
  geom_violin(position="dodge", alpha=0.5)+
  theme_bw()+
  scale_fill_manual(values = wes_palette("Darjeeling1", n=5))+
  xlab("Simulation Model")+
  ylab("Difference between Top 1% and mean proportion of individuals \nwith ROH at a SNP per simulation")+
  theme(
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=10),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold")
    
  )+
  facet_wrap(~Pop_History)+
  theme(legend.position = "none")+
  stat_summary(fun=mean, geom="point",size=2)


####### STATS OF MODELS #######
library(lme4)


m5_mean<-lmer(Mean_pall~Model+(1|Pop_History), data=plot_pall)
summary(m5_mean)

coef(summary(m5_mean))
confint(m5_mean) ## can use these to get the confidence intervals 

m4_mean<-lmer(Mean_pall~Model + (1|Pop_History), data=plot_pall%>%subset(Pop_History!="Severe bottleneck"))
summary(m4_mean)
confint(m4_mean)

##max
m3_max<-lmer(Max_pall~Model + (1|Pop_History), data=plot_pall)
summary(m3_max)
confint(m3_max)

m4_max<-lmer(Max_pall~Model + (1|Pop_History), data=plot_pall%>%subset(Pop_History!="Severe bottleneck"))
summary(m4_max)
confint(m4_max)

## top1 
m3_top1<-lmer(Top_1~Model + (1|Pop_History), data=plot_pall)
summary(m3_top1)
confint(m3_top1)

m4_top1<-lmer(Top_1~Model + (1|Pop_History), data=plot_pall%>%subset(Pop_History!="Severe bottleneck"))
summary(m4_top1)
confint(m4_top1)


### how to test if sim value is within range?


### for mean pall

neu_m<-pall_neutral%>%mutate(Mean_perc=Mean_pall*100)%>%
  ggplot(aes(x=Mean_perc))+
  geom_histogram(binwidth = 1)+
  geom_vline(xintercept=11.304311, linetype=2, color = "black", size = 1)+
  geom_vline(xintercept=6.359613, linetype=2, color = "black", size = 1)+
  geom_vline(xintercept=6.4, linetype=1, color = "red", size = 1)+
  labs(title="M1: Neutral", x="Mean ROH density", tag="C")

quantile(pall_sel_recomb$Mean_pall, probs = c(0.025, 0.975))  
sel_m<-pall_sel_recomb%>%mutate(Mean_perc=Mean_pall*100)%>%
  ggplot(aes(x=Mean_perc))+
  geom_histogram(binwidth = 1)+
  geom_vline(xintercept=10.971928 , linetype=2, color = "black", size = 1)+
  geom_vline(xintercept=6.197386 , linetype=2, color = "black", size = 1)+
  geom_vline(xintercept=6.4, linetype=1, color = "red", size = 1)+
  labs(title="M3: Selection +  v recombination",x="Mean ROH density")

quantile(pall_sel_s0.05$Mean_pall, probs = c(0.025, 0.975))  
str_sel_m<-pall_sel_s0.05%>%mutate(Mean_perc=Mean_pall*100)%>%
  ggplot(aes(x=Mean_perc))+
  geom_histogram(binwidth = 1)+
  geom_vline(xintercept=13.815193 , linetype=2, color = "black", size = 1)+
  geom_vline(xintercept=8.076005 , linetype=2, color = "black", size = 1)+
  geom_vline(xintercept=6.4, linetype=1, color = "red", size = 1)+
  labs(title = "M4: Stronger selection + v recombination",x="Mean ROH density")

quantile(pall_recomb$Mean_pall, probs = c(0.025, 0.975))
recomb_m<-pall_recomb%>%mutate(Mean_perc=Mean_pall*100)%>%
  ggplot(aes(x=Mean_perc))+
  geom_histogram(binwidth = 1)+
  geom_vline(xintercept=12.411037 , linetype=2, color = "black", size = 1)+
  geom_vline(xintercept=5.933127 , linetype=2, color = "black", size = 1)+
  geom_vline(xintercept=6.4, linetype=1, color = "red", size = 1)+
  labs(title = "M2: Varied recombination",x="Mean ROH density")

neu_m+recomb_m+sel_m+str_sel_m
#

### now for threshold 
quantile(pall_neutral$Top_1, probs = c(0.025, 0.975))  
neu_t<-pall_neutral%>%mutate(top1_perc=Top_1*100)%>%
  ggplot(aes(x=top1_perc))+
  geom_histogram(binwidth = 1)+
  geom_vline(xintercept=12.525, linetype=2, color = "black", size = 1)+
  geom_vline(xintercept=28.275 , linetype=2, color = "black", size = 1)+
  geom_vline(xintercept=14.9, linetype=1, color = "red", size = 1)+
  labs(title="M1: Neutral", x = "ROH threshold (99th percentile ROH density)", tag="A")

quantile(pall_sel_recomb$Top_1, probs = c(0.025, 0.975))  
sel_t<-pall_sel_recomb%>%mutate(top1_perc=Top_1*100)%>%
  ggplot(aes(x=top1_perc))+
  geom_histogram(binwidth = 1)+
  geom_vline(xintercept=11.55  , linetype=2, color = "black", size = 1)+
  geom_vline(xintercept=22.45  , linetype=2, color = "black", size = 1)+
  geom_vline(xintercept=14.9, linetype=1, color = "red", size = 1)+
  labs(title="M3: Selection +  v recombination" , x = "ROH threshold (99th percentile ROH density)")

quantile(pall_sel_s0.05$Top_1, probs = c(0.025, 0.975))  
str_sel_t<-pall_sel_s0.05%>%mutate(top1_perc=Top_1*100)%>%
  ggplot(aes(x=top1_perc))+
  geom_histogram(binwidth = 1)+
  geom_vline(xintercept=15.675  , linetype=2, color = "black", size = 1)+
  geom_vline(xintercept=33.975  , linetype=2, color = "black", size = 1)+
  geom_vline(xintercept=14.9, linetype=1, color = "red", size = 1)+
  labs(title = "M4: Stronger selection + v recombination" , x = "ROH threshold (99th percentile ROH density)")

quantile(pall_recomb$Top_1, probs = c(0.025, 0.975))
recomb_t<-pall_recomb%>%mutate(top1_perc=Top_1*100)%>%
  ggplot(aes(x=top1_perc))+
  geom_histogram(binwidth = 1)+
  geom_vline(xintercept=12.625  , linetype=2, color = "black", size = 1)+
  geom_vline(xintercept=23.500  , linetype=2, color = "black", size = 1)+
  geom_vline(xintercept=14.9, linetype=1, color = "red", size = 1)+
  labs(title = "M2: Varied recombination" , x = "ROH threshold (99th percentile ROH density)")

neu_t+recomb_t+sel_t+str_sel_t

###### now max value ####

quantile(pall_neutral$Max_pall, probs = c(0.025, 0.975))  
neu_ma<-pall_neutral%>%mutate(max_perc=Max_pall*100)%>%
  ggplot(aes(x=max_perc))+
  geom_histogram(binwidth = 1)+
  geom_vline(xintercept=13.05 , linetype=2, color = "black", size = 1)+
  geom_vline(xintercept=29.80  , linetype=2, color = "black", size = 1)+
  geom_vline(xintercept=23.4, linetype=1, color = "red", size = 1)+
  labs(title="M1: Neutral", x="Maximum ROH denisty", tag="B")

quantile(pall_sel_recomb$Max_pall, probs = c(0.025, 0.975))  
sel_ma<-pall_sel_recomb%>%mutate(max_perc=Max_pall*100)%>%
  ggplot(aes(x=max_perc))+
  geom_histogram(binwidth = 1)+
  geom_vline(xintercept=12.55   , linetype=2, color = "black", size = 1)+
  geom_vline(xintercept=24.00   , linetype=2, color = "black", size = 1)+
  geom_vline(xintercept=23.4, linetype=1, color = "red", size = 1)+
  labs(title="M3: Selection +  v recombination",x="Maximum ROH denisty")

quantile(pall_sel_s0.05$Max_pall, probs = c(0.025, 0.975))  
str_sel_ma<-pall_sel_s0.05%>%mutate(max_perc=Max_pall*100)%>%
  ggplot(aes(x=max_perc))+
  geom_histogram(binwidth = 1)+
  geom_vline(xintercept=17    , linetype=2, color = "black", size = 1)+
  geom_vline(xintercept=39   , linetype=2, color = "black", size = 1)+
  geom_vline(xintercept=23.4, linetype=1, color = "red", size = 1)+
  labs(title = "M4: Stronger selection + v recombination", x="Maximum ROH denisty")

quantile(pall_recomb$Max_pall, probs = c(0.025, 0.975))
recomb_ma<-pall_recomb%>%mutate(max_perc=Max_pall*100)%>%
  ggplot(aes(x=max_perc))+
  geom_histogram(binwidth = 1)+
  geom_vline(xintercept=13.00   , linetype=2, color = "black", size = 1)+
  geom_vline(xintercept=26.75   , linetype=2, color = "black", size = 1)+
  geom_vline(xintercept=23.4, linetype=1, color = "red", size = 1)+
  labs(title = "M2: Varied recombination", x="Maximum ROH denisty")

neu_ma+recomb_ma+sel_ma+str_sel_ma
