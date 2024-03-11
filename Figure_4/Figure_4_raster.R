library(tidyverse)
library(tools)
library(rstatix)
library(svglite)
library(splitstackshape)
library(ggpubr)
library(plotly)
library(processx)
library(RSelenium)


f_correspond_stx_trace <- function(x) {
  case_when(grepl("2023_09_12_0001.abf", x) ~ "2023_09_12_0004.abf",
            grepl("2023_09_12_0008.abf", x) ~  "2023_09_12_0010.abf",                
            grepl("2023_09_12_0012.abf", x) ~ "2023_09_12_0014.abf",
            grepl("2023_09_12_0015.abf", x) ~  "2023_09_12_0018.abf",
            grepl("2023_09_14_0001.abf", x) ~ "2023_09_14_0002.abf",
            grepl("2023_09_14_0003.abf", x) ~  "2023_09_14_0008.abf",
            grepl("2023_09_14_0010.abf", x) ~ "2023_09_14_0011.abf",
            grepl("2023_09_14_0012.abf", x) ~ "2023_09_14_0013.abf",
            grepl("2023_09_14_0014.abf", x) ~ "2023_09_14_0017.abf",
            grepl("2023_09_14_0021.abf", x) ~  "2023_09_14_0022.abf",                
            grepl("2023_09_14_0025.abf" , x) ~  "2023_09_14_0032.abf",
            grepl("2023_09_15_0000.abf", x) ~  "2023_09_15_0001.abf",
            grepl("2023_09_15_0003.abf", x) ~ "2023_09_15_0007.abf",
            grepl("no ap data", x) ~ "2023_09_15_0011.abf",                
            grepl("2023_09_15_0012.abf", x) ~ "2023_09_15_0013.abf",
            grepl("2023_09_18_0001.abf", x) ~  "2023_09_18_0003.abf",
            grepl("2023_09_18_0005.abf", x) ~ "2023_09_18_0008.abf",
            grepl("2023_09_18_0012.abf", x) ~ "2023_09_18_0013.abf",
            grepl("2023_09_18_0015.abf", x) ~ "2023_09_18_0017.abf",
            grepl("2023_09_18_0022.abf", x) ~ "2023_09_18_0023.abf",
            grepl("2023_09_18_0035.abf", x) ~ "2023_09_18_0036.abf",
            grepl("2023_09_18_0037.abf", x) ~ "2023_09_18_0038.abf",
            
            TRUE ~ "unknown trace")
}




fs_cells_manual <- c(
  "2023_09_12_0001",
  "2023_09_14_0012",
  "2023_09_14_0014",
  "2023_09_14_0018",
  "2023_09_14_0025",
  "2023_09_15_0012",
  "2023_09_18_0015",
  "2023_09_18_0035",
  "2023_09_18_0037"
  )


fs_cell_list <-  paste(fs_cells_manual,collapse ="|")


remove_cell <- "2021"

# Import Data files
WindowsPathName <- WindowsPathName <- ""

setwd(WindowsPathName)

UV_file_list <- dir(WindowsPathName,pattern = "*.txt")



gf_UV_data_raw<- UV_file_list %>%
  set_names(nm = (basename(.) %>% tools::file_path_sans_ext())) %>%
  map_df(read_tsv,
         col_names = T,
         .id = "FileName")
  

matched_ap_stx <- AP_characteristics %>%
  select(c(name,corresponding_Stx))%>%
  mutate(cell_type = ifelse(grepl(fs_cell_list,name),"Fast spiking","Other"))%>%
  select(-name)%>%
  rename(stx_file_name=corresponding_Stx)



gf_UV_data <- gf_UV_data_raw %>%
  rename(stx_file_name = FileNameTbl)%>%
  filter(!grepl("2021",stx_file_name))%>%
 mutate(Time_after_UV_table = AP_time_table-UV_on_time_table)%>%
 filter(between(Time_after_UV_table,-5000,20000))%>%
  right_join(.,matched_ap_stx,by="stx_file_name")

gf_UV_data_for_plot<-gf_UV_data%>%
  mutate(stx_file_name = fct_reorder(stx_file_name, (cell_type)))%>%
  arrange((cell_type))%>%
  mutate(stx_file_name = as.character(stx_file_name))%>%
  mutate(cell_number = as_factor(cumsum(stx_file_name != lag(stx_file_name, default=""))))%>%
mutate(cell_number = fct_rev(cell_number))

raster_plot <- ggplot(gf_UV_data_for_plot,aes(x=Time_after_UV_table,xend=Time_after_UV_table+20,
                                              y=factor(cell_number),
                                              yend=factor(cell_number),color=cell_type))+
  geom_segment(linewidth=10)+
  geom_vline(xintercept = 0,
             linetype="dotted", 
             color = "blue")+
  #scale_x_continuous(breaks = seq(0, 1000, by = 25))+
  theme_classic()+
  theme(axis.text.y=element_text(size=6))+
  labs(title = "Raster of AP firing following UV light exposure",x="Time after start of UV exposure (ms)",y="Cell ID",caption = "Dashed blue light shows when UV light exposure begins")

raster_plot

raster_plot_data_for_export <- gf_UV_data_for_plot%>%
  group_by(stx_file_name) %>% mutate(AP_number = row_number())%>%
  select(-c(FileName,Amplitude_table,Threshold_table,Peak_dv_dt,AP_time_table,UV_on_time_table,cell_number,AP_number_table))%>%
    pivot_wider(names_from = c(stx_file_name,cell_type),values_from = Time_after_UV_table,AP_number)

write.table(raster_plot_data_for_export, "raster_plot_data_for_publication_excel.txt",
            na = "",
            row.names = T,
            col.names = T,
            #append = TRUE,
            sep = "\t")



#ggarrange(raster_plot) %>%
ggsave(raster_plot,filename = "raster_plot.svg")
     
  for_plot <- gf_UV_data %>%
  mutate(Time_after_UV_table = (Time_after_UV_table+5000)/1000)%>%
  filter(!grepl(remove_cell,stx_file_name))

break_vals <-seq(0,25,by=1)
histo <- for_plot %>% 
  mutate(new_bin = cut(Time_after_UV_table, breaks=break_vals))%>%
  mutate(new_bin = gsub("]","",new_bin))%>%
  mutate(new_bin = gsub("(","",new_bin,fixed = T))%>%
  mutate(new_bin = gsub(",","-",new_bin,fixed = T),
         new_bin = as_factor(new_bin),
         stx_file_name = as_factor(stx_file_name))%>%
  drop_na()%>%
  as_tibble()%>%
  select(c(stx_file_name,Time_after_UV_table,new_bin))%>%
  group_by(new_bin,stx_file_name,.drop = FALSE)%>%
  summarise(number=n())%>%
  right_join(.,matched_ap_stx,by="stx_file_name")#%>%
  #mutate(cell_type = ifelse(grepl(fs_cell_list,stx_file_name),"Fast spiking","Other"))

histo_wide <- histo %>%
  pivot_wider(values_from = number,names_from = c(stx_file_name,cell_type),names_sep = ".")%>%
  as_tibble()
  

f_normalize<-function(x){
  x/mean(x[1:5])
}


nomalized_df <- histo_wide%>%
mutate(across(where(is.numeric), f_normalize))%>%
  pivot_longer(-new_bin,values_to = "norm_rate")%>%
  separate_wider_delim(name, delim = ".", names = c("trace", "file_ext","cell_type"))%>%
  select(-file_ext)

##Change histo to time axis
table_for_pop_data <- nomalized_df %>%
  cSplit("new_bin","-",type.convert="as.character")%>%
  rename(time_s = new_bin_1)%>%
  mutate(time_s = as.numeric(time_s))





corrected_type <- table_for_pop_data %>%
  mutate(trace = as_factor(trace),
         time_s = as_factor(time_s),
         cell_type = as_factor(cell_type))%>%
  select(-new_bin_2)%>%
  rename(trace_name = trace)
  

rstatix::anova_test(data=corrected_type,
                    formula = norm_rate~cell_type*time_s,
                    dv=trace_name,
                    within=c(cell_type,time_s))


normalized_plot <- ggplot(table_for_pop_data,aes(x=time_s, y=norm_rate,color=cell_type))+
  stat_summary(fun.data = "mean_se",
               geom = "errorbar",
               width = 0.5) +
  stat_summary(aes(y = norm_rate),
               fun = "mean",
               geom = "point",
               )+
stat_summary(aes(y =norm_rate),
             fun = "mean",
             geom = "line",
)+
  geom_vline(xintercept = 4,
             linetype="dotted", 
             color = "blue")+
  

  guides(x =  guide_axis(angle = 45))+
  theme_classic()+
  theme(
    plot.title=element_text(size=rel(1)),
    axis.text=element_text(size=rel(1)),
    axis.title=element_text(size=rel(1)))+
  labs(x="Time in 1 s bins",y="Average spike number\nin 1 s bin (normalized to control period)",caption = "Dashed line represents point where exposure to UV begins")
 
#coord_cartesian(xlim =c(0, 500))

normalized_plot

#reform for publication data supplement

reformed <- corrected_type %>%
  pivot_wider(names_from = c(trace_name,cell_type),values_from = norm_rate)

write.table(reformed, "population_data_time_course_publication_excel.txt",
            na = "",
            row.names = T,
            col.names = T,
            #append = TRUE,
            sep = "\t")



ggsave(normalized_plot,filename = "average_AP_frequency_vs_light_exposure.svg")
#set_palette(normalized_plot, "jco")

#raster_save_path_name = paste0(save_path_name,"scn8a_raster_varaiance.png")

#ggsave(spike_variance_plot,file = "scn8a_raster_varaiance.png", width = 6, height = 2)


#######Firing frequency vs time of final event 
histo_cum_sum <- histo%>%
  cSplit("new_bin","-",type.convert="as.character")%>%
  rename(time_s = new_bin_1)%>%
  mutate(time_s = as.numeric(time_s))%>%
  select(-starts_with("new"))%>%
  arrange(stx_file_name,-time_s)%>%
  group_by(stx_file_name) %>% 
  filter(number > 0) %>% 
  slice(1)%>%
  mutate(first_with_zero = time_s+1)


control_period_data <- gf_UV_data %>%
  filter(between(Time_after_UV_table,-5000,0))%>%
  group_by(stx_file_name,cell_type)%>%
  summarise(average_frequency = n()/5)

last_AP_time<-gf_UV_data %>%
  group_by(stx_file_name)%>%
  summarise(final_action_potential = max(Time_after_UV_table))

#latency plot for publication supplemental data

latency_supplement <- last_AP_time %>%
  left_join(.,control_period_data,by="stx_file_name")%>%
  mutate(final_action_potential = final_action_potential/1000)%>%
  select(-c(average_frequency,stx_file_name))%>%
  group_by(cell_type )%>%
  mutate(cell_number = row_number())%>%
  pivot_wider(names_from = c(cell_type),values_from = c(final_action_potential))

write.table(latency_supplement, "latency_values_supplemental_data.txt",
            na = "",
            row.names = T,
            col.names = T,
            #append = TRUE,
            sep = "\t")


Frequency

freq_last_AP_time <- control_period_data %>%
  right_join(.,last_AP_time,by = c("stx_file_name"))%>%
  mutate(final_action_potential = final_action_potential/1000)


fast_Spiking_only <- freq_last_AP_time %>%
  filter(cell_type=="Fast spiking")
rr_fs <- round(summary(lm(average_frequency ~ final_action_potential, data=fast_Spiking_only))$r.squared,2)

other_only <- freq_last_AP_time %>%
  filter(cell_type=="Other")
rr_o<-round(summary(lm(average_frequency ~ final_action_potential, data=other_only))$r.squared,2)

caption_text =  paste("Rsquared fast spiking = ",rr_fs,"\nRsquared other = ",rr_o)

last_AP_freq <- ggplot(freq_last_AP_time,aes(x=average_frequency,y=final_action_potential,color=cell_type,group=cell_type))+
  geom_point(size=5)+
  stat_smooth(method = "lm", se=FALSE,linewidth=1,linetype="dashed")+
  # stat_cor(aes(label = paste(..rr.label..)), # adds R^2 value
  #          r.accuracy = 0.01,
  #          label.x = 75, label.y = 15, size = 4)+
  theme_classic()+
  theme(
    plot.title=element_text(size=rel(2)),
    axis.text=element_text(size=rel(2)),
    axis.title=element_text(size=rel(2)))+
  labs(y="Time at which APs stop firing (s)",x = "AP frequency in control period (Hz) ", title="Control action potential frequency vs Time when APs stop firing",caption = caption_text)

last_AP_freq

ggsave(last_AP_freq,filename = "average_AP_frequency_vs_final_AP_time.svg")



#########ISI variance##########
isi_list <- gf_UV_data %>%
  filter(between(Time_after_UV_table,-5000,0))%>%
  select(stx_file_name,Time_after_UV_table)%>%
  group_by(stx_file_name)%>%
  mutate(row_id = row_number())%>%
  pivot_wider(names_from = stx_file_name, values_from = Time_after_UV_table) %>%
  select(-row_id)%>%
  as.list()

f_isi = function(x){
  x - lag(x)
  
}


isi_df <- map_df(isi_list, f_isi)%>%
  mutate(row_id = row_number())%>%
  pivot_longer(-row_id)%>%
  drop_na()%>%
  group_by(name)%>%
  summarise(variance=var(value))

##Add in 

AP_characteristics_renamed <- AP_characteristics%>%
  mutate(corresponding_Stx = f_correspond_stx_trace(name))%>%
  select(-name)%>%
  rename(name=corresponding_Stx)

AP_joined_isi <- isi_df %>%
  right_join(.,AP_characteristics_renamed,by= "name")

AP_chars_complete <- freq_last_AP_time%>%
  rename(name = stx_file_name)%>%
  right_join(.,AP_joined_isi)%>%
  #select(-c(time_s,number))%>%
  drop_na()


fast_Spiking_only <- AP_chars_complete %>%
  filter(cell_type=="Fast spiking")
rr_fs <- round(summary(lm(AP_rise_rate_max_mV_per_ms ~ final_action_potential, data=fast_Spiking_only))$r.squared,2)

other_only <- AP_chars_complete %>%
  filter(cell_type=="Other")
rr_o<-round(summary(lm(AP_rise_rate_max_mV_per_ms ~ final_action_potential, data=other_only))$r.squared,2)

caption_text_dvdt =  paste("Rsquared fast spiking = ",rr_fs,"\nRsquared other = ",rr_o)



last_ap_dv_dt <-ggplot(AP_chars_complete,aes(x=final_action_potential,y=AP_rise_rate_max_mV_per_ms,color=cell_type,group=cell_type))+
  geom_point(size=5)+
  stat_smooth(method = "lm", se=FALSE,linewidth=1,linetype="dashed")+
  theme_classic()+
  theme(
  plot.title=element_text(size=rel(2)),
  axis.text=element_text(size=rel(2)),
axis.title=element_text(size=rel(2)))+
  labs(x="Time at which APs stop firing (s)",y = "Average peak AP dv/dt", title="Average control peak AP dv/dt vs Time when APs stop firing",caption=caption_text_dvdt)

last_ap_dv_dt

ggsave(last_ap_dv_dt,filename = "average_peak_dv_dt_AP_vs_final_AP_time.svg")

#latency to final ap supplement
latency_vs_dv_dt_supplement <- AP_chars_complete %>%
  select(name,cell_type,final_action_potential,AP_rise_rate_max_mV_per_ms)%>%
  group_by(cell_type)%>%
  mutate(cell_number = row_number())%>%
  select(-name)%>%
  pivot_wider(names_from = cell_type,values_from = c(AP_rise_rate_max_mV_per_ms,final_action_potential))

write.table(latency_vs_dv_dt_supplement, "latency_values_vs_dv_dt_supplemental_data.txt",
            na = "",
            row.names = T,
            col.names = T,
            #append = TRUE,
            sep = "\t")
  








rr_fa <- round(summary(lm(AP_duration_ms~average_frequency, data=AP_chars_complete))$r.squared,2)

caption_text_fa =  paste("Rsquared  = ",rr_fa)

freq_dur_plot <-ggplot(AP_chars_complete,aes(x=AP_duration_ms,y=average_frequency,color=cell_type))+
  geom_point(size=5)+
  stat_smooth(method = "lm", se=FALSE,linewidth=1,linetype="dashed")+
  theme_classic()+
  theme(
    plot.title=element_text(size=rel(2)),
    axis.text=element_text(size=rel(2)),
    axis.title=element_text(size=rel(2)))+
  labs(x="AP Half width (ms)",y= "Control AP frequency (Hz)", title="Average frequency vs half width")#,caption=caption_text_fa)

freq_dur_plot

ggsave(freq_dur_plot,filename = "frequency_vs_half_width_color.svg")






plot_ly(AP_chars_complete, x=~AP_duration_ms, y=~AP_rise_rate_max_mV_per_ms, 
             z=~average_frequency, color=~cell_type,colors = c("red","green"),marker=list(size=8)) %>%
  add_markers(size=1)#










AP_chars_complete<-AP_chars_complete %>%
  mutate(cell_type = as_factor(cell_type))%>%
  as_tibble


str(AP_chars_complete)
         
final_ap_stats<- AP_chars_complete%>%
  t_test(final_action_potential~cell_type)%>%
add_significance() %>%
  add_xy_position(x = "cell_type",
                  scales = "free",
                  step.increase = 1
                  )
                  

final_ap_stats




AP_chars_complete


plot_this <- ggplot(AP_chars_complete, aes(cell_type, final_action_potential, group = cell_type)) +
  
  stat_summary(fun.data = "mean_se",
               geom = "errorbar",
               width = 0.5) +
  stat_summary(aes(y = final_action_potential),
               fun = "mean",
               geom = "bar",
               ) +
  geom_jitter(size=3,width = 0.2,
              ) +
  stat_pvalue_manual(
    final_ap_stats,
    tip.length = 0.05,
    hide.ns = F,
    label = "p")+
  #scale_color_manual(values = c("blue", "red")) +
  #scale_x_continuous(limits = c(1986,2014), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0,.20)))+
  theme_classic()+
  theme(legend.position = "none",
        text = element_text(size=25))+
  guides(x =  guide_axis(angle = 45))+
  labs(y="Latency to final action potential (s)",x="")

plot_this




############################################

tres<-ggboxplot(AP_chars_complete, 
          x = "cell_type",
          y = "final_action_potential",
          outlier.shape = NA)+
  geom_jitter(aes(group = cell_type), width = 0.2,size=3)+     
  stat_pvalue_manual(final_ap_stats, hide.ns = F,
                     label = "{p}")+
  labs(y= "Latency to final AP (s)")

tres

ggsave(tres,file = "latency_to_final_action .svg", width = 5, height = 8)
ggsave(tres,file = "latency_to_final_action .png", width = 5, height = 8)

 







AP_chars_complete_for_export <- AP_chars_complete %>%
  arrange(pick("name"))#%>%
  #select(-c(cell_type,first_with_zero))



write.table(AP_chars_complete_for_export,"AP_chars_table_for_cluster.txt",quote = F,row.names = F,sep = "\t")

















dd_boxplot_update <- dd_boxplot +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = 'black'))

dd_boxplot_update

AP_50_ms_save_path_name = paste0(save_path_name,"AP_count_25_ms.png")

ggarrange(dd_boxplot_update) %>%
  ggexport(filename = AP_50_ms_save_path_name,
           width = 1500,
           height = 800,
           res=300)
















#uv times

uv_times <- gf_UV_data %>%
  select(FileNameTbl,UV_on_time_table)%>%
  mutate(UV_on_time_table = UV_on_time_table/1000)%>%
  unique()

write.table(uv_times,file="uv_times.txt",quote=F,sep="\t",row.names = F)


AP_measurements <- AP_chars_complete %>% 
  select(-c(AP_overshoot_mV,AP_trough_mV))%>%
rename(peak_dv_dt = AP_rise_rate_max_mV_per_ms)%>%
pivot_longer(-c("name","cell_type"),names_to = "measurement",values_to = "value")%>%
  na.omit()%>%
  mutate(across(where(is.character), as_factor))

AP_stat_test <-AP_measurements %>%
  group_by(measurement) %>%
  t_test(value ~ cell_type) %>%
  adjust_pvalue() %>%
  #add_significance("p.adj") %>% 
  add_xy_position(x = "cell_type",scales="free",step.increase = 0.50)

  
AP_boxplot <- 
  
  ggboxplot(
    AP_measurements,
    x = "cell_type",
    y = "value",
    add = "jitter",
    color = "cell_type",
    palette = "jco",
    facet.by = "measurement",
    scales = "free",
    strip.position = "left")+
  stat_pvalue_manual(AP_stat_test, tip.length = 0, hide.ns = F,label = "p")+
  scale_y_continuous(expand = expansion(mult = c(0.05, .15)))


AP_boxplot_update <- AP_boxplot +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = 'black'))

AP_boxplot_update


ggsave(AP_boxplot_update,filename = "ap_firing_characteristics.svg",
           width = 11,
           height = 8,
         units = "in")
           


