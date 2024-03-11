#Quantification of the AP properties in the various UV conditions


library(tidyverse)
library(tools)
library(ggpubr)
library(rstatix)
library(svglite)
library(splitstackshape)

f_make_char <- function(x) {
  as.character(x)
}


setwd("/Caged_saxitoxin")

rr <- read_tsv("/advanced_spike_analysis.txt",col_names = T)%>%
  drop_na()
  

rr<-rr%>%
  cSplit("Trace_name_2_Tbl",sep = ":",drop = F)%>%
  unite(Slice_ID,
        Trace_name_2_Tbl_3,
        Trace_name_2_Tbl_4,
        sep = "_",
        remove = T)%>%
  rename(AP_threshold = Threshold_AP_amp_Tbl)%>%
  rename(AP_amplitude = AP_amp_Tbl)%>%
  rename(Peak_dv_dt = dif_amp_Tbl)%>%
  rename(p1_p2_ratio = p1_p2_Tbl)%>%
  rename(p1_amp = p1_amp_Tbl)%>%
  rename(Condition = Condition_Tbl)%>%
  rename(p1_p2_interval = p2_latency_Tbl)
  
  

selected_data <- rr %>%
  select(-c(Sweep_number_Tbl,Trace_name_1_Tbl,Trace_name_2_Tbl,Slice_Tbl,Trace_name_2_Tbl_1,Trace_name_2_Tbl_2,Trace_name_2_Tbl_5))

regular_data_subset <- selected_data %>%
  mutate_at(vars(-Slice_ID),
            f_make_char)%>%
  pivot_longer(-c(Slice_ID,Condition))
  
  

#Remove 
regular_data_subset_2 <-regular_data_subset%>%
  filter((grepl("p1_p2_interval",name)))%>%
  mutate(value = as.numeric(value))%>%
  filter(value > 2)%>%
  select(Slice_ID,name)%>%
  unique()


regular_data_subset_3 <- regular_data_subset_2%>%
  mutate(name = "p1_p2_ratio")
 
regular_data_subset_4 <- rbind(regular_data_subset_2,regular_data_subset_3)
plot_data<-anti_join(regular_data_subset,regular_data_subset_4)%>%
mutate(Condition = fct_relevel(Condition, "Control", "UV","Saxitoxin","Saxitoxin + UV"))

#plot_data<-anti_join(regular_data_subset,regular_data_subset_2)

plot_data<-anti_join(regular_data_subset,regular_data_subset_4)


plot_for_table <-plot_data%>%
  pivot_wider(names_from = c(name,Condition),values_from = value)


#--------------------------------
write.table(
  plot_for_table,
  "Sweep_AP_quantification.txt",
  na = "",
  append = FALSE,
  sep = "\t",
  dec = ".",
  row.names = F,
  col.names = TRUE
)

plot_data<-plot_data%>%
  mutate(Condition = as_factor(Condition))%>%
  mutate(Slice_ID = as_factor(Slice_ID))%>%
  mutate(name = as_factor(name))%>%
  mutate(value = as.numeric(value))
  
  
  stat.test <- plot_data %>%
  group_by(name)%>%
  t_test(value ~ Condition,ref.group = "Control")
  
  stat.test <-stat.test %>%
    add_xy_position(x = "Condition", fun = "max",scales = "free")

  
  ggboxplot(plot_data, x = "Condition", y = "value",add = "jitter") +
    facet_wrap( ~ name, scales = "free") +
    stat_pvalue_manual(stat.test, hide.ns = F,
                       label = "{p.adj}")+
    ggsave("sweep_analysis.svg",width = 30,
           height = 30,
           units = "cm")
                       






  
  
  

