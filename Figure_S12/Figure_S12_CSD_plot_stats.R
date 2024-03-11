library(tidyverse)
library(splitstackshape)
library (data.table)
library(rstatix)
library(ggpubr)
library(svglite)


WindowsPathName <- ""

setwd(WindowsPathName )
tps_files <- dir(WindowsPathName,pattern = ".*.txt")


f_import <- function(x){read_tsv(x,
                                 col_names = T,
                                 #skip = ,
                                 id="FileName")%>%
              pivot_longer(-c(FileName,csd_cond_table))
    }

imported_data<-tps_files%>%
  map_df(.,f_import)

new_order <- c("000UV","100UV","000UVstx","020UVstx","040UVstx","060UVstx","080UVstx","100UVstx1","100UVstx2","100UVstx3","100UVstx4","100UVstx5","100UVstx6","100UVstx7","100UVstx8","100UVstx9","15wash","30wash")  


selected_measurements <- imported_data %>%
  #filter(grepl("total_vol_table",name))%>%
  cSplit("csd_cond_table","_")%>%
  rename(condition = csd_cond_table_2,
         measurement = name)%>%
  cSplit("FileName","_")%>%
  rename(slice_id = FileName_1)%>%
  select(measurement,value,condition,slice_id)%>%
  mutate(condition = factor(condition, levels=new_order))

ggplot(selected_measurements,aes(x=condition,y=value,group=1))+
  geom_point(aes(color=measurement))+
  facet_wrap(~slice_id)+
  #stat_summary(fun = "mean", colour = "red", size = 2, geom = "point")+
  scale_x_discrete(guide = guide_axis(angle = 90))


#220406B remove
conditions_to_remove<-c("100UVstx2","100UVstx3","100UVstx4")
slice_to_remove <- "220406B"

filtered_data <- selected_measurements%>%
  filter(!((condition %in% conditions_to_remove) & (slice_id %in% slice_to_remove)))%>%
  filter(!(grepl("stx6|stx7|stx8|stx9|30wash",condition)))%>%
  filter(!((condition %in% "15wash") & (slice_id %in% "220406C")))%>%
  filter(!(grepl("220408C",slice_id)))
 
ggplot(filtered_data,aes(x=condition,y=value,group=1))+
  geom_point()+
  facet_wrap(~slice_id)+
  scale_x_discrete(guide = guide_axis(angle = 90)) 


f_normalize <- function(x){
  x/x[1]
}


normalized_data <- filtered_data%>%
  pivot_wider(names_from = c(slice_id,measurement),values_from = value,names_sep = ".")%>%
  arrange(condition)%>%
  mutate(across(starts_with("22"),f_normalize))%>%
  #select(-measurement)%>%
  pivot_longer(-condition,names_to = "slice_id")%>%
  cSplit("slice_id",".")%>%
  rename(slice_id = slice_id_1,
         measurement = slice_id_2)%>%
  filter(grepl("csd_source_vol_table|csd_total_postsynaptic_vol_table|csd_total_presynaptic_vol_table",measurement))


  
 
ggplot(normalized_data,aes(x=condition,y=value,group=4))+
  geom_line(aes(color=slice_id))+
  facet_wrap(~measurement)+
  #stat_summary(fun = "mean", colour = "red", size = 2, geom = "point")+
  scale_x_discrete(guide = guide_axis(angle = 90))

table <- normalized_data %>%
  group_by(condition,measurement) %>%
  summarise(avg=mean(value,na.rm = TRUE), std_dev =sd(value,na.rm = TRUE),number = n())%>%
  mutate(std_err = std_dev/sqrt(number-1))#%>%
  filter(grepl("csd_total",measurement))
  

ggplot(table,aes(x=condition,y=avg,group=measurement))+
  geom_line(aes(color=measurement))+
  geom_errorbar(aes(ymin=avg-std_err, ymax=avg+std_err,color=measurement), width=.3)+
  geom_point(size=2,aes(color=measurement))+
  theme_classic()
  #theme(legend.position="none")
  #facet_wrap(~measurement)
ggsave("pop_data2.svg",width = 3,height = 3)


#source sink etc

all_measurements <- imported_data %>%
  filter(grepl("total_vol_table",name))%>%
  cSplit("csd_cond_table","_")%>%
  rename(condition = csd_cond_table_2,
         measurement = name)%>%
  cSplit("FileName","_")%>%
  rename(slice_id = FileName_1)%>%
  select(measurement,value,condition,slice_id)%>%
  mutate(condition = factor(condition, levels=new_order))%>%
  cSplit("measurement","_")%>%
  rename(measurement = measurement_2)%>%
  select(!(starts_with("measurement_")))
  

filtered_all_measurements <- all_measurements%>%
  filter(!((condition %in% conditions_to_remove) & (slice_id %in% slice_to_remove)))%>%
  #filter(!(grepl("stx6|stx7|stx8|stx9|30wash",condition)))%>%
  filter(!(grepl("220408C|220408A",slice_id)))%>%
  filter(!((condition %in% "15wash") & (slice_id %in% "220406C")))%>%
  filter(grepl("000UVstx|100UVstx5|15wash",condition))


anova <- filtered_all_measurements %>% 
  anova_test(value ~ condition)


stat.test <- filtered_all_measurements %>%
  tukey_hsd(value ~ condition) 


ggboxplot(filtered_all_measurements, x = "condition", y = "value",add = "jitter")+
  stat_pvalue_manual(
    stat.test, label = "p.adj", 
    y.position = c(70, 65, 60)
  )

ggsave("plot_with_stats.svg",width = 4,height = 4)


ggplot(filtered_all_measurements,aes(condition,value))+
         geom_boxplot(outlier.shape = NA)+
  #geom_text(aes(label=slice_id))+
  geom_jitter(width = 0.2)+
           #facet_wrap(~condition,scales = "fixed")+
  theme_classic()
ggsave("quantification.svg",width = 5,height = 5)



for_stats <- normalized_data %>%
  filter(grepl("stx2|stx3|stx4",condition))%>%
  drop_na()%>%
  filter(!grepl("total_vol",measurement))


anova <- for_stats %>% 
  anova_test(value ~ condition*measurement)

stat.test <- for_stats %>%
  group_by(condition) %>%
  t_test(value ~ measurement) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test



imported_data%>%
  group_by(name,FileName)%>%
  summarise(count = n())


all_measurements <- imported_data %>%
  #filter(grepl("total_vol_table",name))%>%
  cSplit("csd_cond_table","_")%>%
  rename(condition = csd_cond_table_2,
         measurement = name)%>%
  cSplit("FileName","_")%>%
  rename(slice_id = FileName_1)%>%
  select(measurement,value,condition,slice_id)%>%
  mutate(condition = factor(condition, levels=new_order))%>%
  cSplit("measurement","_")%>%
  unite(measurement,measurement_2,measurement_3)%>%
  select(!(starts_with("measurement_")))


filtered_all_measurements <- all_measurements%>%
  filter(!((condition %in% conditions_to_remove) & (slice_id %in% slice_to_remove)))%>%
  filter(!(grepl("220408C|220408A",slice_id)))%>%
  filter(!((condition %in% "15wash") & (slice_id %in% "220406C")))%>%
  filter(grepl("000UVstx|100UVstx5|15wash",condition))%>%
  filter(grepl("min_table|sink_vol|source_vol|total_vol",measurement))

for_supplemental_data <- filtered_all_measurements%>%
  pivot_wider(names_from = measurement, values_from = value)

write.table(for_supplemental_data, "csd_scatter_supplemental_data.dat",
            na = "",
            row.names = T,
            col.names = T,
            #append = TRUE,
            sep = "\t")




ggplot(filtered_all_measurements,aes(x=condition,y=value))+
  geom_jitter()+
  facet_wrap(~measurement,scales="free_y")
  


anova <- filtered_all_measurements %>% 
  anova_test(value ~ condition)


stat.test <- filtered_all_measurements %>%
  tukey_hsd(value ~ condition) 


ggboxplot(filtered_all_measurements, x = "condition", y = "value",add = "jitter",facet.by = "measurement",scales="free_y")#+
  stat_pvalue_manual(
    stat.test, label = "p.adj", 
    y.position = c(70, 65, 60)
  )

ggsave("plot_with_stats.svg",width = 4,height = 4)



for_timecourse_csd <-filtered_data%>%
  pivot_wider(names_from = c(slice_id,measurement),values_from = value,names_sep = ".")%>%
  arrange(condition)%>%
  #select(-measurement)%>%
  pivot_longer(-condition,names_to = "slice_id")%>%
  cSplit("slice_id",".")%>%
  rename(slice_id = slice_id_1,
         measurement = slice_id_2)%>%
  filter(grepl("csd_source_vol_table|csd_total_postsynaptic_vol_table|csd_total_presynaptic_vol_table",measurement))%>%
  pivot_wider(names_from = measurement,values_from = value)


write.table(for_timecourse_csd, "csd_time_course_supplemental_data.dat",
            na = "",
            row.names = T,
            col.names = T,
            #append = TRUE,
            sep = "\t")











