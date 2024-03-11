library(tidyverse)
library(tools)
library(rstatix)
library(svglite)


setwd("C:/Users/")

rr <- read_csv("stx_quantification.csv",col_names = F)%>%
  drop_na()%>%
  rename(., FileNameTbl = X1, Treatment = X2)%>%
  mutate(FileNameTbl = paste0(FileNameTbl,".abf"))


Folder_name<-"C:/Users/Damian/Dropbox/IGM projects/Caged_saxitoxin"

setwd(Folder_name)

file_list <- list.files(path=Folder_name,recursive = T, full.names = T, include.dirs = FALSE,pattern = "\\.txt$")

Raw_data <- file_list %>%
  set_names(nm = (basename(.) %>%
                    tools::file_path_sans_ext())) %>%
  map_df(read_tsv,
         col_names = T,
         .id = "File.name")


average_by_run <- Raw_data %>%
  group_by(FileNameTbl,File.name)%>%
  summarise(APs = mean(NumberAPsTbl))

joined <- inner_join(average_by_run,rr,by = "FileNameTbl")%>%
  mutate(Treatment = as_factor(Treatment))%>%
  mutate(Treatment = fct_relevel(Treatment, "Control", "UV","Saxitoxin","Saxitoxin + UV","Wash"))
 


ggplot(joined,aes(Treatment,APs,group = File.name))+
  geom_point()+
  geom_line()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  
  #scale_x_continuous(limits = c(0, 15:15))+
  #facet_wrap(~File.name,scales = "free")+
  ggsave("line_plot.png",width = 8,
         height = 8,
         units = "cm")

ggplot(joined,aes(Treatment,APs,label=File.name))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.1,alpha = 0.5,size = 3)+
  #geom_text(position=position_jitter(width=0.2))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x ="", y = "# APs")+
  #scale_y_continuous(limits = c(-1, 1))+
  
  ggsave("box_plot.svg",width = 10,
         height = 8,
         units = "cm")


###-----------------------------redrawn_with_stats

joined<-joined%>%
  mutate(FileNameTbl = as_factor(FileNameTbl))%>%
  mutate(File.name = as_factor(File.name))%>%
  mutate(Treatment = as_factor(Treatment))%>%
  mutate(APs = as.numeric(APs))

AP_table <- joined %>%
  pivot_wider(id_cols = -c(FileNameTbl),names_from = Treatment,values_from = APs )



write.table(
  AP_table,
  "Sweep_AP_count.txt",
  na = "",
  append = FALSE,
  sep = "\t",
  dec = ".",
  row.names = F,
  col.names = TRUE
)


stat_table = joined%>%
  ungroup()%>%
  filter(!(grepl("2021-06-08",File.name)))%>%
  select(-FileNameTbl)#%>%
  #mutate(APs = ifelse(APs==0,2,APs))#%>%
  #pivot_wider(names_from = Treatment,values_from = APs )
  
stat_table <- stat_table %>%
  mutate(File.name = droplevels(File.name))%>%
  mutate(Treatment = droplevels(Treatment))
  
  
stat.test <- stat_table %>%
  t_test(APs~Treatment, ref.group = "Control")%>%
  add_xy_position(x = "Treatment", fun = "max",scales = "free")

ggboxplot(stat_table, x = "Treatment", y = "APs",add = "jitter") +
  stat_pvalue_manual(
    stat.test, label = "{p.adj}"
  )+

ggsave("C:/Users/dw2471/Dropbox/TEMP/AP_firing.svg",width = 30,
         height = 30,
         units = "cm")

