library(tidyverse)
library(ggalluvial)
library(viridis)
dat_lmtp <- read_rds(here::here("data/derived/dat_final.rds")) %>%
  mutate(I_00 = case_when(hypoxia_ed == 1 & I_00 == 0 ~ 1, TRUE ~ I_00))

dat_lmtp %>%
  select(id, event, cr, starts_with("I_")) %>%
  pivot_longer(cols = starts_with("I_")) %>%
  mutate(day = parse_number(name),
    status = case_when(value == 0 ~ "No Supp O2",
                       value == 1 ~ "Supp O2",
                       value == 2 ~ "Intubated",
                       event == 1 ~ "AKI",
                       cr == 1 ~ "Deceased",
                       TRUE ~ "Discharged"
                       )) %>%
  group_by(day, status)

outcome_day <- 15
padded_days <- str_pad(0:(outcome_day-1), 2, pad = "0")
library(ggsci)
vars_to_group <- paste0("I_",padded_days)
type_colors <- c("#FF0000","#00A08A","#F98400","navy","#5BBCD6","gray")
type_colors <- c("black","#FF0000","dodgerblue3","#00A08A","navy")
type_colors <- c("#CC2936","#0B0033","#157F1F","#3F88C5","#70EE9C")

library(extrafont)
font_import()
loadfonts(device="win") 

alluv <- dat_lmtp %>%
  select(id, event, cr, starts_with("I_")) %>%
  mutate(across(starts_with("I_"),
                            ~ case_when(.x == 0 ~ "No supplemental oxygen",
                            .x == 1 ~ "Non-IMV supplemental oxygen",
                            .x == 2 ~ "IMV",
                            event == 1 ~ "AKI",
                            cr == 1 ~ "Deceased",
                            TRUE ~ "Discharged"
         ))) %>%
  group_by_at(vars_to_group) %>%
  count() %>%
  to_lodes_form(key = "day", axes = 1:15) %>%
  mutate(day = factor(parse_number(as.character(day))),
        stratum = fct_rev(fct_relevel(stratum,
                                      "Discharged", "No supplemental oxygen", "Non-IMV supplemental oxygen", "IMV", "AKI","Deceased"))) %>%
  # filter(stratum != "Discharged") %>%
  ggplot(
    aes(x = day,
        stratum = stratum, 
        col = stratum,
        alluvium = alluvium,
        y=n,
        fill = stratum)
  ) +
  #geom_flow(col="white") +
  #geom_stratum(col="white") +
  geom_flow() +
  geom_stratum() +
  #scale_color_viridis(discrete=TRUE) +
  #scale_fill_viridis(discrete=TRUE) +
  #scale_fill_brewer(type = "qual", palette = "Set1") +
  #scale_color_brewer(type = "qual", palette = "Set1") +
  #scale_fill_manual(values = type_colors) +
  #scale_color_manual(values = type_colors) +
  scale_fill_jama() +
  scale_color_jama( ) +
  theme_classic() +
  scale_x_discrete(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Study Day", y = "Number of Patients", title = "Number of Patients per Exposure-Outcome Status by Day",
       fill = "Status", col = "Status") +
  theme(text=element_text(family="Times", size=11),
        #legend.spacing.y = unit(1, 'cm'),
        #legend.position = "bottom",
        legend.text = element_text(#margin = margin(t = 5, b=5),
          size=9))
        #legend.key.size = unit(1, "cm"))
        
        #legend.position = c(.8,.8),
        #legend.box.background = element_rect(color="gray", size=.7))
        #plot.title=element_text( size=14, face="bold"),
        #axis.title.x =element_text(face="bold"),
        #legend.title =element_text(face="bold"),
        #axis.title.y =element_text(face="bold")
alluv
ggsave( "output/figure_alluvial.pdf", alluv, width=8.5, height=4)

#+
  #theme(
  ##  panel.background = element_rect(fill = "black",
   #                                 colour = "black"))
  #geom_text(stat = "stratum", size = 3)
  