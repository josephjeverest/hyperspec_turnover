# 00 - ggplot themes created for FuncDiv outputs
# Joseph Everest
# February 2021, adapted February 2022

# PACKAGES ----

library(ggplot2)
library(ggthemes)


# THEME MANUSCRIPT: For actual final figures

theme_ms <- function(){
  theme_bw() +
    theme(text = element_text(family = "Arial"),
          axis.text = element_text(size = 18), 
          axis.title = element_text(size = 18, face = "plain"),
          axis.line.x = element_line(color = "black"), 
          axis.line.y = element_line(color = "black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(0.1, 0.1, 0.75, 0.1), units = , "cm"),
          plot.title = element_text(size = 20, vjust = 1, hjust = -0.5, face = "plain"),
          plot.subtitle = element_text(size = 22, vjust = 1, hjust = 0.5, face = "bold"),
          plot.caption = element_text(size = 22, vjust = 0, hjust = 0.5, face = "bold"),
          legend.text = element_text(size = 18),          
          legend.title = element_text(size = 18, face = "bold", hjust = 0.5),
          legend.position = c(0.9, 0.9), 
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 2, linetype = "blank"))
}


# THEME 1: MODEL RUNS (based on theme by IMS) ----

theme_1 <- function(){
  theme_bw() +
    theme(text = element_text(family = "Helvetica Light"),
          axis.text = element_text(size = 16), 
          axis.title = element_text(size = 18),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          plot.title = element_text(size = 18, vjust = 1, hjust = 0, face = "bold"),
          plot.subtitle = element_text(size = 14, hjust = 0, face = "plain"),
          plot.caption = element_text(size = 12, hjust = 1, face = "plain"),
          legend.text = element_text(size = 12),          
          legend.title = element_text(size = 10, face = "bold", hjust = 0.5),
          legend.position = c(0.9, 0.9), 
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 2, linetype = "blank"))
}


# THEME 2: MODEL RUNS with same title and subtitle (based on theme by IMS) ----

theme_2 <- function(){
  theme_bw() +
    theme(text = element_text(family = "Helvetica Light"),
          axis.text = element_text(size = 16), 
          axis.title = element_text(size = 18),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          plot.title = element_text(size = 14, vjust = 1, hjust = 0, face = "plain"),
          plot.subtitle = element_text(size = 12, hjust = 0, face = "plain"),
          plot.caption = element_text(size = 12, hjust = 1, face = "plain"),
          legend.text = element_text(size = 12),          
          legend.title = element_text(size = 10, face = "bold", hjust = 0.5),
          legend.position = c(0.9, 0.9), 
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 2, linetype = "blank"))
}


# THEME 3: MAPS with same title and subtitle (based on theme by GD) ----

theme_3_map <- function(){
  theme_map() +
    theme(text = element_text(family = "Helvetica Light"),
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 14),
          axis.line.x = element_line(color = "black"), 
          axis.line.y = element_line(color = "black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          plot.title = element_text(size = 16, vjust = 1, hjust = 0, face = "bold"),
          plot.subtitle = element_text(size = 14, hjust = 0, face = "plain"),
          plot.caption = element_text(size = 12, hjust = 1, face = "plain"),
          legend.text = element_text(size = 12),          
          legend.title = element_text(size = 10, face = "bold", hjust = 0.5),
          legend.position = "right",
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 2, linetype = "blank"),
          plot.background = element_rect(fill = "white", color = NA))
}


# THEME 4: MAPS produced using ggOceanMaps ----

theme_4_map <- function(){
  theme(text = element_text(family = "Helvetica Light"),
        plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
        plot.title = element_text(size = 16, vjust = 1, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 14, hjust = 0.5, face = "plain"),
        plot.caption = element_text(size = 12, hjust = 1, face = "plain"),
        legend.text = element_text(size = 12),          
        legend.title = element_text(size = 10, face = "bold", hjust = 0.5),
        legend.position = "right", 
        legend.key = element_blank(),
        legend.background = element_rect(color = "black", 
                                         fill = "transparent", 
                                         size = 2, linetype = "blank"),
        plot.background = element_rect(fill = "white", color = NA))
}
