library(tidyverse)
library(wesanderson)
theme_set(theme_bw()+ 
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                strip.background = element_blank(),
                panel.border = element_rect(size=0.5,colour="black"),
                axis.text = element_text(color="black")
                )
          )
