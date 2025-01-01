library(tidyverse)
library(patchwork)
theme_set(theme_bw(base_size = 10)) 


patchwork <- (gamplt + mmplt)/ (gamdoyplt +  lgvarplt)

patchwork + plot_annotation(tag_levels = 'A', tag_suffix = ')')
  
  
(fig5b + fig5c) + plot_annotation(tag_levels = list(c('B', 'C')), tag_suffix = ')')

  
  
  
  
ggsave(filename = "fig5BC.png", path = "C:/Users/tmwepprich/Desktop/voltinism_revision", 
      device='png', dpi=600)

  
  
  