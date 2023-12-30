library(tidyverse)
library(gplots)
library("RColorBrewer")
library(ggpubr)
library(ggthemes)

setwd("/Users/dorota/Dropbox/UC_DAVIS/PROMISE/Methionine_experiment/DEG_v2/")
#background
background <- read_csv("Methionine_cpmExpression.csv")
universe <- nrow(background)
universe

#DE lists
Met <- read.csv("DiffExprs/GeneLists/Met_effect.csv", header = T,stringsAsFactors = F)
Met_up <- subset(Met, logFC>=1)
Met_down <- subset(Met, logFC<=-1)

#Met_up
#enriched categories
Metup.GO <- read_csv("GO enrichments/GO_FC1_Met_up.csv")
#DEGs
Met_upDEGs <- Met_up 

numAllDE <- nrow(Met_up)
numAllDE

Metup <- Metup.GO  %>% 
  mutate("DEG proportion"=numDEInCat/numInCat) %>%
  mutate (DEGsPerc=numDEInCat/numAllDE) %>% mutate (GOsPerc=numInCat/universe) %>% mutate (odds.ratio=DEGsPerc/GOsPerc) %>%
  filter(odds.ratio>1) %>%
  filter(numInCat>1)
write_csv(Metup, "GO enrichments/Met up GOs FC1-processed.csv")

#Met_down
#enriched categories
Metdown.GO <- read_csv("GO enrichments/GO_FC1_Met_down.csv")
#DEGs
Met_downDEGs <- Met_down

numAllDE <- nrow(Met_down)
numAllDE

Metdown <- Metdown.GO  %>% 
  mutate("DEG proportion"=numDEInCat/numInCat) %>%
  mutate (DEGsPerc=numDEInCat/numAllDE) %>% mutate (GOsPerc=numInCat/universe) %>% mutate (odds.ratio=DEGsPerc/GOsPerc) %>%
  filter(odds.ratio>1)%>%
  filter(numInCat>1)
write_csv(Metdown, "GO enrichments/Met down GOs FC1-processed.csv")


###Visualise it!


pal <- wes_palette("Zissou1", n = 4)

Metdown_vis <- ggplot(Metdown, aes(x = reorder(term, -odds.ratio), y = odds.ratio, fill = over_represented_pvalue)) + 
  geom_bar(stat = "identity") + 
  coord_flip()+
  theme_classic()+
  theme(legend.position="bottom")+
  scale_fill_gradientn(colours = pal, trans = 'reverse')+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),  
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 12)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 70))+
  ggtitle("GO FC1 categories enriched within methionie down-regulated genes")

Metup_vis <- ggplot(Metup, aes(x = reorder(term, -odds.ratio), y = odds.ratio, fill = over_represented_pvalue)) + 
  geom_bar(stat = "identity") + 
  coord_flip()+
  theme_classic()+
  theme(legend.position="bottom")+
  scale_fill_gradientn(colours = pal, trans = 'reverse')+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),  
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 12)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 70))+
  ggtitle("GO FC1 categories enriched within methionie up-regulated genes")


