
setwd("C:/Users/Volpato/OneDrive/Sorgo MTME/Analises_prior/test_corr")
corr<-read.csv("corr_VCOV_Bayes.csv")
corr

library(tidyverse)
library(corrr)
library(igraph)
library(ggraph)

set.seed(778)
graph_cors <- corr %>%
  graph_from_data_frame(directed = FALSE)

graph_cors

ggraph(graph_cors) +
  geom_edge_link() +
  geom_node_point() +
  geom_node_text(aes(label = name))

set.seed(778)
ggraph(graph_cors) +
  geom_edge_link(aes(edge_alpha = abs(r), edge_width = abs(r), color = r)) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_colour_gradientn(limits = c(-1, 1),  colors = c("firebrick4","firebrick2", "gray", "dodgerblue2", "dodgerblue4")) +
  geom_node_point(color = "seashell2", size = 5) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_graph() 


#############FREQUENTISTA###########

setwd("C:/Users/Volpato/OneDrive/Sorgo MTME/Analises_prior/test_corr")
corr<-read.csv("corr_VCOV_Freq.csv")
corr

library(tidyverse)
library(corrr)
library(igraph)
library(ggraph)

set.seed(799)
graph_cors <- corr %>%
  graph_from_data_frame(directed = FALSE)

graph_cors

ggraph(graph_cors) +
  geom_edge_link() +
  geom_node_point() +
  geom_node_text(aes(label = name))

set.seed(799)
ggraph(graph_cors) +
  geom_edge_link(aes(edge_alpha = abs(r), edge_width = abs(r), color = r)) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_colour_gradientn(limits = c(-1, 1),  colors = c("firebrick4","firebrick2", "gray", "dodgerblue2", "dodgerblue4")) +
  geom_node_point(color = "seashell2", size = 5) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_graph() 


