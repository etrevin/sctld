require(ggplot2)
library(scales)

#resultKS_bp_ac_low_high_down             <- runTest(GOdata_bp_ac_low_high_down,             algorithm = "parentchild", statistic = "fisher")  #cc
#allRes_cc_ac_low_high_down   <- GenTable(GOdata_bp_ac_low_high_down,  ParentChild_Fisher = resultKS_bp_ac_low_high_down,   topNodes = 10, orderBy = ParentChild_Fisher, numChar = 99)
allRes_mf_vt_up2 = allRes_mf_ac_low_high_down
#allRes_mf_vt_up2[2,]=NaN
#allRes_mf_vt_up2[26,]=NaN
allRes_mf_vt_up2$ParentChild_Fisher<- as.numeric(allRes_mf_vt_up2$ParentChild_Fisher)
allRes_mf_vt_up2 <- allRes_mf_vt_up2[allRes_mf_vt_up2$ParentChild_Fisher < 0.05,]
allRes_mf_vt_up2 <- allRes_mf_vt_up2[, c("GO.ID","Term","ParentChild_Fisher","Annotated","Significant")]
#ggdata <- allRes_mf_vt_up2[1:10,] 
ggdata2 <- allRes_mf_vt_up2[1:2,]
ggdata2$Annotated <- c(1,100) 
ggdata2$Significant <- c(1,20) 
ggdata2$Term <- c("ano","sig") 
#ggdata = rbind(allRes_mf_vt_up2,ggdata2)
ggdata = allRes_mf_vt_up2
ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) # fixes order
gg1 <- ggplot(ggdata,
              aes(x = Term, y = ParentChild_Fisher, size = Significant, fill = Annotated)) +
  
  expand_limits(y = 0.05) +
  geom_point(shape = 21) +
  scale_size(range = c(2,12.5)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4')  +
    xlab('') + ylab('Enrichment score') +
  labs(
    title = 'GO Molecular Function ',
    subtitle = 'Top 5 terms ordered by Kolmogorov-Smirnov p-value',
    caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001') +
  
  geom_hline(yintercept = c(0.05, 0.01, 0.001),
             linetype = c("dotted", "longdash", "solid"),
             colour = c("black", "black", "black"),
             size = c(0.5, 1.5, 3)) +
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 8, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 8, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 8, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 8, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 8, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 8, face = 'bold'),
    axis.title.x = element_text(size = 8, face = 'bold'),
    axis.title.y = element_text(size = 8, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 8, face = "bold"), # Text size
    title = element_text(size = 8, face = "bold")) +
  
  coord_flip()
gg1

