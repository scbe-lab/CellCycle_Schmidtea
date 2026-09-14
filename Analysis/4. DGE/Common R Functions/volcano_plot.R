require(ggplot2)
require(ggrastr)
volcano_plot <- function(d,fc_thres = 1.5, p_thres = 0.1, main = "Volcano plot", ylim = NULL){
  
  if(is.null(ylim)) ylim = max(-log(d$padj))
  
  p <-
    d %>% arrange(deg) %>%
    ggplot(aes(x = log2FoldChange, y = -log(padj), fill = deg, colour = deg))+
    geom_point_rast(shape = 21)+
    scale_fill_manual(values = setNames(c("#FF6347","#24242430"),c("DEG","none")))+
    scale_colour_manual(values = colorspace::darken(setNames(c("tomato","#24242460"),c("DEG","none")),.3))+
    theme_classic()+
    geom_vline(xintercept=c(-fc_thres), color="#a6c8e2",linetype = "longdash")+
    geom_vline(xintercept=c(fc_thres), color="#a6c8e2",linetype = "longdash")+
    geom_hline(yintercept=-log(p_thres), color="#a6c8e2",linetype = "longdash")+
    ylim(0,ylim)+
    ggtitle(main)
  
  p
}