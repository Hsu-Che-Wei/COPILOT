#' Gene histogram
#' @param lnc_gg Datafrome that contains log (number of genes), rank and the whether the cell is selected for plotting
#' @param save Output directory for the plot
#' @import ggplot2
#' @import scales

Gene_hist_plot <- function(lng_gg,save){

  library(ggplot2)
  library(scales)

  png(save, width = 7, height = 5, units = 'in', res=300)
  print({
    ggplot(lng_gg, aes(lng, fill = select)) +
      geom_histogram(bins = 100)+
      scale_fill_manual(values = c("#003B6D","#EBEDF3","#d96459","#f2e394"))+
      xlab("Number of Genes")+
      ylab("Freqency")+
      scale_y_continuous(labels = plain)+
      scale_x_log10(labels = plain, breaks = trans_breaks("log10", function(x) round(10^x, 0))) +
      theme(panel.background=element_blank(),plot.background=element_blank(),legend.title = element_blank(), axis.line = element_line(colour = "black"),
            axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
            axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),legend.position="none",
            plot.margin=unit(c(0,1,0,0),"cm"))
  })
  dev.off()
}
