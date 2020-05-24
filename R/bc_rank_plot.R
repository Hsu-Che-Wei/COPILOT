#' Barcode Rank Plot
#' @param lnc_gg Datafrome that contains log (number of counts), rank and the whether the cell is selected for plotting
#' @param save Output directory for the plot
#' @import ggplot2

bc_rank_plot <- function(lnc_gg,save){
  png(save, width = 7, height = 5, units = 'in', res=300)
  print({
    ggplot(lnc_gg, aes(x=rank, y=lnc, col = select, alpha=select)) +
      geom_point(size=3) +
      scale_x_log10(labels = plain, breaks = trans_breaks("log10", function(x) round(10^x, 0))) +
      scale_y_log10(labels = plain, breaks = trans_breaks('log10', function(x) floor(10^x)))+
      scale_color_manual(values = c("#003B6D","#EBEDF3","#d96459","#f2e394"), name = NULL, guide = guide_legend(reverse = TRUE)) +
      scale_alpha_manual(values = c(0.8,0.8,0.8,0.8)) +
      labs(x = 'Barcodes', y = 'UMI Counts') +
      guides(alpha = FALSE, colour = guide_legend(reverse = TRUE, override.aes=list(size = 5))) +
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
            axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
            legend.background = element_rect(fill = 'transparent'),
            legend.position = "none",legend.key=element_blank(),
            plot.margin=unit(c(0,1,0,0),"cm"))
  })
  dev.off()
}
