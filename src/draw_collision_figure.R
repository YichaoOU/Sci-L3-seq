
library(ggplot2)
library(readr)
library(dplyr) 
args <- commandArgs(trailingOnly=TRUE)
inputFile = args[1]
outputFile = args[2]

plot_collision <- function(df, threshold=0.95) {
  xmax <- max(df$HUMAN_READS * 1.1)
  ymax <- max(df$MOUSE_READS * 1.1)
  df %>%
  mutate(CELL=case_when(HUMAN_PERCENT >= threshold ~ "Human",
                        MOUSE_PERCENT >= threshold ~ "Mouse",
                        TRUE ~ "Collision")) %>%
    mutate(CELL=factor(CELL, levels=c("Human", "Collision", "Mouse"))) %>%
    ggplot(aes(HUMAN_READS, MOUSE_READS, color=CELL)) +
    geom_point(alpha=0.3) +
    xlab("# of human unique insertions") +
    ylab("# of mouse unique insertions") +
    scale_color_manual(values=c("Human"="orangered", "Mouse"="dodgerblue", "Collision"="grey80")) +
    scale_x_continuous(labels=scales::format_format(scientific = FALSE),
                       breaks=scales::pretty_breaks(n=4)) +
    scale_y_continuous(labels=scales::format_format(scientific = FALSE),
                       breaks=scales::pretty_breaks(n=4)) +
    coord_cartesian(xlim=c(0, xmax), ylim=c(0, ymax)) +
    theme_bw() +
    theme(legend.title=element_blank(),
          legend.justification=c(1, 1),
          legend.position=c(0.95, 0.95),
          legend.background=element_rect(colour="black", size=0.5),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
}

df <- read_delim(inputFile, col_names=c("CELL_ID", "HUMAN_READS", "MOUSE_READS", "TOTAL_READS", "HUMAN_PERCENT", "MOUSE_PERCENT"), delim=" ")


p <- plot_collision(df)
ggsave(outputFile)