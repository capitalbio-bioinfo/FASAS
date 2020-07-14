#!/usr/bin/env Rscript

##########################################################################################################################
#    Part of the FASAS project
#    File name: step7-SamplingCoveragePlot.r
#    Auther: Ke Zhang
#    Version: 1.0.0
#    Data: 2019.10.10
##########################################################################################################################

##########################################################################################################################
#    Script of plotting Sequence Coverage
#    Copyright (C) 2018-2020  CapitalBio Corporation
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
##########################################################################################################################

args <- commandArgs(T)
inputfile <- args[1]
cutoff <- args[2]
outputfile <- args[3]

library(ggplot2)

read.table(inputfile, header = F) -> data

as.data.frame(table(data)) -> data2
data.MAX.y <- max(data2$Freq)
Anno.y <- as.integer(data.MAX.y * 0.75)
cutoff <- as.numeric(cutoff)
Anno.x <- cutoff + 15
Anno.label <- paste0("f(x)=", cutoff)

p <- ggplot() + geom_histogram(data = data, aes(x = V1), binwidth = 0.2, fill = "white", colour = "red") +
    geom_vline(xintercept = cutoff, linetype = "dotted") +
    xlab('The length of Sequence') + ylab('Freq.') +
    theme(axis.line = element_line(colour = 'black'),
        axis.text = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 60, hjust = 1),
        axis.ticks = element_line(colour = 'black'),
        axis.title = element_text(colour = 'black', family = 'Helvetica', face = 'bold'),
        plot.background = element_rect(colour = NA),
        panel.grid.minor =element_line(colour = NA),
        panel.background =element_rect(fill = 'white', colour = NA),
        legend.title=element_blank(),
        legend.background = element_rect(colour = NA),
        legend.text=element_text(family = 'Helvetica', face = 'bold'),
    ) + scale_x_continuous(breaks = seq(50,300,10), expand = c(0,0)) + scale_y_continuous(expand = c(0, 0))
p <- p + annotate("text", x = Anno.x, y = Anno.y, label = Anno.label, size = 4)
pdf(outputfile)
p
dev.info <- dev.off()
