#!/usr/bin/env Rscript

##########################################################################################################################
#    Part of the FASAS project
#    File name: Tools-Len2Num_CumulativeCurve.r
#    Auther: Ke Zhang
#    Version: 1.0.0
#    Data: 2019.10.10
##########################################################################################################################

##########################################################################################################################
#    Script of plotting contig length cumulative curve
#    Copyright (C) 2019 CapitalBio Corporation
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
outputfile <- args[2]

library(ggplot2)
library(scales)

read.table(inputfile, header = T, sep = '\t', check.names = F) -> data
#cul max Y
int.y <- ceiling(max(data$Freq.) / 1000)
max.y <- int.y * 1000
if (max.y <= 20000){
    y.interval <- 1000
}else if (max.y <= 40000){
    y.interval <- 2000
}else if (max.y <= 60000){
    y.interval <- 3000
}else{
    y.interval <- 6000
}

pdf(outputfile, width = 10, height = 8)
ggplot() + geom_path(data = data, aes(x = Length, y = Freq., group = Sample, color = Sample)) +
    theme(axis.line = element_line(colour = 'black'),
        axis.text = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 60, hjust = 1),
        axis.ticks = element_line(colour = 'black'),
        axis.title = element_text(colour = 'black', face = 'bold', family = 'Helvetica'),
        plot.background = element_rect(colour = NA),
        panel.grid.major.y = element_line(linetype = 3, color="grey"),
        panel.grid.major.x = element_line(linetype = 3, color="grey"),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = 'white', colour = NA),
        legend.title = element_text(colour = 'black', face = 'bold', family = 'Helvetica'),
        legend.background = element_rect(colour = NA),
        legend.text=element_text(family = 'Helvetica', face = 'plain', size = 8),
    ) + scale_x_continuous(limits = c(300, 1800), breaks = seq(300, 1800, 100), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, max.y, y.interval), expand = c(0, 0))
dev.info <- dev.off()
