#!/usr/bin/env Rscript

##########################################################################################################################
#    Part of the FASAS project
#    File name: step8-ContigLengthDistrabutionPlot.r
#    Auther: Ke Zhang
#    Version: 1.0.0
#    Data: 2019.10.10
##########################################################################################################################

##########################################################################################################################
#    Script of poltting contig length distrabution
#    Copyright (C) 2019  CapitalBio Corporation
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

library(ggplot2)

args <- commandArgs(T)
inputfile <- args[1]
outputfile <- args[2]

read.table(inputfile, sep='\t') -> data
colnames(data) <- c("Length","Freq.")

pdf(outputfile)
ggplot()+geom_area(data = data,aes(x=Length,y=Freq.), fill = NA, colour = 'yellowgreen') +
    xlab('The Length of the Contig') +
    ylab('Num.') +
    theme(axis.line = element_line(colour = 'black'),
        axis.text = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 60, hjust = 1),
        axis.ticks = element_line(colour = 'black'),
        axis.title = element_text(colour = 'black', family = 'Helvetica', face = 'bold'),
        plot.background = element_rect(colour = NA),
        panel.grid.major.y = element_line(linetype = 3, color="grey"),
        panel.grid.major.x = element_line(linetype = 3, color="grey"),
        panel.grid.minor =element_line(colour = NA),
        panel.background =element_rect(fill = 'white', colour = NA),
    ) + scale_x_continuous(breaks = seq(300, 1800, 100), limits = c(300, 1800), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))
dev.info <- dev.off()
