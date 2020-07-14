#!/usr/bin/env Rscript

##########################################################################################################################
#    Part of the FASAS project
#    File name: step7-SamplingCoveragePlot.r
#    Auther: Ke Zhang
#    Version: 1.0.0
#    Data: 2019.10.10
##########################################################################################################################

##########################################################################################################################
#    Script of poltting sequence coverage
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
inputfile1 <- args[1] #R1
inputfile2 <- args[2] #R2
outputfile <- args[3] #output

#read csv
varstat.r1 <- read.csv( inputfile1, header = F)
varstat.r2 <- read.csv( inputfile2, header = F)

#varstat.ref <- read.csv( var.file.path, header = F)
var.region <- c(rep(0,70),rep(160,30),rep(0,37),rep(160,106),rep(0,190),rep(160,65),rep(0,78),rep(160,107),rep(0,139),rep(160,58),rep(0,106),rep(160,58),rep(0,73),rep(160,57),rep(0,69),rep(160,52),rep(0,140),rep(160,31),rep(0,135))
varstat.ref <- data.frame(V1=seq(0,1600),V2=var.region)


#get max y
varstat.r1.max <- max( varstat.r1[,2])
varstat.r2.max <- max( varstat.r2[,2])
max.y.value <- max( varstat.r1.max, varstat.r2.max) + 1

#varstat.ref 80:40000
varstat.ref.ratio <- max.y.value / 40000
varstat.ref$V2 <- floor( varstat.ref$V2 * varstat.ref.ratio)

pdf(outputfile,width = 8,height = 6)

ggplot() + geom_area(data = varstat.r1, aes(x = V1, y = V2), fill = NA, colour = 'orange') +
    geom_area(data = varstat.r2, aes(x = V1, y = V2), fill = NA, colour = 'mediumblue') +
    geom_area(data = varstat.ref, aes(x = V1, y = V2), fill = 'red', size = 1.1) +
    scale_y_continuous(limits = c(0, max.y.value), expand = c(0,0)) +
    scale_x_continuous(breaks = seq(0, 1600, 100), limits = c(0, 1600), expand = c(0, 0)) +
    xlab('V region') +
    ylab('Reads Num.') +
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
    )
dev.info <- dev.off()
