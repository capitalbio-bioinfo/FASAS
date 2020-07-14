#!/usr/bin/env Rscript

##########################################################################################################################
#    Part of the FASAS project
#    File name: step4-UMICumulativeCurve.r
#    Auther: Ke Zhang
#    Version: 1.0.0
#    Data: 2019.10.10
##########################################################################################################################

##########################################################################################################################
#    Script of plotting cumulative curve
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

#usage: Rscript step4-UMIFrequencyDistrabutionPlot.r StatUMI_LtoR.tsv StatUMI_RtoL.tsv StatUMI_freq.pdf
library(ggplot2)

args <- commandArgs(T)
inputfile1 <- args[1]
inputfile2 <- args[2]
outputfile <- args[3]

#read data
read.table(inputfile1, sep = "\t", fill = T, header = F) -> data
data[1:2] -> U1.col
read.table(inputfile2, sep = "\t", fill = T, header = F) -> data2
data2[1:2] -> U2.col

#build data
k <- rep(0,length(U1.col[,2]))
for (i in 1:length(U1.col[,2])){
    if (i == 1){
        k[i] <- U1.col[i,2]
    }else{
        k[i] <- k[i-1] + U1.col[i,2]
    }
}
l <- rep(0,length(U2.col[,2]))
for (i in 1:length(U2.col[,2])){
    if (i == 1){
        l[i] <- U2.col[i,2]
    }else{
        l[i] <- l[i-1] + U2.col[i,2]
    }
}

#calculate range
if (k[length(k)] > l[length(l)]){
    Max.y.value <- k[length(k)]
}else{
    Max.y.value <- l[length(l)]
}
Step.y.value <- as.integer(Max.y.value / 10)
if (length(k) > length(l)){
    Max.x.value <- length(k)
}else{
    Max.x.value <- length(l)
}
Step.x.value <- as.integer(Max.x.value / 10)


#convert data.frame
data.frame(k,rep('Left',length(k)),1:length(k)) -> k.data
data.frame(l,rep('Right',length(l)),1:length(l)) -> l.data
colnames(k.data) <- c('Value','Group','Index')
colnames(l.data) <- c('Value','Group','Index')
rbind(k.data,l.data) -> all.data

#ggplot2 plot
pdf(outputfile)

ggplot(data = all.data)+geom_line(aes(Index,Value,colour = Group)) +
    xlab('Index') +
    ylab('UMI Number') +
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
        legend.text=element_text(family = 'Helvetica',face = 'bold'),
    ) + scale_x_continuous(breaks = seq(0, Max.x.value, Step.x.value), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, Max.y.value, Step.y.value))

dev.info <- dev.off()
