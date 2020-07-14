#!/usr/bin/env Rscript

##########################################################################################################################
#    Part of the FASAS project
#    File name: step4-UMIFrequencyDistrabutionPlot.r
#    Auther: Ke Zhang
#    Version: 1.0.0
#    Data: 2019.10.10
##########################################################################################################################

##########################################################################################################################
#    Script of poltting UMI Frequency
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

#usage: Rscript step4-UMIFrequencyDistrabutionPlot.r StatUMI_1to2.log StatUMI_2to1.log StatUMI_freq.pdf
library(ggplot2)

args <- commandArgs(T)
inputfile1 <- args[1]
inputfile2 <- args[2]
outputfile <- args[3]

#read data
read.table(inputfile1, sep = "\t", fill = T, header = F) -> U1.col
read.table(inputfile2, sep = "\t", fill = T, header = F) -> U2.col

#change data
colnames(U1.col) <- c("U1","U1.Num")
colnames(U2.col) <- c("U2","U2.Num")

U1.Freq <- data.frame(table(U1.col$U1.Num))
U2.Freq <- data.frame(table(U2.col$U2.Num))
U1.Freq$Var1 <- as.numeric(as.character(U1.Freq$Var1))
U2.Freq$Var1 <- as.numeric(as.character(U2.Freq$Var1))

#calculate range
U1.Freq.Var1.Max <- max(U1.Freq$Var1)
U2.Freq.Var1.Max <- max(U2.Freq$Var1)
if (U1.Freq.Var1.Max > U2.Freq.Var1.Max){
    Max.Var1 <- U1.Freq.Var1.Max
}else{
    Max.Var1 <- U2.Freq.Var1.Max
}

if (Max.Var1 > 1000){
    Spac.x <- as.integer(Max.Var1/50)
}
if (Max.Var1 > 200){
    Spac.x <- as.integer(Max.Var1/20)
}else{
    Spac.x <- 10
}

U1.Freq.Freq.Max <- max(U1.Freq$Freq)
U2.Freq.Freq.Max <- max(U2.Freq$Freq)
if (U1.Freq.Freq.Max > U2.Freq.Freq.Max){
    Max.Freq <- U1.Freq.Freq.Max
}else{
    Max.Freq <- U2.Freq.Freq.Max
}

#get anno coordinate
Anno.x <- as.integer(Max.Var1 * 0.75)
Anno.y <- as.integer(Max.Freq * 0.75)

#plot by ggplot
p <- ggplot() + geom_area(data = U1.Freq, aes(x=Var1,y=Freq), fill = NA, colour = 'green') +
    geom_area(data = U2.Freq,aes(x=Var1, y=Freq), fill = NA, colour = 'red') +
    xlab('The number of UMI(>2)') +
    ylab('Freq.') +
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
    ) + scale_x_continuous(breaks = seq(0, Max.Var1, Spac.x), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))

p <- p + annotate("text", x = Anno.x, y = Anno.y, label = "U1 is green, U2 is red", size = 4)

#output
pdf(outputfile)
p
dev.info <- dev.off()
