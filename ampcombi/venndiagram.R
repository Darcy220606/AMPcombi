#!/usr/bin/env Rscript

##############################
# Rscript to generate VennDiagrams on AMP-containing contigs by detection tool ####
##############################
# Date ####
# October, 19 2022
##############################
# Authors ####
# Anan Ibrahim - ananhamido@hotmail.com - @darcy220606
# Louisa Perelo - louperelo@gmail.com - @louperelo
##############################
# Working_directory ####
#setwd(getwd())
setwd('/Users/lp1/Documents/nextflow/funcscan/ampcombi/')
##############################
# Libraries used + arguments ####
if (!require("dplyr")) install.packages('dplyr')
if (!require("DT")) install.packages('DT')
if (!require("optparse")) install.packages('optparse')
if (!require("htmlwidgets")) install.packages('htmlwidgets')
if (!require("VennDiagram")) install.packages('VennDiagram')
if (!require("eulerr")) install.packages('eulerr')

#if (!require("limma")) install.packages('limma')

library("dplyr")
library("DT")
library("optparse")
library("htmlwidgets")
library("VennDiagram")
library("eulerr")
#library("limma")

#option_list = list(
#  make_option(c("-f", "--file"), type="character", default="AMPcombi_summary.csv",
#              help="AMpcombi complete summary table [default= %default]", metavar="character"),
#  make_option(c("-o", "--out"), type="character", default="VennDiagram",
#              help="Provide the name of the output file [default= %default]", metavar="character"));
# Turns warnings off
#options(warn=-1)
#opt_parser = OptionParser(option_list=option_list);
#opt = parse_args(opt_parser);

##############################
#Generate HTML interactive files ####
table <-
  #readr::read_csv(opt$file,show_col_types = FALSE) %>%
  readr::read_csv("AMPcombi_summary.csv",show_col_types = FALSE) %>%
  unique()

# matching columns should contain 'prob_' or 'evalue_hmmer'
#TODO: select columns to be included (toMatch) by clicking on the column headers
# in the interactive table
toMatch <- c('prob_', 'evalue_hmmer')
#TODO: include selection of sample name, otherwise use all samples in counts

# identify matching columns indices
tool_cols <- grep(paste(toMatch,collapse="|"), colnames(table))
# new df containing data for VennDiagram (sample, contig_id, toMatch columns)
venn_df <- table[,c(1,2,tool_cols)]
# replace values with 1 == AMP identified by tool on contig
venn_df[,tool_cols] <- replace(venn_df[,tool_cols], venn_df[,tool_cols]>0, 1)

#Determine how many variables will enter the VennDiagram (max 5!)
venn_num <- sum(grepl(paste(toMatch,collapse="|"), colnames(table)))
#venn_num <- 2

# create VENN diagram
if (venn_num >=1)
{
  A = venn_df[,tool_cols][1]
}
if (venn_num >1)
{
  B = venn_df[,tool_cols][2]
}
if (venn_num >=2)
{
  C = venn_df[,tool_cols][3]
}
if (venn_num >=3)
{
  D = venn_df[,tool_cols][4]
}
if (venn_num ==5)
{
  E = venn_df[,tool_cols][5]
}

grid.newpage()
# single Venn
if (venn_num ==1){
  draw.single.venn(sum(A), 
                   colnames(venn_df[0,tool_cols[1:venn_num]]),
                   fill = c("orange"));
}
# pairwise Venn
if (venn_num ==2){
  draw.pairwise.venn(area1 = sum(A),
                     area2 = sum(B),
                     cross.area = sum(A == 1 & B == 1),
                     category = c(colnames(venn_df[0, tool_cols[1:venn_num]])),
                     fill = c("orange", "red"),
                     cat.col = c("orange", "red")
                    )
}
# triple Venn
if (venn_num ==3){
  draw.triple.venn(area1 = sum(A),
                      area2 = sum(B),
                      area3 = sum(C),
                      n12 = sum(A == 1 & B == 1),
                      n13 = sum(A == 1 & C == 1),
                      n23 = sum(B == 1 & C == 1),
                      n123 = sum(A == 1 & B == 1 & C == 1),
                      category = c(colnames(venn_df[0, tool_cols[1:venn_num]])),
                      fill = c("orange", "red", "green"),
                      cat.col = c("orange", "red", "green")
  )
}
# quadruple Venn
if (venn_num ==4){
  draw.quad.venn(area1 = sum(A),
                      area2 = sum(B),
                      area3 = sum(C),
                      area4 = sum(D),
                      n12 = sum(A == 1 & B == 1),
                      n13 = sum(A == 1 & C == 1),
                      n14 = sum(A == 1 & D == 1),
                      n23 = sum(B == 1 & C == 1),
                      n24 = sum(B == 1 & D == 1),
                      n34 = sum(C == 1 & D == 1),
                      n123 = sum(A == 1 & B == 1 & C == 1),
                      n124 = sum(A == 1 & B == 1 & D == 1),
                      n134 = sum(A == 1 & C == 1 & D == 1),
                      n234 = sum(B == 1 & C == 1 & D == 1),
                      n1234 = sum(A == 1 & B == 1 & C == 1 & D == 1),
                      category = c(colnames(venn_df[0, tool_cols[1:venn_num]])),
                      fill = c("orange", "red", "green", "blue"),
                      cat.col = c("orange", "red", "green", "blue")
  )
}
# quintuple Venn
if (venn_num ==5){
  draw.quintuple.venn(area1 = sum(A),
                      area2 = sum(B),
                      area3 = sum(C),
                      area4 = sum(D),
                      area5 = sum(E),
                      n12 = sum(A == 1 & B == 1),
                      n13 = sum(A == 1 & C == 1),
                      n14 = sum(A == 1 & D == 1),
                      n15 = sum(A == 1 & E == 1),
                      n23 = sum(B == 1 & C == 1),
                      n24 = sum(B == 1 & D == 1),
                      n25 = sum(B == 1 & E == 1),
                      n34 = sum(C == 1 & D == 1),
                      n35 = sum(C == 1 & E == 1),
                      n45 = sum(D == 1 & E == 1),
                      n123 = sum(A == 1 & B == 1 & C == 1),
                      n124 = sum(A == 1 & B == 1 & D == 1),
                      n125 = sum(A == 1 & B == 1 & E == 1),
                      n134 = sum(A == 1 & C == 1 & D == 1),
                      n135 = sum(A == 1 & C == 1 & E == 1),
                      n145 = sum(A == 1 & D == 1 & E == 1),
                      n234 = sum(B == 1 & C == 1 & D == 1),
                      n235 = sum(B == 1 & C == 1 & E == 1),
                      n245 = sum(B == 1 & D == 1 & E == 1),
                      n345 = sum(C == 1 & D == 1 & E == 1),
                      n1234 = sum(A == 1 & B == 1 & C == 1 & D == 1),
                      n1235 = sum(A == 1 & B == 1 & C == 1 & E == 1),
                      n1245 = sum(A == 1 & B == 1 & D == 1 & E == 1),
                      n1345 = sum(A == 1 & C == 1 & D == 1 & E == 1),
                      n2345 = sum(B == 1 & C == 1 & D == 1 & E == 1),
                      n12345 = sum(A == 1 & B == 1 & C == 1 & D == 1 & E == 1),
                      category = c(colnames(venn_df[0, tool_cols[1:venn_num]])),
                      fill = c("orange", "red", "green", "blue", "grey"),
                      cat.col = c("orange", "red", "green", "blue", "grey")
  )
}
# if more than 5 tools are to be compared, try to draw an EULER diagram
# the euler diagram may not be able to show all overlaps: " if the euler model can not describe all intersections the count will not be ploted."
# https://stackoverflow.com/questions/46629193/create-5-way-venn-diagram-from-csv-file-in-r
if (venn_num >5){
  set.seed(10) #this seed changes the orientation of the sets         
  plot(euler(venn_df[,tool_cols]), counts = T, fontface = 1, quantities = TRUE)
}


             
