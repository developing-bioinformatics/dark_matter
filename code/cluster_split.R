#!/usr/bin/R
library(rBLAST)
library(Biostrings)
library(ggplot2)
library(ShortRead)
library(dplyr)
library(stringr)
library(DECIPHER)
library(taxonomizr)
library(dendextend)
library(parallel)
options(scipen=999, digits=1) 

function(clusterMax, subtrees=50){

  dend_list <- get_subdendrograms(clusterMax, subtrees)
  return(dend_list)
}