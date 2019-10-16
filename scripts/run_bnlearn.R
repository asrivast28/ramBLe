#!/usr/bin/env Rscript
##
# @file run_bnlearn.R
# @brief Script for running bnlearn with the same options as csl. 
# @author Ankit Srivastava <asrivast@gatech.edu>
#/

library('optparse')
parser <- OptionParser()
parser <- add_option(parser, c('--nvars', '-n'), type='integer', help='Number of variables in the dataset.')
parser <- add_option(parser, c('--nobs', '-m'), type='integer', help='Number of observations in the dataset.')
parser <- add_option(parser, c('--file', '-f'), type='character', help='Name of the file from which dataset is to be read.')
parser <- add_option(parser, c('--colobs', '-c'), action='store_true', default=FALSE, help='The file contains observations in columns.')
parser <- add_option(parser, c('--separator', '-s'), type='character', default=',', help='Delimiting character in the file.')
parser <- add_option(parser, c('--varnames', '-v'), action='store_true', default=FALSE, help='Read variable names from the first row of the file.')
parser <- add_option(parser, c('--algorithm', '-a'), type='character', default='gs', help='Name of the algorithm to be used.')
parser <- add_option(parser, c('--target', '-t'), type='character', help='Name of the target variable.')
parser <- add_option(parser, c('--blanket', '-b'), action='store_true', help='Find MB instead of PC for the target var.')
parser <- add_option(parser, c('--output', '-o'), type='character', help='Name of the file to which the learned network should be written.')
parser <- add_option(parser, c('--directed', '-d'), action='store_true', default=FALSE, help='Orient the edges in the learned network.')
parser <- add_option(parser, c('--walltime', '-w'), action='store_true', default=FALSE, help='Time the top level operations.')
args <- parse_args(parser, args=commandArgs(trailing=TRUE))

tRead <- proc.time()
data <- read.table(file=args$file, sep=args$separator, header=args$varnames)
if (args$colobs) {
        data <- as.data.frame(t(data))
}
data <- as.data.frame(lapply(data, as.factor))
#print(data)
tRead <- proc.time() - tRead
if (args$walltime) {
        cat('Time taken in reading the file:', tRead['elapsed'], 'sec\n')
}
library('bnlearn')
tNetwork <- proc.time()
network <- eval(parse(text=args$algorithm))(data, undirected=!args$directed)
tNetwork <- proc.time() - tNetwork

if (!is.null(args$target)) {
        tNeighborhood <- proc.time()
        neighbors <- c() 
        if (args$blanket) {
                neighbors <- mb(network, args$target)
        }
        else {
                neighbors <- pc(network, args$target)
        }
        tNeighborhood <- proc.time() - tNeighborhood
        if (args$walltime) {
                cat("Time taken in getting the neighborhod:", tNeighborhood['elapsed'], 'sec\n') 
        }
}

if (args$walltime) {
        cat('Time taken in getting the network:', tNetwork['elapsed'], 'sec\n')
}

if (!is.null(args$output)) {
        tWrite <- proc.time()
        write.dot(args$output, network)
        tNetwork <- proc.time() - tWrite 
        if (args$walltime) {
                cat('Time taken in writing the network:', tWrite['elapsed'], 'sec\n')
        }
}
