#!/usr/bin/env Rscript

##
# @file csl_bnlearn.R
# @brief Script for running bnlearn with the same options as csl.
# @author Ankit Srivastava <asrivast@gatech.edu>
#

library('optparse')
parser <- OptionParser()
parser <- add_option(parser, c('--nvars', '-n'), type='integer', help='Number of variables in the dataset.')
parser <- add_option(parser, c('--nobs', '-m'), type='integer', help='Number of observations in the dataset.')
parser <- add_option(parser, c('--file', '-f'), type='character', help='Name of the file from which dataset is to be read.')
parser <- add_option(parser, c('--colobs', '-c'), action='store_true', default=FALSE, help='The file contains observations in columns.')
parser <- add_option(parser, c('--separator', '-s'), type='character', default=',', help='Delimiting character in the file.')
parser <- add_option(parser, c('--varnames', '-v'), action='store_true', default=FALSE, help='The file contains variable names.')
parser <- add_option(parser, c('--indices', '-i'), action='store_true', default=FALSE, help='The file contains observation indices.')
parser <- add_option(parser, c('--algorithm', '-a'), type='character', default='gs', help='Name of the algorithm to be used.')
parser <- add_option(parser, c('--target', '-t'), type='character', help='Name of the target variable.')
parser <- add_option(parser, c('--blanket', '-b'), action='store_true', default=FALSE, help='Find MB instead of PC for the target var.')
parser <- add_option(parser, c('--learn', '-l'), action='store_true', default=FALSE, help='Force learn the network.')
parser <- add_option(parser, c('--output', '-o'), type='character', help='Name of the file to which the learned network should be written.')
parser <- add_option(parser, c('--directed', '-d'), action='store_true', default=FALSE, help='Orient the edges in the learned network.')
parser <- add_option(parser, c('--alpha', '-p'), type='double', default=0.05, help='Threshold p-value.')
parser <- add_option(parser, c('--conditioning', '-g'), type='integer', help='Maximum size of conditioning sets.')
parser <- add_option(parser, c('--log'), type='character', help='Level of logging.')
args <- parse_args(parser, args=commandArgs(trailing=TRUE))

if (!args$learn && is.null(args$target) && is.null(args$output)) {
        cat("At least one of --target, --learn, or --output should be specified.\n")
        quit(status=1)
}

tRead <- proc.time()
data <- NULL
if (args$indices) {
        data <- read.table(file=args$file, sep=args$separator, header=args$varnames, row.names=1)
} else {
        data <- read.table(file=args$file, sep=args$separator, header=args$varnames)
}
if (args$colobs) {
        data <- as.data.frame(t(data))
}
data <- as.data.frame(lapply(data, as.factor))
tRead <- proc.time() - tRead
if (!((ncol(data) == args$nvars) && (nrow(data) == args$nobs))) {
        cat('Read dimensions:', nrow(data), 'x', ncol(data), '\n')
        stop('Read file did not match the expected dimensions.')
}
cat('Time taken in reading the file:', tRead['elapsed'], '\n')

library('bnlearn')
network <- NULL
tNetwork <- NULL
if (args$learn || !is.null(args$target) || !is.null(args$output)) {
        tNetwork <- proc.time()
        network <- eval(parse(text=args$algorithm))(data, alpha=args$alpha, max.sx=args$conditioning, debug=!is.null(args$log), undirected=!args$directed)
        tNetwork <- proc.time() - tNetwork
}

if (!is.null(args$target)) {
        tNeighborhood <- proc.time()
        neighbors <- c()
        if (args$blanket) {
                neighbors <- mb(network, args$target)
        } else {
                neighbors <- pc(network, args$target)
        }
        tNeighborhood <- proc.time() - tNeighborhood
        cat("Time taken in getting the neighborhod:", tNeighborhood['elapsed'], '\n')
}

if (!is.null(tNetwork)) {
        cat('Time taken in getting the network:', tNetwork['elapsed'], '\n')
}

if (!is.null(args$output)) {
        tWrite <- proc.time()
        write.dot(args$output, network)
        tNetwork <- proc.time() - tWrite
        cat('Time taken in writing the network:', tWrite['elapsed'], 'sec\n')
}
