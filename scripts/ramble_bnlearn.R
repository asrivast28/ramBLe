#!/usr/bin/env Rscript

##
# @file ramble_bnlearn.R
# @brief Script for running bnlearn with the same options as ramble.
# @author Ankit Srivastava <asrivast@gatech.edu>
#
# Copyright 2020 Georgia Institute of Technology
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

cat('Using bnlearn from', find.package('bnlearn'), '\n')
library('parallel')
library('bnlearn')
library('optparse')


if (!exists("argv")) {
        argv = commandArgs(trailing=TRUE)
}

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
parser <- add_option(parser, c('--nprocs'), type='integer', default=1, help='Number of processes to use.')
parser <- add_option(parser, c('--ppn'), type='integer', default=16, help='Number of processes per node to use.')
args <- parse_args(parser, args=argv)

if (!args$learn && is.null(args$target) && is.null(args$output)) {
        cat("At least one of --target, --learn, or --output should be specified.\n")
        quit(status=1)
}

cl <- NULL
if (args$nprocs > 1) {
        tCluster <- proc.time()
        if (ceiling(args$nprocs / args$ppn) == 1) {
                cl <- makeCluster(args$nprocs)
        }
        else {
                nf <- Sys.getenv('PBS_NODEFILE')
                nodes <- readLines(nf)
                nnodes <- floor(args$nprocs / args$ppn)
                procs <- sort(rep(head(unique(nodes), n=nnodes), times=args$ppn))
                remaining <- args$nprocs - length(procs)
                if (remaining > 0) {
                        nodes <- c(nodes, rep(unique(nodes)[nnodes], times=remaining))
                }
                cat('Creating a PSOCK cluster with', length(procs), 'processes\n')
                cl <- makePSOCKcluster(procs)
        }
        tCluster <- proc.time() - tCluster
        cat('Created', show(cl), '\n')
        cat('Time taken in creating the cluster:', tCluster['elapsed'], '\n')
}

tRead <- proc.time()
data <- NULL
if (args$indices) {
        data <- read.table(file=args$file, sep=args$separator, check.names=!args$varnames, header=args$varnames, row.names=1)
} else {
        data <- read.table(file=args$file, sep=args$separator, check.names=!args$varnames, header=args$varnames)
}
if (args$colobs) {
        data <- as.data.frame(t(data), optional=args$varnames)
}
data <- as.data.frame(lapply(data, as.factor), optional=args$varnames)
tRead <- proc.time() - tRead
if (!((ncol(data) == args$nvars) && (nrow(data) == args$nobs))) {
        cat('Read dimensions:', nrow(data), 'x', ncol(data), '\n')
        stop('Read file did not match the expected dimensions.')
}
cat('Time taken in reading the file:', tRead['elapsed'], '\n')

network <- NULL
tNetwork <- NULL
if (args$learn || !is.null(args$target) || !is.null(args$output)) {
        tNetwork <- proc.time()
        network <- eval(parse(text=args$algorithm))(data, alpha=args$alpha, max.sx=args$conditioning, debug=!is.null(args$log), undirected=!args$directed, cluster=cl)
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
        tWrite <- proc.time() - tWrite
        cat('Time taken in writing the network:', tWrite['elapsed'], 'sec\n')
}
