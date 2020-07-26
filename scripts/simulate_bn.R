#!/usr/bin/env Rscript

##
# @file simulate_bn.R
# @brief Script for simulating a Bayesian network using pcalg 
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

library('pcalg')
library('optparse')


if (!exists('argv')) {
        argv = commandArgs(trailing=TRUE)
}

parser <- OptionParser()
parser <- add_option(parser, c('--seed'), type='integer', help='PRNG seed.')
parser <- add_option(parser, c('--nvars', '-n'), type='integer', help='Number of variables in the simulated network.')
parser <- add_option(parser, c('--prob', '-p'), type='double', default=0.05, help='Threshold p-value used for generating the graph.')
parser <- add_option(parser, c('--mbsize'), action='store_true', default=FALSE, help='Print out the average MB size.')
parser <- add_option(parser, c('--bn', '-b'), type='character', help='Name of the dot file to which the network is written.')
parser <- add_option(parser, c('--nobs', '-m'), type='integer', help='Number of observations in the simulated dataset.')
parser <- add_option(parser, c('--datafile', '-d'), type='character', help='Name of the file to which dataset is written.')
parser <- add_option(parser, c('--colobs', '-c'), action='store_true', default=FALSE, help='The file contains observations in columns.')
parser <- add_option(parser, c('--sep', '-s'), type='character', default=' ', help='Delimiting character in the dataset file.')
args <- parse_args(parser, args=argv)


if (!is.null(args$seed)) {
        set.seed(args$seed)
}

tGenerate <- proc.time()
nodes <- c()
for (v in seq(1, args$nvars)) {
        nodes <- c(nodes, paste('V', v, sep=''))
}
# Graph of type graphNEL
dag <- randomDAG(args$nvars, prob=args$prob, V=nodes)
show(dag)
tGenerate <- proc.time() - tGenerate
cat('Time taken in generating the network:', tGenerate['elapsed'], 'sec\n')

if (args$mbsize) {
        cat('Using bnlearn from', find.package('bnlearn'), '\n')
        library('bnlearn')
        tMBSize <- proc.time()
        # Convert to BN
        bn <- as.bn(dag)
        avgmb <- mean(sapply(nodes, function(n) { length(bn$nodes[[n]]$mb) }))
        cat('Average MB size is', avgmb, '\n')
        tMBSize <- proc.time() - tMBSize
        cat('Time taken in getting the MB sizes:', tMBSize['elapsed'], 'sec\n')
}

if (!is.null(args$bn)) {
        tWrite <- proc.time()
        write.dot(args$bn, bn)
        tWrite <- proc.time() - tWrite
        cat('Time taken in writing the network:', tWrite['elapsed'], 'sec\n')
}

if (!is.null(args$nobs)) {
        tData <- proc.time()
        data <- rmvDAG(args$nobs, dag, use.node.names=TRUE)
        if (!args$colobs) {
                data <- t(data)
        }
        tData <- proc.time() - tData
        cat('Time taken in getting the data:', tData['elapsed'], 'sec\n')
        tWrite <- proc.time()
        write.table(data, file=args$datafile, sep=args$sep, row.names=!args$colobs, col.names=args$colobs)
        tWrite <- proc.time() - tWrite
        cat('Time taken in writing the dataset:', tWrite['elapsed'], 'sec\n')
}
