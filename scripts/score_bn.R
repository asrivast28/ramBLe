#!/usr/bin/env Rscript

##
# @file score_bn.R
# @brief Script for scoring the given Bayesian network using bnlearn.
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
library('bnlearn')
library('optparse')
library('stringr')


if (!exists("argv")) {
        argv = commandArgs(trailing=TRUE)
}

parser <- OptionParser()
parser <- add_option(parser, c('--bn', '-b'), type='character', help='Name of the dot file from which the network is to be read.')
parser <- add_option(parser, c('--nvars', '-n'), type='integer', help='Number of variables in the dataset.')
parser <- add_option(parser, c('--nobs', '-m'), type='integer', help='Number of observations in the dataset.')
parser <- add_option(parser, c('--file', '-f'), type='character', help='Name of the file from which dataset is to be read.')
parser <- add_option(parser, c('--colobs', '-c'), action='store_true', default=FALSE, help='The file contains observations in columns.')
parser <- add_option(parser, c('--separator', '-s'), type='character', default=',', help='Delimiting character in the file.')
parser <- add_option(parser, c('--varnames', '-v'), action='store_true', default=FALSE, help='The file contains variable names.')
parser <- add_option(parser, c('--indices', '-i'), action='store_true', default=FALSE, help='The file contains observation indices.')
parser <- add_option(parser, c('--algorithm', '-a'), type='character', default='gs', help='Name of the algorithm used for generating the graph.')
parser <- add_option(parser, c('--alpha', '-p'), type='double', default=0.05, help='Threshold p-value used for generating the graph.')
parser <- add_option(parser, c('--loglevel'), type='character', help='Level of logging.')
args <- parse_args(parser, args=argv)

lines <- scan(args$bn, what='character', sep='\n')
nodes <- NULL
edges <- NULL
for (l in lines) {
        # Skip lines which are not statements
        if (!grepl(';', l)) {
                next
        }
        # Remove extraneous whitespaces
        # and quotes from the string
        l <- gsub('\"', '', str_squish(l), fixed=TRUE)
        s <- sub('\\s*(.*?)\\s*;', '\\1', l, perl=TRUE)
        if (startsWith(s, 'edge')) {
                s <- str_trim(sub('edge', '', s, fixed=TRUE))
                directed <- NULL
                if (startsWith(s, '[dir=none]')) {
                        s <- str_trim(sub('[dir=none]', '', s, fixed=TRUE))
                        directed <- FALSE
                }
                else if (startsWith(s, '[dir=forward]')) {
                        s <- str_trim(sub('[dir=forward]', '', s, fixed=TRUE))
                        directed <- TRUE
                }
                else {
                        stop('Unable to read the line:\n', l, '\n')
                }
                e <- paste(sub('(.*) -[>-]? (.*)', '\\1 \\2', s), directed)
                edges <- c(edges, e)
        }
        else { # Must be a node
                nodes <- c(nodes, s)
        }
}
cat('Read', length(nodes), 'nodes\n')
cat('Read', length(edges), 'edges\n')

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
if (!((ncol(data) == args$nvars) && (nrow(data) == args$nobs))) {
        cat('Read dimensions:', nrow(data), 'x', ncol(data), '\n')
        stop('Read file did not match the expected dimensions.')
}
if (!all(nodes == names(data))) {
    stop('The node names are different in the graph and the data file.')
}

net <- empty.graph(nodes)
for (e in edges) {
        edge <- unlist(strsplit(e, ' ', fixed=TRUE))
        if (edge[3] == 'FALSE') { # Undirected edge
                set.edge(net, c(edge[1]), edge[2], debug=!is.null(args$log))
        }
        else { # Directed edge
                set.arc(net, edge[1], edge[2], debug=!is.null(args$log))
        }
}
directed.net <- net
if (!directed(net)) {
        cat('The given network is not fully directed\n')
        cat('Directing it using the node ordering in the dataset\n')
        # Convert undirected edges to directed edges
        # assuming the same topological ordering as the
        # order of nodes in the dataset
        directed.net <- pdag2dag(net, nodes)
}
score.net <- score(directed.net, data)
cat('Score of the network is', score.net, '\n')
