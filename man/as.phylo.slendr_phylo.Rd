% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree-sequences.R
\name{as.phylo.slendr_phylo}
\alias{as.phylo.slendr_phylo}
\title{Convert an annotated \code{slendr_phylo} object to a \code{phylo} object}
\usage{
\method{as.phylo}{slendr_phylo}(x)
}
\arguments{
\item{x}{Tree object of the class \code{slendr_phylo}}
}
\description{
This function servers as a workaround around a ggtree error:
\code{Error in UseMethod("as.phylo") :
    no applicable method for 'as.phylo' applied to an object of class
    "c('phylo', 'slendr_phylo')"}
}
