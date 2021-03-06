% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/name_conversion.R
\name{convert_gene_to_fbgn}
\alias{convert_gene_to_fbgn}
\title{Convert Gene Names to FlyBase Gene Numbers}
\usage{
convert_gene_to_fbgn(genes, db)
}
\arguments{
\item{genes}{a character vector containing gene names}

\item{db}{a list from \code{make_database()}}
}
\value{
a character vector containing FlyBase gene numbers
}
\description{
\code{convert_gene_to_fbgn()} takes a vector of gene names and query the
synonym table from FlyBase to convert them to a vector of FlyBase gene
numbers. Due to potential multiple mapping (one alias corresponding to many
FlyBase gene numbers, like Cha can be ChAT or Charaff), it's important to
note that when there's warning message, the vector
/code{convert_gene_to_fbgn()} returns will not retain the same order as the
input vector of gene names.
}
\examples{
# Prepare the database (using dummy local data)
dummypath <- system.file("extdata", "dummy.gtf",
package = "genofeatutil")
testdb <- make_database(species = "test",
gtf.path = dummypath)

# Convert gene names
query <- c("gene_1", "alias_2.1", "mt:dummy", "alias_3.3", "gene_21")
convert_gene_to_fbgn(query[1:4], db = testdb)
}
