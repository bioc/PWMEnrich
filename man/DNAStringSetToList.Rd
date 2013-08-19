\name{DNAStringSetToList}
\alias{DNAStringSetToList}
\title{Convert DNAStringSet to list of DNAString objects}
\usage{
  DNAStringSetToList(x)
}
\arguments{
  \item{x}{an object of class DNAStringSet}
}
\description{
  as.list doesn't seem to always work for DNAStringSets, so
  implementing this ourselves.
}

