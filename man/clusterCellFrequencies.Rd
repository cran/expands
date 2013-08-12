\name{clusterCellFrequencies}
\alias{clusterCellFrequencies}
\title{Clustering of cellular frequency probability distributions}
\description{Calculates overrepresented cell frequencies using a two-step clustering procedure. Based on the assumption that passenger mutations occur within a cell prior to the driver event that initiates the expansion, each clonal expansion should be marked by multiple mutations. Thus mutations and copy number variations that took place in a cell prior to a clonal expansion should be present in a similar fraction of cells and leave a similar trace in the subsequent clonal expansion.}
\usage{clusterCellFrequencies(densities, precision, plotF=0, label=NA, nrep=30)}
\arguments{
  \item{densities}{
	 Matrix as obtained by \code{\link{computeCellFrequencyDistributions}.}Each row corresponds to a mutation and each column corresponds to a cellular frequency. Each value \eqn{densities[i,j]} represents the probability that mutation \eqn{i} is present in a fraction \eqn{f} of cells, where \eqn{f} is given by: \eqn{colnames(densities[,j]).}
}
  \item{precision}{
	 Precision with which subpopulation size is predicted, a small value reflects a high resolution and can trigger a higher number of predicted subpopulations (recommended: 0.1/log(n/7), where n = # mutations).
}
  \item{plotF}{
	 Value of 0 indicates no plot. If plotF > 0, a 3D plot will be generated to display the clustered probability distributions.
}
  \item{label}{
	 The sample name, used as title of the 3D plot. Only used if plotF > 0.
}
  \item{nrep}{
	 Positive integer giving the number of times to repeat the algorithm (default: 30).
}
}
\details{In the first step, mutations with similar cellular frequencies are grouped together by hierarchical cluster analysis of the probability distributions using the Kullback-Leibler divergence as a distance measure. The cell frequency   at each cluster-maxima denotes the size of the subpopulation that harbors the clustered mutations. In the second step, each cluster is extended by members with similar distributions in an interval around the cluster-maxima.}
\value{
	\item{SPs}{
	Matrix of predicted subpopulations. Each row corresponds to a subpopulation and each column contains information about that subpopulation, such as the size in the sequenced tumor bulk (column \bold{Mean Weighted}) and the confidence with which the subpopulation has been detected (column \bold{score}).}
}
\references{Noemi Andor, Julie Harness, Hans Werner Mewes and Claudia Petritsch. (2013) ExPANdS: Expanding Ploidy and Allele Frequency on Nested Subpopulations. Bioinformatics. In Review.}
\author{Noemi Andor}

