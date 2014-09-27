\name{buildPhylo}
\alias{buildPhylo}
\title{Relations between subpopulations}
\description{Predicts phylogenetic relations between subpopulations from subpopulation specific copy number profiles.}
\usage{buildPhylo(ploidy, outF, treeAlgorithm="bionjs")}
\arguments{
  \item{ploidy}{
       Ploidy-matrix in which each row corrsponds to a copy number segment. Has to contain at least one column for each predicted subpopulation, with columnname labeled SP_xx, where xx is the size of the corresponding SP. Ploidy-matrix can be obtained by calling \code{\link{assignQuantityToSP}}. 
}
  \item{outF}{
       Prefix of file to which phylogeny will be saved.
}
  \item{treeAlgorithm}{
       Neighbor joining algorithm used for phylogeny reconstruction (from library ape). Options: bionjs (default), njs.
}
}
\details{
Reconstructs phylogenetic relationships between subpopulations using the BIONJ algorithm of Gascuel. Pairwise distances between subpopulations are calculated as the number of copy number segments for which both subpopulations have the same copy number, divided by the total number of copy number segments for which both subpopulations have available copy number information. Subpopulations with insufficient ploidy information are excluded from phylogeny.
}
\value{An object of class "phylo" (library ape).
}

\author{Noemi Andor}


\seealso{
\code{\link{assignQuantityToSP}}
}

