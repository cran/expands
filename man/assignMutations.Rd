\name{assignMutations}
\alias{assignMutations}
\title{Mutation Assignment}
\description{Assigns mutations to previously predicted subpopulations.}
\usage{assignMutations(dm, finalSPs, densities, max_PM=6)}
\arguments{
  \item{dm}{
       Matrix in which each row corrsponds to a mutation. Has to contain at least the following columnnames:\cr 
	\bold{chr} - the chrmosome on which each mutation is located; \cr
	\bold{startpos} - the genomic position of each mutation;\cr
	\bold{AF_Tumor} - the allele-frequency of each mutation;\cr
	\bold{PN_B} - the ploidy of the B-allele in normal (non-tumor) cells (binary variable: 1 if the mutation is a germline variant, 0 if somatic).
}
  \item{finalSPs}{
       Matrix in which each row corresponds to a subpopulation, as computed by \code{\link{clusterCellFrequencies}.}
}
  \item{densities}{
       The probability density distribution of cellular frequencies computed for each mutation via \code{\link{computeCellFrequencyDistributions}.} Has to contain the same number of rows as \code{dm.}
}
  \item{max_PM}{
       Upper threshold for the number of amplicons per mutated cell (default: 6). See also \code{\link{cellfrequency_pdf}.}
}
}
\details{
Each mutated locus \eqn{l} is assigned to the subpopulation \eqn{C}, whose size is closest to the maximum likelyhood cellular frequency of \eqn{l}:\cr
 \eqn{C:=argmin_C |argmax_f P_l(f)-f^C|}, where \eqn{P_l(f)} is the probability distribution of cellular frequencies as computed by \code{\link{cellfrequency_pdf}} and \eqn{f^C} is the size of subpopulation \eqn{C}. The mutated loci assigned to each subpopulation cluster represent the genetic profile of each predicted subpopulation.
}
\value{A list with two fields:
\item{dm}{
       The input matrix with two additional columns:\cr
	\bold{SP} - the subpopulation to which the mutation has been asssigned;\cr
	\bold{\%maxP} - confidence of the assignment.
} 
\item{finalSPs}{
       The input matrix of subpopulations with the column \bold{nMutations} updated according to the total number of mutations assigned to each subpopulation.
}
}

\author{Noemi Andor}


\seealso{
\code{\link{clusterCellFrequencies}}
}

