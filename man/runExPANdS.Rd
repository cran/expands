\name{runExPANdS}
\alias{runExPANdS}
\title{Main Function}
\description{Given a set of mutations, ExPANdS predicts the number of clonal expansions in a tumor, the size of the resulting subpopulations in the tumor bulk and which mutations accumulate in a cell prior to its clonal expansion. Input-parameters SNV and CBS hold the paths to tabdelimited files containing the mutations and the copy numbers respectively. Alternatively SNV and CBS can be read into the workspace and passed to \code{runExPANdS  } as a numeric matrix. The robustness of the subpopulation predictions by ExPANdS increases with the number of mutations provided. It is recommended that SNV contains at least 200 mutations to obtain stable results.}
\usage{runExPANdS(SNV, CBS, maxScore=2.5, max_PM=6, precision=NA, plotF=1,snvF="out.expands")}
\arguments{
  \item{SNV}{
	Matrix in which each row corresponds to a mutation. Columns in SNV must include: \cr
	\bold{chr} - the chrmosome on which each mutation is located;\cr
	\bold{startpos} - the genomic position of each mutation; \cr
	\bold{AF_Tumor} - the allele-frequency of each mutation;\cr
	\bold{PN_B} - the ploidy of the B-allele in normal cells (binary variable: 1 if the SNV is a germline variant, 0 if somatic).
}
  \item{CBS}{
	Matrix in which each row corresponds to a copy number segment. CBS is typically the output of a circular binary segmentation algorithm. Columns in CBS must include:\cr
	\bold{chr} - chromosome; \cr
	\bold{startpos} - the first genomic position of a copy number segment;\cr 
	\bold{endpos} - the last genomic position of a copy number segment; \cr
	\bold{CN_Estimate} - the copy number estimated for each segment.
}
  \item{maxScore}{
	Upper threshold for the confidence of subpopulation detection. Only subpopulations identified at a score below \eqn{maxScore} (default 2.5) are kept.
}
  \item{max_PM}{
	Upper threshold for the number of amplicons per mutated cell (default: 6). Increasing the value of this variable is not recommended unless extensive depth and breadth of coverage underly the measurements of copy numbers and allele frequencies. See also \code{\link{cellfrequency_pdf}.}
}
  \item{precision}{
	Precision with which subpopulation size is predicted, a small value reflects a high resolution and can trigger a higher number of predicted subpopulations (default 0.1/log(n/7), where n = # mutations).
}
  \item{plotF}{
	Option for displaying a visual representation of the identified SPs (0 - no display; 1 - display subpopulation size; 2 - display subpopulation size and cell-frequency probability clusters; default: 1).}
  \item{snvF}{
	The name of the file from which mutations have been read.}
}
\value{List with three fields:
        \item{finalSPs}{Matrix of predicted subpopulations. Each row corresponds to a subpopulation and each column contains information about that subpopulation, such as the size in the sequenced tumor bulk (column \bold{Mean Weighted}) and the confidence with which the subpopulation has been detected (column \bold{score}).}
        \item{dm}{Matrix containing the input mutations with at least two additional columns:\cr
	\bold{SP} - the subpopulation to which the mutation has been asssigned; \cr
	\bold{\%maxP} - the confidence of assignment.}
	\item{densities}{Matrix as obtained by \code{\link{computeCellFrequencyDistributions}.} Each row corresponds to a mutation and each column corresponds to a cellular frequency. Each value \eqn{densities[i,j]} represents the probability that mutation \eqn{i} is present in a fraction \eqn{f} of cells, where \eqn{f} is given by: \eqn{colnames(densities[,j]).}}
}
\references{Noemi Andor, Julie Harness, Hans Werner Mewes and Claudia Petritsch. (2013) ExPANdS: Expanding Ploidy and Allele Frequency on Nested Subpopulations. Bioinformatics. In Review.}
\author{Noemi Andor}
\examples{
data(snv);
data(cbs);
maxScore=2.5;
set.seed(4); idx=sample(1:nrow(snv), 60, replace=FALSE);
#out= runExPANdS(snv[idx,], cbs, maxScore);
}
