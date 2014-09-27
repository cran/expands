\name{runExPANdS}
\alias{runExPANdS}
\title{Main Function}
\description{Given a set of mutations, ExPANdS predicts the number of clonal expansions in a tumor, the size of the resulting subpopulations in the tumor bulk and which mutations accumulate in a cell prior to its clonal expansion. Input-parameters SNV and CBS hold the paths to tabdelimited files containing the mutations and the copy numbers respectively. Alternatively SNV and CBS can be read into the workspace and passed to \code{runExPANdS}  as numeric matrices. The robustness of the subpopulation predictions by ExPANdS increases with the number of mutations provided. It is recommended that SNV contains at least 200 mutations to obtain stable results.}
\usage{runExPANdS(SNV, CBS, maxScore=2.5, max_PM=6, min_CellFreq=0.1, precision=NA,
 plotF=2,snvF="out.expands",maxN=8000,region=NA)}
\arguments{
  \item{SNV}{
	Matrix in which each row corresponds to a mutation. Only mutations located on autosomes should be included. Columns in SNV must be labeled and must include: \cr
	\bold{chr} - the chromosome on which each mutation is located;\cr
	\bold{startpos} - the genomic position of each mutation; \cr
	\bold{AF_Tumor} - the allele-frequency of each mutation;\cr
	\bold{PN_B} - ploidy of B-allele in normal cells. A value of 0 indicates that the mutation has only been detected in the tumor sample (i.e. somatic mutation). A value of 1 indicates that the mutation is also present in the normal (control) sample, albeit at reduced allele frequency (i.e. mutation is consequence of LOH). Mutations, for which the allele frequency in the tumor sample is lower than the corresponding allele frequency in the normal sample, should not be included.
}
  \item{CBS}{
	Matrix in which each row corresponds to a copy number segment. CBS is typically the output of a circular binary segmentation algorithm. Columns in CBS must be labeled and must include:\cr
	\bold{chr} - chromosome; \cr
	\bold{startpos} - the first genomic position of a copy number segment;\cr 
	\bold{endpos} - the last genomic position of a copy number segment; \cr
	\bold{CN_Estimate} - the absolute copy number estimated for each segment.
}
  \item{maxScore}{
	Upper threshold for the confidence of subpopulation detection. Only subpopulations identified at a score below \eqn{maxScore} (default 2.5) are kept.
}
  \item{max_PM}{
	Upper threshold for the number of amplicons per mutated cell (default: 6). Increasing the value of this variable is not recommended unless extensive depth and breadth of coverage underly the measurements of copy numbers and allele frequencies. See also \code{\link{cellfrequency_pdf}.}
}
  \item{min_CellFreq}{
	Lower threshold for the prevalence of a mutated cell (default: 0.1); mutations for which allele frequency x copy number are below this value, are excluded from further computation. Decreasing the value of this variable is not recommended unless extensive depth and breadth of coverage underly the measurements of copy numbers and allele frequencies.
}
  \item{precision}{
	Precision with which subpopulation size is predicted, a small value reflects a high resolution and can trigger a higher number of predicted subpopulations (default 0.1/log(n/7), where n = # mutations).
}
  \item{plotF}{
	Option for displaying a visual representation of the identified SPs (0 - no display; 1 - display subpopulation size; 2 - display subpopulation size and phylogeny; default: 2).}
  \item{snvF}{
	The name of the file from which mutations have been read.}
  \item{maxN}{
	Upper limit for # SNVs during clustering. If number of user supplied SNVs exceeds <maxN>, the clustering of cellular frequency distributions will be restricted to SNVs found within <region> (default: 8000; increasing value of this parameter not recommended as complexity at least O(n^2)).}
  \item{region}{
	Regional boundary for mutations included during clustering. 
	Matrix in which each row corresponds to a genomic segment. Columns must include: \cr
	\bold{chr} - the chromosome of the segment ;\cr
	\bold{start} - the first genomic position of the segment;\cr
	\bold{end} - the last genomic position of the segment.\cr
	Default: SureSelectExome_hg19, comprising ca. 468 MB centered on the human exome. User supplied regions should also be coding regions, as the seletive pressure is higher and biological noise affecting cellular frequencies is reduced compared to non-coding regions.
}
}
\value{List with fields:
        \item{finalSPs}{Matrix of predicted subpopulations. Each row corresponds to a subpopulation and each column contains information about that subpopulation, such as the size in the sequenced tumor bulk (column \bold{Mean Weighted}) and the confidence with which the subpopulation has been detected (column \bold{score}).}
        \item{dm}{Matrix containing the input mutations with at least four additional columns:\cr
	\bold{SP} - the subpopulation to which the mutation has been asssigned; \cr
	\bold{\%maxP} - the confidence of assignment.
	\bold{f} - Deprecated. The maximum likelyhood cellular prevalence of this mutation, before it has been assigned to a SP. This value is based on the copy number and allele frequency of the mutation exclusively and is independent of other mutations. Column SP is less sensitive to noise and considered the more accurate estimation of cellular mutation prevalence.
	\bold{PM} - the total ploidy of all alleles at the mutated genomic locus, in the assigned subpopulation.
	\bold{PM_B} - the ploidy of the B-allele at the mutated genomic locus, in the assigned subpopulation.}	     
	\item{densities}{Matrix as obtained by \code{\link{computeCellFrequencyDistributions}.} Each row corresponds to a mutation and each column corresponds to a cellular frequency. Each value \eqn{densities[i,j]} represents the probability that mutation \eqn{i} is present in a fraction \eqn{f} of cells, where \eqn{f} is given by: \eqn{colnames(densities[,j]).}}
	\item{ploidy}{Matrix as obtained by \code{\link{assignQuantityToSP}.} Each row corresponds to a copy number segment as obtained by CBS. Includes one additional column for each predicted SP, holding the ploidy of each segment in the corresponding SP.}
	\item{tree}{An object of class "phylo" (library ape) as obtained by \code{\link{buildPhylo}.} Holds the inferred phylogenetic relationships between subpopulations.}
}
\references{Noemi Andor, Julie Harness, Sabine Mueller, Hans Werner Mewes and Claudia Petritsch. (2013) ExPANdS: Expanding Ploidy and Allele Frequency on Nested Subpopulations. Bioinformatics.}
\author{Noemi Andor}
\examples{
data(snv);
data(cbs);
maxScore=2.5;
set.seed(4); idx=sample(1:nrow(snv), 60, replace=FALSE);
#out= runExPANdS(snv[idx,], cbs, maxScore);
}
