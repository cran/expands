\name{cellfrequency_pdf}
\alias{cellfrequency_pdf}
\title{Computes the probability distribution of cellular frequencies for a single mutation.}
\description{Calculates \eqn{P} - the probability density distribution of cellular frequencies for one single mutation. For each \eqn{f}, the value of \eqn{P(f)} reflects the probability that the mutation is present in a fraction \eqn{f} of cells.
}
\usage{cellfrequency_pdf(af, cnv, pnb, freq, max_PM=6)}
\arguments{
  \item{af}{
	  The allelic frequency at which the muation has been observed.
}
  \item{cnv}{
	  The ploidy of the locus in which the mutation is embedded.
}
  \item{pnb}{
	  The ploidy of the B-allele in normal cells (binary variable: 1 if the mutation is a germline variant, 0 if somatic).
}
  \item{freq}{
	  Array of cellular frequencies at which the probabilities will be calculated.
}
  \item{max_PM}{
	  Upper threshold for the number of amplicons per mutated cell (default: 6). \eqn{max\_PM } is the maximum number of amplicons above which solutions are rejected in the cell-frequency estimation step described below, i.e. \eqn{PM} <= \eqn{max\_PM. } The choice of \eqn{max\_PM }should depend on genomic depth of coverage and on the fraction of the genome sequenced: the higher the quality and abundance of data, the higher \eqn{max\_PM.}
}
}
\details{
	We consider two types of molecular mechanisms that convert a locus into its mutated state: copy number variation (CNV) inducing events and single nucleotide variation (SNV) inducing events. We assume that a normal state is defined by a total ploidy of two and B allele ploidy below two, whereas a mutated state has an increased fraction of B alleles. The conditions defining these states for each locus \eqn{l} are as follows: i) \eqn{PM^B, PN^B, PM, PN \in N}; ii) \eqn{PM^B \geq 1; PN^B \leq 1; PN = 2}; iii) \eqn{\frac{PM^B}{PM} \geq \frac{PN^B}{PN} }. \cr	\cr         
\eqn{PM^B} and \eqn{PN^B}  denote the ploidy of the B allele in each cell type: mutated cells and normal cells, respectively. The value of \eqn{PN^B} is one if \eqn{l} has a germline variant, zero otherwise. \eqn{PM, PN}  are the total ploidy of mutated cells and normal cells.  \eqn{PM} is required to be between one and \eqn{max\_PM}, that is, we exclude solutions for which the maximum number of amplicons per cell exceeds the user defined constant \eqn{max\_PM}. \cr
The function returns the probability distribution, \eqn{P_l(f)}, that the mutation at locus \eqn{l} is present in a fraction \eqn{f} of cells, where \eqn{f \in [min_CellFreq,1.1]}. At default settings the interval starts at 0.1 because cellular frequencies below 0.1 are typically detected at very low allele-frequencies (<0.05), which in turn are often artifacts at moderate sequencing coverage.The interval ends at 1.1 because a local maxima at \eqn{f \approx 1.0} implies monotonically decreasing funtion-values before and after 1.0. So the extended threshold of 1.1 hat the purpose of allowing for this monotonically decreasing phase to the right of a peak close to 1.0.
}
\value{List with three components:
	  \item{p}{
	         The probability that the mutation is present in a fraction \eqn{f} of cells, for each input frequency \eqn{f}.
}         
	  \item{bestF}{
                 The cellular frequency that best explains the observed allele frequency and ploidities.
}
	   \item{errors}{
                Errors encountered during the density estimation step. 
}
}
\references{Noemi Andor, Julie Harness, Sabine Mueller, Hans Werner Mewes and Claudia Petritsch. (2013) ExPANdS: Expanding Ploidy and Allele Frequency on Nested Subpopulations. Bioinformatics.}
\author{Noemi Andor}
\examples{
freq=seq(0.1,1.1,by=0.01);
cfd=cellfrequency_pdf(af=0.26,cnv=2.13,pnb=0,freq=freq, max_PM=6)
plot(freq,cfd$p,type="l",xlab="f",ylab="P(f)");
}
