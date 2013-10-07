\name{snv}
\alias{snv}
\docType{data}
\title{Single Nucleotide Variations}
\description{Somatic mutations and Loss of Heterozygosity (LOH) of a Glioblastoma tumor (TCGA-06-0152-01)}
\usage{data(snv)}
\format{
	Numeric matrix with 773 rows (one per mutation) and 7 columns:\cr
	 \bold{chr} - the chromosome \cr
	 \bold{startpos} - genomic position \cr
	 \bold{endpos} - same as above \cr
	 \bold{REF} - ASCII code of the reference nucleotide (in hg18/hg19) \cr
	 \bold{ALT} - ASCII code of the B-allele nucleotide \cr
	 \bold{AF_Tumor} - allele frequeny of B-allele \cr
	 \bold{PN_B} - ploidy of B-allele in normal cells (1 if B-allele is a germline variant, 0 if the B-allele is a consequence of a somatic mutation)
}

\source{Data derived from The Cancer Genome Atlas (TCGA)}


\keyword{datasets}
