\name{assignQuantityToSP}
\alias{assignQuantityToSP}
\title{Quantity assignment (ploidy) to subpopulations}
\description{Assigns quantities to predicted SPs. Currently, the only assignable quantity are SP specific ploidies for the input genome segments (obtained by CBS).
}
\usage{assignQuantityToSP(cbs, dm)}
\arguments{
  \item{dm}{
	Matrix in which each row corrsponds to a mutation. Has to contain at least the following columnnames: \cr 
	\bold{chr} - the chrmosome on which each mutation is located;\cr
	\bold{startpos} - the genomic position of each mutation;\cr
	\bold{SP} - the SP to which the mutation has been assigned;\cr
	\bold{PM} - the ploidy of the SP for the segment in which the mutation is embedded.
}
  \item{cbs}{
	Matrix in which each row corrsponds to a copy number fragment as computed by a circular binary segmentation algorithm. Has to contain at least the following columnnames:\cr
	\bold{chr} - chromosome;\cr
	\bold{startpos} - the first genomic position of a copy number segment; \cr
	\bold{endpos} - the last genomic position of a copy number segment; \cr
	\bold{CN_Estimate} - the copy number estimated for each segment. 
}
}

\value{
\item{cbs}{
       The input matrix with one additional column for each predicted SP:\cr
	\bold{SP_size} - the ploidy of each segment in the corresponding SP; Value <NA> indicates that no ploidy could be inferred for the segment in the corresponding SP (either because the SP had zero mutations within the segment, or because the SP had multiple, non-uniform mutations within the segment).
}
}

\author{Noemi Andor}
