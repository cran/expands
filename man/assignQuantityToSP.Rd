\name{assignQuantityToSP}
\alias{assignQuantityToSP}
\title{Quantity assignment (ploidy) to subpopulations}
\description{Assigns quantities to predicted SPs. Currently, the only assignable quantity are SP specific ploidies for the input genome segments (obtained by CBS).
}
\usage{assignQuantityToSP(cbs, dm, colName = "PM", keepAmbigSeg=FALSE)}
\arguments{
  \item{dm}{
	Matrix in which each row corrsponds to a mutation. Has to contain at least the following columnnames: \cr 
	\bold{chr} - the chromosome on which each mutation is located;\cr
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
  \item{colName}{
	The SP specific value assigned to each copy number fragment. Possible values: PM, PM_B. Default: PM.
}
\item{keepAmbigSeg}{
        Whether to assign ploidy to a subpopulation, SP_i, for a segment containing multiple SP_i specific SNVs, at least two of which have distinct ploidies. If set to TRUE, the median SNV ploidy is assigned as segment ploidy. Setting this parameter to TRUE is not recommended as the output will include segment-assignements where subpopulation specific ploidy is ambiguous.Recommend repeating circular binary segmentation with less stringent parameters instead, to reduce segment length and thus the prevalence of ambiguous assignements. Default: FALSE.
}
}

\value{
\item{cbs}{
       The input matrix with one additional column for each predicted SP:\cr
	\bold{SP_size} - the ploidy of each segment in the corresponding SP; Value <NA> indicates that no ploidy could be inferred for the segment in the corresponding SP (either because the SP had zero mutations within the segment, or because the SP had multiple, ambiguous mutation-ploidies within the segment).
}
}

\author{Noemi Andor}
