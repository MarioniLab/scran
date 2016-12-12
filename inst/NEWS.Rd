\name{scrannews}
\title{scran News}
\encoding{UTF-8}

\section{Version 1.3.3}{\itemize{
\item
Added run argument to selectorPlot().

\item
Minor bug fix to DM() when extrapolation is required.

\item
Moved spike-in specification to the scater package.

\item
Deprecated isSpike<- to avoid confusion over input/output types.

\item
Streamlined examples when mocking up data sets.
}}

\section{Version 1.2.0}{\itemize{
\item
Transformed correlations to a metric distance in quickCluster().

\item 
Removed normalize() in favour of scater's normalize().

\item 
Switched isSpike()<- to accept a character vector rather than a logical vector, to enforce naming of spike-in sets.
Also added warning code when the specified spike-in sets overlap.

\item
Allowed compute*Factors() functions to directly return the size factors.

\item 
Added selectorPlot() function for interactive plotting.

\item 
Switched to a group-based weighted correlation for one-way layouts in correlatePairs() and correlateNull(),
and to a correlation of residuals for more complex design matrices.

\item
Added phase assignments to the cyclone() output.

\item
Implemented Brennecke et al.'s method in the technicalCV2() function.

\item
Updated convertTo() to store spike-in-specific size factors as offsets.

\item
Moved code and subsetting into C++ to improve memory efficiency.

\item
Switched to loess-based trend fitting as the default in trendVar(), replaced polynomial with semi-loess fitting.

\item
Added significance statistics to output of decomposeVar(), with only the p-values replaced by NAs for spike-ins.

\item
Updated documentation and tests.
}}

\section{Version 1.0.0}{\itemize{
\item
New package scran, for low-level analyses of single-cell RNA sequencing data.
}}

