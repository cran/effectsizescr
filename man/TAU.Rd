\name{TAU}
\alias{TAU}
\title{Parker's TAU-U}
\usage{
TAU(data1, nameTime = "TIME", namePhase = "DUMMYPHASE",

	nameDV = "DV", Aphase = 0, Bphase = 1)
}
\arguments{
\item{data1}{a matrix or a dataframe with time,outcome and phases columns}

\item{nameTime}{the name of the time column}

\item{namePhase}{the name of the phase column}

\item{nameDV}{the name of the dependent variable}

\item{Aphase}{how phase A was coded in the phases column}

\item{Bphase}{how phase B was coded in the phases column}
}
\value{
A list including partition matrix, full matrix and TAU-U analysis
}
\description{
Nonoverlap and trend for single case research: the TAU-U function of Parker et al. (2011).
}
\references{
Parker, R. I., Vannest, K. J., Davis, J. L., & Sauber, S. B. (2011). Combining nonoverlap and trend for single-case research: Tau-U. Behavior Therapy, 42(2), 284-299, doi: 10.1016/j.beth.2010.08.006
}
\examples{
data=cbind(rnorm(16),1:16,c(rep(0,8),rep(1,8)))
colnames(data)=c("DV","TIME","PHASE")
TAU(data1=data,nameTime = "TIME",namePhase = "PHASE",
nameDV = "DV")
}
