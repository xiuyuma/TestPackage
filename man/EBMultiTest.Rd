\name{EBMultiTest}
\alias{EBMultiTest}
\title{
Using EM algorithm to calculate the posterior probabilities 
of interested patterns in a multiple condition study
}
\description{
'EBMultiTest' is built based on the assumption of NB-Beta Empirical Bayes model. It utilizes the EM algorithm to give the posterior probability of the interested patterns.
}
\usage{
EBMultiTest(Data, NgVector = NULL, Conditions, sizeFactors, uc, 
    Alpha = NULL, Beta = NULL, Qtrm = 1, QtrmCut = 0, maxround = 50, 
    step1 = 1e-06, step2 = 0.01, thre = log(2), sthre = 0.001, 
    filter = 10, stopthre = 0.001, nequal = 2)
}
\arguments{

  \item{Data}{A data matrix contains expression values for each transcript (gene or isoform level). In which rows should be transcripts and columns should be samples.}
  
  \item{NgVector}{A vector indicates the uncertainty group assignment of each isoform.
e.g. if we use number of isoforms in the host gene to define the uncertainty groups,
suppose the isoform is in a gene with 2 isoforms, Ng of this isoform should be 2.
The length of this vector should be the same as the number of rows in Data.
If it's gene level data, Ngvector could be left as NULL.}

  \item{Conditions}{A vector indicates the condition in which each sample belongs to. }

  \item{sizeFactors}{The normalization factors. It should be a vector with lane specific numbers (the length of the vector should be the same as the number of samples, with the same order as the columns of Data).}
  
  \item{uc}{number of unceratin positions, unit level}
  
  \item{Alpha}{start value of hyper parameter alpha}

  \item{Beta}{start value of hyper parameter beta}


\item{Qtrm, QtrmCut}{
Transcripts with Qtrm th quantile < = QtrmCut  will be removed before testing. The default value is Qtrm = 1 and QtrmCut=0.
By default setting, transcripts with all 0's
won't be tested.
}

  \item{maxround}{Number of iterations. The default value is 50.
Users should always check the convergency by looking at the Alpha and
Beta in output. If the hyper-parameter estimations are not converged
in 50 iterations, larger number is suggested.}

\item{step1}{stepsize for gradient ascent of alpha}

\item{step2}{stepsize for gradietn ascent of beta}

\item{thre}{threshold for determining the state of a position}

\item{sthre}{shrinkage threshold for iterative pruning during the EM updates}

\item{filter}{filterthreshold for low expression units}

\item{stopthre}{stopping threshold for EM}

\item{nequal}{when there is a chain of equal states with the number of equal states bigger than nequal, equalhandle algorithm will be used to further checking the homogeneity between the group means}
}


\value{
\item{Alpha}{Fitted parameter alpha of the prior beta distribution. }
\item{Beta}{Fitted parameter beta of the prior beta distribution.}
\item{P}{Global proportion of DE patterns.}
\item{RList}{The fitted values of r for each transcript.}
\item{MeanList}{The mean of each transcript (across conditions).}
\item{VarList}{The variance of each transcript (across conditions).}
\item{QList}{The fitted q values of each transcript within the two conditions}
\item{Mean}{The mean of each transcript within the two conditions (adjusted by normalization factors).}
\item{Var}{The estimated variance of each transcript within the two conditions (adjusted by normalization factors).}
\item{PoolVar}{The variance of each transcript (The pooled value of within condition EstVar).}
\item{DataNorm}{Normalized expression matrix.}
\item{AllZeroIndex}{The transcript with expression 0 for all samples (which are not tested).}
\item{PPMat}{The Posterior Probability of following each pattern (columns) for each transcript (rows). 
Transcripts with expression 0 for all samples are not shown in this matrix.}
\item{AllParti}{selected patterns}
\item{PPMatWith0}{The Posterior Probability of following each pattern (columns) for each transcript (rows). Transcripts with expression 0 for all samples are shown in this matrix with PP(any_pattrn)=NA.
The transcript order is exactly the same as the order of the input data.}
\item{Conditions}{The input conditions.}
\item{NumUC}{The number of uncertain positions at each unit}
}
\references{
Ning Leng, John A. Dawson, James A. Thomson, Victor Ruotti, Anna I. Rissman, Bart M.G. Smits, Jill D. Haag, Michael N. Gould, Ron M. Stewart, and Christina Kendziorski. EBSeq: An empirical Bayes hierarchical model for inference in RNA-seq experiments. Bioinformatics (2013)
}
\author{
Ning Leng, Xiuyu Ma
}

\seealso{
EBTest, GetMultiPP, GetMultiFC
}
\examples{
data(MultiGeneMat)
Conditions = c("C1","C1","C2","C2","C3","C3")
MultiSize = MedianNorm(MultiGeneMat)
MultiOut = EBMultiTest(MultiGeneMat,Conditions=Conditions,uc = 2,
                     sizeFactors=MultiSize)
MultiPP = GetMultiPP(MultiOut)
}
\keyword{ DE }
\keyword{ Multiple Condition }% __ONLY ONE__ keyword per line
