\name{EBTest}
\alias{EBTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Using EM algorithm to calculate the posterior probabilities of being DE
}
\description{
Base on the assumption of NB-Beta Empirical Bayes model, the EM algorithm is used to get the posterior probability of being DE.
}
\usage{
EBTest(Data, NgVector = NULL, Conditions, sizeFactors, uc,
    Alpha = NULL, Beta = NULL, Qtrm = 1, QtrmCut = 0, maxround = 50, 
    step1 = 1e-06, step2 = 0.01, thre = log(2), sthre = 0.001, 
    filter = 10, stopthre = 0.001) 
}
\arguments{

  \item{Data}{A data matrix contains expression values for each transcript (gene or isoform level). In which rows should be transcripts and columns should be samples.}
  
  \item{NgVector}{A vector indicates the uncertainty group assignment of each isoform. 
e.g. if we use number of isoforms in the host gene to define the uncertainty groups, suppose the isoform is in a gene with 2 isoforms, Ng of this isoform should be 2. The length of this vector should be the same as the number of rows in Data. If it's gene level data, Ngvector could be left as NULL.}

  \item{Conditions}{A factor indicates the condition which each sample belongs to. }
  
  \item{sizeFactors}{The normalization factors. It should be a vector with lane specific numbers (the length of the vector should be the same as the number of samples, with the same order as the columns of Data).}
  
  \item{uc}{number of unceratin positions, unit level}
  
  \item{Alpha}{start value of hyper parameter alpha}

  \item{Beta}{start value of hyper parameter beta}

\item{Qtrm, QtrmCut}{
Transcripts with Qtrm th quantile < = QtrmCut  will be removed before testing. The default value is Qtrm = 1 and QtrmCut=0.
By default setting, transcripts with all 0's 
won't be tested.
}

 \item{maxround}{Number of iterations. The default value is 50. Users should always check the convergency by looking at the Alpha and Beta in output. If the hyper-parameter estimations are not converged in 50 iterations, larger number is suggested.}

\item{step1}{stepsize for gradient ascent of alpha}

\item{step2}{stepsize for gradietn ascent of beta}

\item{thre}{threshold for determining the state of a position}

\item{sthre}{shrinkage threshold for iterative pruning during the EM updates}

\item{filter}{filterthreshold for low expression units}

\item{stopthre}{stopping threshold for EM}

}

\details{For each transcript gi within condition, the model assumes:
X_{gis}|mu_{gi} ~ NB (r_{gi0} * l_s, q_{gi})
q_gi|alpha, beta^N_g ~ Beta (alpha, beta^N_g)
In which the l_s is the sizeFactors of samples.

The function will test "H0: q_{gi}^{C1} = q_{gi}^{C2}" and "H1: q_{gi}^{C1} != q_{gi}^{C2}."
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
\item{PPMat}{A matrix contains posterior probabilities of being EE (the first column) or DE (the second column). 
Rows are transcripts. 
Transcripts with expression 0 for all samples are not shown in this matrix.}
\item{AllParti}{selected patterns}
\item{PPMatWith0}{A matrix contains posterior probabilities of being EE (the first column) or DE (the second column). 
Rows are transcripts. 
Transcripts with expression 0 for all samples are shown as PP(EE) = PP(DE) = NA in this matrix. 
The transcript order is exactly the same as the order of the input data.}
\item{Conditions}{The input conditions.
}
}
\references{
Ning Leng, John A. Dawson, James A. Thomson, Victor Ruotti, Anna I. Rissman, Bart M.G. Smits, Jill D. Haag, Michael N. Gould, Ron M. Stewart, and Christina Kendziorski. EBSeq: An empirical Bayes hierarchical model for inference in RNA-seq experiments. Bioinformatics (2013)
}
\author{
Ning Leng, Xiuyu Ma
}


\seealso{
EBMultiTest, PostFC, GetPPMat
}
\examples{
data(GeneMat)
str(GeneMat)
GeneMat.small = GeneMat[c(1:10,511:550),]
Sizes = MedianNorm(GeneMat.small)
EBOut = EBTest(Data=GeneMat, Conditions=as.factor(rep(c("C1","C2"),each=5)),
       uc = 1,sizeFactors = Sizes)
PP = GetPPMat(EBOut)
}
\keyword{ DE }
\keyword{ Two condition }% __ONLY ONE__ keyword per line
