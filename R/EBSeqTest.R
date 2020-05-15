

EBSeqTest <- function(data,conditions,uc, iLabel = 1,sizefactor = 1,
iter = 50,alpha = 0.4, beta = 0, step1 = 1e-6,step2 = 0.01,
thre = log(2), sthre = 0.001, filter = 10, stopthre = 1e-3, nequal = 2) {
    
    if(!is.matrix(data))
    {
        stop("data must be a numerical matrix")
    }
    if(length(conditions) != ncol(data))
    {
        stop("incorrect length of conditions")
    }
    if(length(beta) > 1)
    {
        if(length(beta) != nrow(data)){stop("incorrect length of hyper parameters")}
    }
    if(beta == 0)
    {
        beta = rep(2,nrow(data))
    }
    if(length(iLabel) == 1 && iLabel == 1)
    {
        iLabel = 1:nrow(data)
    }
    if(length(iLabel) != nrow(data))
    {
        stop("incorrect length of isoform label")
    }
    if(length(sizefactor) == 1 && sizefactor == 1)
    {
        sizefactor = rep(1,ncol(data))
    }
    if(length(sizefactor) != ncol(data))
    {
        stop("incorrect length of size factor")
    }
    
    
    
    .Call('EBSeq',
    scExpMatrix = data,
    groupLabel = conditions,
    iLabel = iLabel,
    sizeFactor = sizefactor,
    iter = iter,
    alpha = alpha,
    beta = beta,
    step1 = step1,
    step2 = step2,
    uc = uc,
    thre = thre,
    sthre = sthre,
    filter = filter,
    stopthre = stopthre,
    nequal = nequal)
    
    
}



