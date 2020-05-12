#pragma once


#include "EBSeq.hpp"
#include "agglomerativeClustering.hpp"
#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <set>
//#include <tbb/tbb.h>

namespace EBS
{
  
    class NB:public EBSeq
    {
    public:
        
        
        NB(COUNTS& scRNAexpMatrix, std::vector<int>& cellCluster, Eigen::VectorXd& sizeFactor) : EBSeq(scRNAexpMatrix, cellCluster, sizeFactor)
        {
            size_t nk = (_clusinfo.size).size();
            
            if(scRNAexpMatrix.cols() > nk)
            {
                // method of moments to estimate size factor r of NB
                _var = aggregate::groupVar(scRNAexpMatrix, _mean, _clusinfo, sizeFactor);
                
                size_t G = scRNAexpMatrix.rows();
                
                _poolVar = _var.rowwise().mean();
                
                COUNTS q(G,1), I(G,1), mn;
                
                I.fill(1);
                
                mn = _mean.rowwise().mean();
                
                _q.resize(G,_mean.cols());
                
                for(size_t i = 0; i < G; i++){
                    
                    if(fabs(_poolVar(i,0) - 0) < 0.0001)
                        _poolVar(i,0) = 1;
                    
                    if(_poolVar(i,0) <= mn(i,0))
                        q(i,0) = 0.99;
                    else
                        q(i,0) = mn(i,0) / _poolVar(i,0);
                    
                    for(size_t j = 0; j < _mean.cols(); j++)
                    {
                        _q(i,j) = _mean(i,j) / _var(i,j);
                    }
                }
                
                
                // aggregate sum of size factor
                auto ssum = aggregate::sum(sizeFactor, _clusinfo);
                
                _r = ((mn.cwiseProduct(q)).array() / (I - q).array()).matrix();
                
                _rsum =  _r * ssum;
                
            }else
            {
                
                Eigen::ArrayXd use(_sum.rows());
                
                use.fill(1);
                
                Eigen::ArrayXd add(_sum.rows());
                
                add.fill(0);
                
                // case when there is no replicates
                for(size_t i = 0; i < nk - 1; i++)
                {
                    for(size_t j = i + 1; j < nk; j++)
                    {
                        Eigen::ArrayXd ratioFilter = (_sum.col(i)).array() / (_sum.col(j)).array();
                        
                        std::vector<Float> rmNA;
                        
                        for(size_t iter = 0; iter < _sum.rows(); iter++)
                        {
                            if(!std::isnan(ratioFilter[iter]))
                                rmNA.push_back(ratioFilter[iter]);
                        }
                        
                        if(rmNA.empty())
                        {
                            continue;
                        }
                        
                        size_t Q1 = rmNA.size() / 4;
                        size_t Q2 = rmNA.size() * 3 / 4;
                        
                        std::nth_element(rmNA.begin(),          rmNA.begin() + Q1, rmNA.end());
                        
                        std::nth_element(rmNA.begin() + Q1 + 1,          rmNA.begin() + Q2, rmNA.end());
                        
                        Float q25 = rmNA[Q1];
                        
                        Float q75 = rmNA[Q2];
                        
                        for(size_t iter = 0; iter < _sum.rows(); iter++)
                        {
                            if(std::isnan(ratioFilter[iter]) || ratioFilter[iter] < q25 || ratioFilter[iter] > q75)
                            {
                                ratioFilter[iter] = 0;
                            }
                            else
                            {
                                ratioFilter[iter] = 1;
                            }
                        }
                        
                        use *= ratioFilter;
                        
                        add += ratioFilter;
                    }
                }
                
                
                if(use.sum() == 0)
                {
                    //use top 3 most frequent valid genes
                    Eigen::ArrayXd::Index maxRow, maxCol;
                    
                    auto maxholder = add.maxCoeff(&maxRow, &maxCol);
                    
                    use(maxRow,maxCol) = 1;
                    
                    add(maxRow,maxCol) = 0;
                    
                    maxholder = add.maxCoeff(&maxRow, &maxCol);
                    
                    use(maxRow,maxCol) = 1;
                    
                    add(maxRow,maxCol) = 0;
                    
                    maxholder = add.maxCoeff(&maxRow, &maxCol);
                    
                    use(maxRow,maxCol) = 1;
                }
                
                COUNTS dataUse(int(use.sum()),_sum.cols());
                
                size_t rowCount = 0;
                
                for(size_t useRow = 0; useRow < use.size(); useRow++)
                {
                    if(use[useRow]>0)
                    {
                        dataUse.row(rowCount) = _sum.row(useRow);
                        rowCount++;
                    }
                }
                
                Eigen::VectorXd mean_use = dataUse.rowwise().mean();
                
                COUNTS var_use = (dataUse.colwise() - mean_use).rowwise().squaredNorm() / (nk - 1);
                
                auto wholeMean = _sum.rowwise().mean().array();
                
                auto phi = (var_use - mean_use).array() / (mean_use.array() * mean_use.array());
                
                Float rmNEG = 0;
                int cnn = 0;
                
                for(size_t iter = 0; iter < phi.size(); iter++)
                {
                    if(phi(iter) >= 0)
                    {
                        rmNEG += phi(iter);
                        cnn++;
                    }
                }
                
                Float phiUse = rmNEG / cnn;
                
                auto varEst = wholeMean * (1 + wholeMean * phiUse);
                
                _q = (wholeMean / varEst).matrix();
                
                COUNTS I(_sum.rows(),1);
                
                I.fill(1);
                
                _r = ((wholeMean.matrix().cwiseProduct(_q)).array() / (I - _q).array()).matrix();
                
                auto ssum = aggregate::sum(sizeFactor, _clusinfo);
                
                _rsum = _r * ssum;
                
            }
            
            // mean and sum needs to be transformed back to raw
            // size factor only affects r
            
            _sum = aggregate::sum(scRNAexpMatrix, _clusinfo);
            
        }
        
        NB(COUNTS& scRNAexpMatrix, std::vector<int>& cellCluster, Eigen::VectorXd& sizeFactor, Eigen::VectorXd& r) : EBSeq(scRNAexpMatrix, cellCluster, sizeFactor)
        {
            // aggregate sum of size factor
            auto ssum = aggregate::sum(sizeFactor, _clusinfo);
            
            _rsum = r * ssum;
        }
        
        void init(Float alpha, Eigen::VectorXd beta, std::vector<int> iLabel, std::vector<Float> lrate, int UC, Float thre, Float sthre, Float filter, int nequal)
        {
            // uncertatinty should be smaller than number of subtypes - 1 (number of in/equalities)
            assert(UC < _sum.cols());
            
            // scalar of alpha for beta prior, shared by genome
            _alpha = alpha;
            
            // vector of beta for beta prior, gene specific
            _beta = beta;
            
            // label for isoform.
            _isoLabel = iLabel;
            
            // length of equal chain that need check
            _nequal = nequal;
            
            _ng = *(std::max_element(iLabel.begin(),iLabel.end()));
            
            // only when there is nontrivial isoform label, then init isoPos
            if(_ng < _sum.rows())
            {
                for(size_t lab = 1; lab < _ng + 1; lab++)
                {
                    std::vector<int> tmpPos;
                    
                    for(size_t pos = 0; pos < _isoLabel.size(); pos++)
                    {
                        if(_isoLabel[pos] == lab)
                        {
                            tmpPos.push_back(pos);
                        }
                    }
                    
                    _isoPos.push_back(tmpPos);
                }
            }
            
            // stepsize of gradient updates for alpha and beta
            _lrate = lrate;
            
            // number of uncertatinty of in/equalities between subtypes
            _uncertainty = UC;
            
            // threshold for distinguish certain and uncertain
            _threshold = thre;
            
            // threshold for shrinkage
            _shrinkthre = sthre;
            
            // filter threshold for low expression subtypes
            _filter = filter;
            
            // resize mde matrix as K by K
            //_mde.resize(_sum.cols(),_sum.cols());
            
            //_mde.fill(1.0 / 2);
            
            // get promising DE pattern
            DEpat();
            
            // error checking, number of promising DE patterns must > 0
            size_t n = _dep.size();
            
            assert(n > 0);
            
            // inita prop for each DE pattern with equal proportion
            _p.resize(n);
            
            _p.fill(1.0 / n);
            
        }
        
        inline Float lbeta(Float x,Float y)
        {
            return boost::math::lgamma(x) + boost::math::lgamma(y) - boost::math::lgamma(x + y);
        }
        
        inline COUNTS lbeta(COUNTS& A, COUNTS& B)
        {
            return A.unaryExpr<Float(*)(Float)>(&boost::math::lgamma) + B.unaryExpr<Float(*)(Float)>(&boost::math::lgamma) - (A + B).unaryExpr<Float(*)(Float)>(&boost::math::lgamma);
        }
    
        Float OBJ(Eigen::VectorXd& p)
        {
            setP(p);
            
            return (_kernel * _p).sum();
        }
        
//        void oneRunUpdate()
//        {
//            //updateMDE();
//
//            // then given p and dep update alpha and beta
//            gradientAscent();
//
//            kernel();
//
//            posterior();
//
//            updateP();
//        }
        
        
        size_t DEPsize()
        {
            return _dep.size();
        }
        
        std::vector<std::vector<int>> getDEP()
        {
            return _dep;
        }
        
        COUNTS getR()
        {
            return _r;
        }
        
        
        COUNTS getVar()
        {
            return _var;
        }
        
        COUNTS getQ()
        {
            return _q;
        }
        
        COUNTS getPoolVar()
        {
            return _poolVar;
        }
        
        std::vector<size_t> getGUC()
        {
            return _guc;
        }
        
        std::vector<COUNTS> getPAT()
        {
            return _pat;
        }
        
        
        COUNTS getKernel()
        {
            return _kernel;
        }
        
//        COUNTS getMDE()
//        {
//            return _mde;
//        }
        
        Float getALP()
        {
            return _alpha;
        }
        
        Eigen::VectorXd getBETA()
        {
            return _beta;
        }
        
        COUNTS getPOST()
        {
            return _post;
        }
        
        Float getOBJ()
        {
            return OBJ(_p);
        }
        
        void setAlphaBeta(Float alpha, Eigen::VectorXd& beta)
        {
            _alpha = alpha;
            
            _beta = beta;
        }
        
        void setP(Eigen::VectorXd& p)
        {
            _p = p;
        }
        
        Eigen::VectorXd getP()
        {
            return _p;
        }
        
    private:
        
        void Estep()
        {
            // update kernel given alpha and beta
            kernel();
            
            // update posterior with p
            posterior();
        }
        
        void Mstep()
        {
            // update alpha and beta
            gradientAscent();

            // update p
            updateP();
        }
        
        // only to be called in init genelevel depat
        
        typedef Eigen::Block<COUNTS,1> ROW;
//        std::vector<std::vector<int>> DEpatGene(ROW&& mean, ROW&& cs, ROW&& rs,
//                                                std::vector<int>& sz, size_t& K,
//                                                Float& alpha, Float& beta, Float& filter,
//                                                Float& threshold, int& uncertainty)
//        {
//
//            std::vector<Float> abslogRatio(K - 1);
//
//            std::vector<int> baseClus(K);
//
//            std::vector<std::vector<int>> res;
//
//            auto ord = helper::sortIndexes<ROW>(mean);
//
//            auto ord2 = helper::sortIndexes2<ROW>(mean);
//
//            baseClus[0] = 1;
//
//            for(size_t j = 1; j < K; j++)
//            {
//                // position for j-1th and jth subtypes
//                auto o1 = ord[j - 1];
//
//                auto o2 = ord[j];
//
//                Float tmp = kernel2case(cs(o1),cs(o2),rs(o1),rs(o2)
//                                        ,sz[o1],sz[o2],alpha,beta,filter);
//
//                abslogRatio[j - 1] = abs(tmp);
//
//                //  more favorable for equal mean
//                if(tmp > 0)
//                {
//                    baseClus[j] = baseClus[j - 1];
//                }
//                else
//                {
//                    // DE start a new cluster
//                    baseClus[j] = baseClus[j - 1] + 1;
//                }
//
//            }
//
//
//            auto tmpOrd = helper::sortIndexes<std::vector<Float>>(abslogRatio);
//
//            auto baseBit = partition::mapToBit(baseClus);
//
//
//
//            equalHandle(baseBit, abslogRatio, cs, rs,
//                        alpha, beta, _threshold, _filter, sz);
//
//
//
//
//            int localUC = 0;
//
//            for(auto score:abslogRatio)
//            {
//                if(score < threshold)
//                {
//                    localUC++;
//                }
//            }
//
////            _guc.push_back(localUC);
//
//            localUC = (localUC < uncertainty) ? localUC : uncertainty;
//
//            if(localUC < 1)
//            {
//
//                auto newClusOr = baseClus;
//
//                for(size_t iter = 0; iter < ord2.size(); iter++)
//                {
//                    newClusOr[iter] = baseClus[ord2[iter]];
//                }
//
//                auto newClusOrd = partition::reorder(newClusOr);
//
//                res.push_back(newClusOrd);
//
//                return res;
//            }
//
//            auto pBit = partition::genBit(localUC);
//
//            // get promising DE pattern
//            for(auto x:pBit)
//            {
//
//                auto newBit = baseBit;
//
//                for(int t = 0; t < _uncertainty; t++)
//                {
//                    auto tmpJ = tmpOrd[t];
//
//                    newBit[tmpJ] = x[t];
//                }
//
//                auto newClus = partition::bitToPart(newBit);
//
//                auto newClusOr = newClus;
//
//                for(size_t iter = 0; iter < ord2.size(); iter++)
//                {
//                    newClusOr[iter] = newClus[ord2[iter]];
//                }
//
//                std::vector<int> newClusOrd = partition::reorder(newClusOr);
//
////                auto sClus = partition::toString<std::vector<int>>(newClusOrd);
//
////                if(dep.count(sClus) < 1)
////                {
////                    dep.insert(sClus);
////
////                    _dep.push_back(newClusOrd);
////
////                    _pat.push_back(partition::toMatrix(newClusOrd));
////                }
//
//                res.push_back(newClusOrd);
//
//            }
//            return res;
//        }
        
        // only to be called in init
        void DEpat()
        {
            size_t G = _mean.rows();
            
            size_t K = _mean.cols();
            
            // DE patterns to be considered
            std::set<std::string> dep;
            
//            tbb::parallel_for(blocked_range<size_t>(0,G), DEpatGene());
//
//            for(size_t i = 0; i < G; i++)
//            {
//
//                Rcpp::Rcout << "gene " << i << "\n";
//
//                auto res = DEpatGene(_mean.row(i), _sum.row(i), _rsum.row(i),
//                                     _clusinfo.size, K, _alpha, _beta(i),_filter,
//                                     _threshold,_uncertainty);
//
//                for(auto newClusOrd:res)
//                {
//                    std::string sClus = partition::toString<std::vector<int>>(newClusOrd);
//
//                    if(dep.count(sClus) < 1)
//                    {
//                        dep.insert(sClus);
//
//                        _dep.push_back(newClusOrd);
//
//                        _pat.push_back(partition::toMatrix(newClusOrd));
//                    }
//                }
//            }
            
            std::vector<Float> abslogRatio(K - 1);
            
            std::vector<int> baseClus(K);
          
            for(size_t i = 0; i < G; i++)
            {
                auto ord = helper::sortIndexes<ROW>(_mean.row(i));
                
                auto ord2 = helper::sortIndexes2<ROW>(_mean.row(i));
                
                baseClus[0] = 1;
                
                for(size_t j = 1; j < K; j++)
                {
                    // position for j-1th and jth subtypes
                    auto o1 = ord[j - 1];
                    
                    auto o2 = ord[j];
                    
                    Float s1 = _sum(i,o1);
                    
                    Float s2 = _sum(i,o2);
                    
                    int n1 = _clusinfo.size[o1];
                    
                    int n2 = _clusinfo.size[o2];
                    
                    Float r1 = _rsum(i,o1);
                    
                    Float r2 = _rsum(i,o2);
                    
                    // ratio of EE and DE prior predictive functions
                    //                    Float use, marginal;
                    //
                    //                    if(o1 < o2)
                    //                        use = _mde(o1,o2);
                    //                    else
                    //                        use = _mde(o2,o1);
                    
                    //                    marginal = use / (1 - use);
                    
                    //                    Float tmp = kernel2case(s1,s2,r1,r2,n1,n2,i) + log(marginal);
                    
                    Float tmp = kernel2case(s1,s2,r1,r2,n1,n2,_alpha,_beta(i),_filter);
                    
                    abslogRatio[j - 1] = fabs(tmp);
                    
                    //  more favorable for equal mean
                    if(tmp > 0)
                    {
                        baseClus[j] = baseClus[j - 1];
                    }
                    else
                    {
                        // DE start a new cluster
                        baseClus[j] = baseClus[j - 1] + 1;
                    }
                    
                }
                // basebit and abslogRatio is ordered by mean expression
                
                auto tmpOrd = helper::sortIndexes<std::vector<Float>>(abslogRatio);
                
                auto baseBit = partition::mapToBit(baseClus);
                
                // use vector avoid value ref by eigen
                Eigen::VectorXd csOrd(K);
                
                Eigen::VectorXd rsOrd(K);
                
                auto sz = _clusinfo.size;
                
                auto szOrd = sz;
                
                
                for(size_t j = 0; j < K; j++)
                {
                    csOrd(j) = _sum(i,ord[j]);
                    
                    rsOrd(j) = _rsum(i,ord[j]);
                    
                    szOrd[j] = sz[ord[j]];
                }
                
                
                
                equalHandle(baseBit, abslogRatio, csOrd, rsOrd,
                            _alpha, _beta(i), _threshold, _filter, szOrd, _nequal, i);
                
                
                
                int localUC = 0;
                
                for(auto score:abslogRatio)
                {
                    if(score < _threshold)
                    {
                        localUC++;
                    }
                }
                
                //                while(abslogRatio[localUC] < _threshold)
                //                {
                //                    localUC++;
                //                }
                
                _guc.push_back(localUC);
                
                localUC = std::min(localUC , _uncertainty);
                
                if(localUC < 1)
                {
                    
                    auto newClusOr = baseClus;
                    
                    for(size_t iter = 0; iter < ord2.size(); iter++)
                    {
                        newClusOr[iter] = baseClus[ord2[iter]];
                    }
                    
                    auto newClusOrd = partition::reorder(newClusOr);
                    
                    auto sClus = partition::toString<std::vector<int>>(newClusOrd);
                    
                    if(dep.find(sClus) == dep.end())
                    {
                        dep.insert(sClus);
                        
                        _dep.push_back(newClusOrd);
                        
                        _pat.push_back(partition::toMatrix(newClusOrd));
                    }
                }
                
                auto pBit = partition::genBit(localUC);
                
                // get promising DE pattern
                for(auto x:pBit)
                {
                    
                    auto newBit = baseBit;
                    
                    for(int t = 0; t < _uncertainty; t++)
                    {
                        auto tmpJ = tmpOrd[t];
                        
                        newBit[tmpJ] = x[t];
                    }
                    
                    auto newClus = partition::bitToPart(newBit);
                    
                    auto newClusOr = newClus;
                    
                    for(size_t iter = 0; iter < ord2.size(); iter++)
                    {
                        newClusOr[iter] = newClus[ord2[iter]];
                    }
                    
                    std::vector<int> newClusOrd = partition::reorder(newClusOr);
                    
                    auto sClus = partition::toString<std::vector<int>>(newClusOrd);
                    
                    if(dep.count(sClus) < 1)
                    {
                        dep.insert(sClus);
                        
                        _dep.push_back(newClusOrd);
                        
                        _pat.push_back(partition::toMatrix(newClusOrd));
                    }
                    
                }
                
            }
            
        }
        
        
        void equalHandle(std::vector<bool>& baseBit, std::vector<double>& logRatio, Eigen::VectorXd& cs, Eigen::VectorXd& rs,
                         Float& alpha, Float& beta, Float& th1, Float& th2, std::vector<int>& sz, int nequal, int id)
        {
            // corner case when lots of contiguous subtypes having big bayes factor supporting
            // equal expression while the head and tail be quite different
            int start = 0;
            int pos = 0;
            int counter = 0;
            
            bool allequal = true;
            
            while(pos < baseBit.size())
            {
                //equality and greater than threshold
                if(baseBit[pos] < 1 && logRatio[pos] > _threshold)
                {
                    pos++;
                    counter++;
                }
                else
                {
                    allequal = false;
                    if(counter > nequal)
                    {
                        
                        // hclust
                    	ALGO::hclust<Eigen::VectorXd>(cs,rs,logRatio,start,pos,
                                     alpha,beta,th1,th2,sz);
                    }
                    
                    // reset
                    pos++;
                    start = pos;
                    counter = 0;
                }
            }
            if(allequal)
            {
                if(counter > nequal)
                {
                    ALGO::hclust<Eigen::VectorXd>(cs,rs,logRatio,start,pos,
                                      alpha,beta,th1,th2,sz);
                }
            }
            
            
        }
        
       
        
        
        void shrinkage()
        {
            std::vector<std::vector<int>> newDep;
            
            std::vector<COUNTS> newPat;
            
            
            for(size_t i = 0; i < _dep.size(); i++)
            {
                if(_p(i) > _shrinkthre)
                {
                    newDep.push_back(_dep[i]);
                    
                    newPat.push_back(_pat[i]);
                }
            }
            
            Eigen::VectorXd newP;
            
            newP.resize(newDep.size());
            
            size_t start = 0;
            
            for(size_t i = 0; i < _dep.size(); i++)
            {
                if(_p(i) > _shrinkthre)
                {
                    newP(start) = _p(i);
                    
                    start++;
                }
            }
            
            newP = newP / newP.sum();
            
            std::swap(_dep,newDep);
            
            std::swap(_pat,newPat);
            
            _p = newP;
        }
        
        Float kernel2case(Float& s1, Float& s2, Float& r1, Float& r2,
                          int& n1, int& n2, Float& alpha, Float& beta, Float& filter)
        {
            // if too small mean, assume they are the same
            if(s1 / n1 < filter && s2 / n2 <filter )
            {
                return INT_MAX;
            }
            
            Float res = lbeta(alpha + r1 + r2, beta + s1 + s2) + lbeta(alpha, beta) - lbeta(alpha + r1, beta + s1) - lbeta(alpha + r2, beta + s2);
            
            return res;
        }
        
        void kernel()
        {
            // adjust dim of kernel matrix
            _kernel.resize(_sum.rows(),_pat.size());
            
            for(size_t i = 0; i < _pat.size(); i++)
            {
                COUNTS _csum = _sum * _pat[i];
                
                COUNTS rsum = _rsum * _pat[i];
    
                COUNTS A = (rsum.array() + _alpha).matrix();
                
                COUNTS B = _csum.colwise() + _beta;
                
                COUNTS res = lbeta(A,B);
                
                res = (res.array() - boost::math::lgamma(_alpha)).matrix();
                
                res =  res.colwise() - (_beta.unaryExpr<Float(*)(Float)>(&boost::math::lgamma) + (_alpha + _beta.array()).matrix().unaryExpr<Float(*)(Float)>(&boost::math::lgamma));
    
                _kernel.col(i) = res.rowwise().sum();
                
            }
            
        }
        
        void gradientAscent()
        {
            size_t G = _sum.rows();

            size_t npat = _pat.size();

            COUNTS alpDRV(G,npat);

            COUNTS betaDRV(G,npat);

            for(size_t i = 0; i < npat; i++)
            {
                COUNTS _csum = _sum * _pat[i];
                
                COUNTS rsum = _rsum * _pat[i];

                COUNTS A = (rsum.array() + _alpha).matrix();

                COUNTS B = _csum.colwise() + _beta;

                COUNTS C = (A + B).unaryExpr<Float(*)(Float)>(&(boost::math::digamma));

                Eigen::VectorXd D = (_alpha + _beta.array()).matrix().unaryExpr<Float(*)(Float)>(&(boost::math::digamma));

                COUNTS resAlpha = ((A.unaryExpr<Float(*)(Float)>(&(boost::math::digamma)) - C).array() - boost::math::digamma(_alpha)).matrix();

                resAlpha = resAlpha.colwise() + D;

                COUNTS resBeta = B.unaryExpr<Float(*)(Float)>(&(boost::math::digamma)) - C;

                resBeta = resBeta.colwise() - (_beta.unaryExpr<Float(*)(Float)>(&(boost::math::digamma)) + D);

                alpDRV.col(i) = resAlpha.rowwise().sum();

                betaDRV.col(i) = resBeta.rowwise().sum();
            }

            auto tmp1 = _alpha + _lrate[0] * (alpDRV * _p).sum();
            
            //auto tmp2 = _beta + _lrate[1] * (betaDRV.array() * _post.array()).matrix().rowwise().sum();
            
            Eigen::VectorXd tmp2;
            
            if(_ng < G)
            {
                // has to make extra copy
                Eigen::VectorXd bt = betaDRV * _p;
                
                for(size_t localIter = 0; localIter < _ng; localIter++)
                {
                    Float localGrad = 0;
                    
                    for(auto localPos:_isoPos[localIter])
                    {
                        localGrad += bt(localPos);
                    }
                    
                    for(auto localPos:_isoPos[localIter])
                    {
                        bt(localPos) = localGrad;
                    }
                }
                
                tmp2 = _beta + _lrate[1] * bt;
            }else
            {
                tmp2 = _beta + _lrate[1] * (betaDRV * _p);
            }
            
            
            // check validity
            if(tmp1 > 0)
                _alpha = tmp1;

            for(size_t i = 0; i < G; i++)
            {
                if(tmp2(i) > 0)
                {
                    _beta(i) = tmp2(i);
                }
            }
        }
        
        void posterior()
        {
            assert(fabs(_p.sum() - 1) < 0.0001);
            
            Eigen::VectorXd M = _kernel.rowwise().maxCoeff();
            
            auto rmMax = _kernel.colwise() - M;

            _post = rmMax.unaryExpr<Float(*)(Float)>(& exp);
            
            Eigen::VectorXd total = _post * _p;
        
            total = (1 / total.array()).matrix();
           
            //outer product of total and p
            COUNTS div = total * _p.transpose();
            
            _post = (_post.array() * div.array()).matrix();
            
        }
        
        void updateP()
        {
            _p = _post.colwise().sum() / _post.sum();
        }
        
//        void updateMDE()
//        {
//            _mde.fill(0);
//
//            size_t K = _sum.cols();
//
//            for(size_t i = 0; i < K; i++)
//            {
//                for(size_t j = i; j < K; j++)
//                {
//                    for(size_t p = 0; p < _dep.size(); p++)
//                    {
//                        if(_dep[p][i] == _dep[p][j])
//                        {
//                            _mde(i,j) += _p[p];
//                        }
//                    }
//
//                }
//            }
//        }
        
        
        
    private:
        // get type of row in eigen matrix
        //typedef decltype(_mean.row(0)) ROW;
        
        int _nequal;
        
        //  var per group
        COUNTS _var;
        
        // q per group
        COUNTS _q;
        
        // pool var
        COUNTS _poolVar;
        
        // hyper parameter r can be estimated by MM
        COUNTS _r,_rsum;
        
        // alpha for the beta prior shared by genome, estimated by EM
        Float _alpha;
        
        // beta for the beta prior specified by each gene, estimated by EM
        Eigen::VectorXd _beta;
        
        // label for isoform.
        std::vector<int> _isoLabel;
        
        // number of isoforms
        size_t _ng;
        
        // position for each label
        std::vector<std::vector<int> > _isoPos;
        
        // step size for update hyper parameters
        std::vector<Float> _lrate;
        
        // prop of each nonzero pattern
        Eigen::VectorXd _p;
        
        // upper bound of unsure "<" and "=", controlling number of patterns DE
        int _uncertainty;
        
        // at least one group mean should be no smaller than this value to do DE comparison
        Float _filter;
        
        // positve threshold to decide how many uncertain patterns
        Float _threshold;
        
        // positive threshold for iterative shrinkage size of de patterns
        Float _shrinkthre;
        
        // gene level uncertainty
        std::vector<size_t> _guc;
        
        // matrix for DE pattern
        std::vector<COUNTS> _pat;
        
        // vector for DE pattern
        std::vector<std::vector<int>> _dep;
        
        // kernel matrix
        COUNTS _kernel;
        
        // marginal for two subtypes being ED or DD
//        COUNTS _mde;
        
        // posterior prob
        COUNTS _post;
        
    };
    
};
