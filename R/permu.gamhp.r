#' Permutation Test of Hierarchical Partitioning for GAM Analysis

#' @param  mod gam model generated by mgcv::gam()
#' @param  type The type of total explained variation, either "dev" or "adjR2", in which "dev" is deviance explained and "adjR2" is adjusted R-square, the default is "adjR2". 

#' @param  permutations An integer; Number of permutations for computing p value of individual contribution for the randomized dataset.

#' @details This function is a permutation test of hierarchical partitioning for gam analysis. It returns a matrix of I values (the individual contribution towards total explained variation) for all values from permutations randomizations. For each permutation, the values in each variable (i.e each column of iv) are randomized independently, and gam.hp is run on the randomized iv. As well as the randomized I matrix, the function returns a summary table listing the observed I values, the p value of I for the randomized dataset.

#' @return a data.frame containing a summary table listing the observed individual contribution, the p value of individual contribution for the randomized dataset

#' @author {Jiangshan Lai} \email{lai@njfu.edu.cn}


#'@export
#'@examples
#'library(mgcv)
#'mod1 <- gam(Sepal.Length ~ s(Petal.Length) + s(Petal.Width) + Sepal.Width,data = iris)
#'permu.gamhp(mod1,type="dev",permutations=10)


permu.gamhp <- function(mod=NULL,type='dev',permutations=10){
  #cat("\nPlease wait: running", permutations-1, "permutations \n") 
  
  #print(mod)
  #if (mod$method=='GCV') {method_name='GCV.Cp'} else if (mod$method=='GACV') {method_name='GACV.Cp'} else{method_name=mod$method}
  #print(paste('method name is ',method_name,sep=''))
  y_name=as.character(mod$formula)[2]
  #print(y_name)
  x_name=names(mod$var.summary)
  d_y=as.data.frame(mod$y)
  colnames(d_y) <- y_name
  #print(d_y)
  
  iv <- as.data.frame(mod$model[,x_name])
  #print(iv)
  obs <- gam.hp(mod=mod,type=type,commonality = FALSE)
  r2q_p <- obs$hierarchical.partitioning[,3]
  n <- dim(iv)[1]
  nvar <- dim(iv)[2]
  for(i in seq(1,permutations-1,1)){
    #print(i)
    newiv<-iv
    for(j in 1:nvar){
      perms <- sample(1:n,n)
      newiv[,j] <-iv[,j][perms]}
    #print(head(newiv))

    #print(method_name)
    data_p <- as.data.frame(cbind(d_y,newiv))
    #print(head(data_p))
    # print(mod$family)
    # print(method_name)
    mod_p<- stats::update(object = mod, formula. = .~., data = data_p)
    #mod_p <- gam(mod$formula,data=data_p,family=mod$family,method=method_name,select=FALSE)
    # print(mod_p)
    # summary(mod_p)
    # print(class(mod_p))
    
    
    
    simu_p=gam.hp(mod=mod_p,type=type,commonality = FALSE)
    # print('simu')
    # print(simu_p)
    r2q_p=cbind(r2q_p,simu_p$hierarchical.partitioning[,3])
  }
  
  Signi <- function(x){pval=round(1-ecdf(x)(x[1])+1/(permutations+1),nchar(permutations))
  if (pval <= 0.001) {
    return(noquote(paste(pval,"***",sep=" ")))
  }
  if (pval <= 0.01) {
    return(noquote(paste(pval," **",sep=" ")))
  }        
  if (pval <= 0.05) {
    return(noquote(paste(pval,"  *",sep=" ")))
  }
  else {return(noquote(paste(pval,"   ",sep=" ")))
  }
  }
  
  p.R2_p=apply(r2q_p,1,Signi)
  result=data.frame(obs$hierarchical.partitioning[,3],Pr=p.R2_p)
  colnames(result)<-c("Individual","Pr(>I)")
  return(result)
}
