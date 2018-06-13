#' Intersect samples between datasets.
#'
#' In order to make all datasets comparable, we have to make sure they describe
#' the same set of samples. This function takes a list of datasets (i.e. data matrices),
#' takes the intersect of all rownames, and returns a list of datasets with only
#' those samples.
#'
#' @param data A list of data matrices. The data matrices need to have rownames.
#' @return A list with of data matrices, all with the same set of samples.
#' @examples
#' set.seed(2)
#' n = 100
#' p = 100
#' x1 = matrix(rnorm(n*p), n, p)
#' x2 = matrix(rnorm(n*p), n, p)
#' rownames(x1) = rownames(x2) = paste0("X",1:n)
#' data = list(x1=x1[1:90,], x2=x2[10:100,])
#' data = intersect.samples(data)
#' @export
intersect.samples = function(data) {
  if(!is.list(data)) {
    stop("The parameter data should be a list.")
  }
  if(any(sapply(data, function(x){is.null(rownames(x))}))) {
    stop("At least one of the datasets does not have rownames.")
  }
  ind = Reduce(intersect, lapply(data, rownames))
  for(i in 1:length(data)) {
    data[[i]] = data[[i]][ind,,drop=F]
  }
  return(data)
}

#' Inner product similarity.
#'
#' Computes the inner product between x and y.
#'
#' @param x A vector of numbers.
#' @param y A vector of numbers.
#' @return The inner product similarity between x and y.
#' @examples
#' set.seed(2)
#' n = 100
#' x = rnorm(n)
#' y = rnorm(n)
#' inner.product(x, y)
#' @export
inner.product = function(x,y) {
  if(length(x)!=length(y)) {
    stop("The vectors x and y should be of the same length.")
  }
  if(!is.numeric(x) | !is.numeric(y)) {
    stop("The vectors x and y should be numeric.")
  }
  return(t(x)%*%y)
}

#' Jaccard similarity.
#'
#' Computes the Jaccard similarity between x and y. When both x and
#' y only contain zeroes, the Jaccard similarity it not defined. This
#' function returns zero for that specific case.
#'
#' @param x A vector of zeroes and ones.
#' @param y A vector of zeroes and ones.
#' @return The Jaccard similarity between x and y.
#' @examples
#' set.seed(2)
#' n = 100
#' x = rbinom(n, 1, 0.5)
#' y = rbinom(n, 1, 0.5)
#' jaccard(x, y)
#' @export
jaccard = function(x,y) {
  if(length(x)!=length(y)) {
    stop("The vectors x and y should be of the same length.")
  }
  if(!all(unique(c(x,y)) %in% c(0,1))) {
    stop("The vectors x and y should contain only 0s and 1s.")
  }
  a = sum(x==0 & y==0)
  b = sum(x==1 & y==0)
  c = sum(x==0 & y==1)
  d = sum(x==1 & y==1)

  if(all(c(b,c,d)==0)) {
    return(0)
  } else {
    return(d/(b+c+d))
  }
}

#' Process a custom configuration matrix.
#'
#' This function can be used to process a custom-made configuration matrix (i.e. similarity matrix) for use with the RV coefficient.
#' The function can perform two tasks: centering and preparation for the modified RV coefficient, both of which we will briefly explain here.
#'
#' The RV coefficient often results in values very close to one when both datasets are not centered around zero, even for orthogonal data.
#' For inner product similarity and Jaccard similarity, we recommend using centering. However, for some other similarity measures, centering
#' may not be beneficial (for example, because the measure itself is already centered, such as in the case of Pearson correlation). For more information on
#' centering of binary (and other non-continuous) data, for which we used kernel centering of the configuration matrix, we refer to our manuscript: Aben et al., 2018, doi.org/10.1101/293993.
#'
#' The modified RV coefficient was proposed for high-dimensional data, as the regular RV coefficient would result in values close to one even for
#' orthogonal data. We recommend always using the modified RV coefficient.
#'
#' @param S A configuration matrix.
#' @param center Should the configuration matrix be centered using kernel centering?
#' @param mod.rv Should the configuration matrix be prepared for the modified RV coefficient?
#' @return The processed configuration matrix.
#' @examples
#' set.seed(2)
#' n = 100
#' p = 100
#' x = matrix(rnorm(n*p)+10, n, p)
#' S = x%*%t(x)
#' S_dash = process.custom.config.matrix(S, center=TRUE, mod.rv=TRUE)
#' @export
process.custom.config.matrix = function(S, center=TRUE, mod.rv=TRUE) {
  if(!is.matrix(S) | !is.numeric(S)) {
    stop("S should be a numeric matrix.")
  }
  if(nrow(S) != ncol(S)) {
    stop("S should be square matrix.")
  }
  if(!is.logical(center)) {
    stop("center should be a boolean.")
  }
  if(!is.logical(mod.rv)) {
    stop("mod.rv should be a boolean.")
  }

  n = nrow(S)
  if(center) {
    S = S - (1/n)*(S%*%rep(1,n))%*%rep(1,n) - (1/n)*rep(1,n)%*%(rep(1,n)%*%S) + (1/n^2)*sum(S)
  }
  if(mod.rv) {
    diag(S) = 0
  }
  return(S)
}

#' Computes a configuration matrix
#'
#' Given a data matrix, this function computes the configuration matrix for the corresponding dataset. You'll typically won't need to call
#' this function directly, but should use compute.config.matrices() instead, as it will make determining partial RV coefficients, p-values and confidence
#' intervals easier later on.
#'
#' @param x Data matrix.
#' @param similarity_fun A function pointer to the similarity function to be used (default=inner.product).
#' @param center A boolean indicating whether centering should be used (default=TRUE).
#' @param mod.rv A boolean indicating whether the modified RV coefficient should be used (default=TRUE).
#' @return A configuration matrix.
#' @examples
#' set.seed(2)
#' n = 100
#' p = 100
#' x1 = matrix(rnorm(n*p), n, p)
#' x2 = x1 + matrix(rnorm(n*p), n, p)
#' S1 = compute.config.matrix(x1)
#' S2 = compute.config.matrix(x1)
#' rv.coef(S1, S2)
#' @export
compute.config.matrix = function(x, similarity_fun=inner.product, center=TRUE, mod.rv=TRUE) {
  if(!is.matrix(x)) {
    stop("x should be a matrix.")
  }
  if(!is.function(similarity_fun)) {
    stop("similarity_fun should be a function pointer.")
  }
  if(!is.logical(center)) {
    stop("center should be a boolean.")
  }
  if(!is.logical(mod.rv)) {
    stop("mod.rv should be a boolean.")
  }

  n = nrow(x)
  S = matrix(NA, n, n)
  for(i in 1:n) {
    for(j in i:n) {
      S[i,j] = similarity_fun(x[i,], x[j,])
    }
  }
  S = as.matrix(Matrix::triu(S)) + Matrix::t(Matrix::triu(S, k=1))
  S = as.matrix(S)
  S = process.custom.config.matrix(S, center, mod.rv)
  return(S)
}

#' Compute configuration matrices
#'
#' Given a list of n data matrices (corresponding to n datasets), this function computes the configuration matrix for each of these
#' configuration matrices. By default inner product similarity is used, but other similarity (such as Jaccard similarity for binary data)
#' can also be used (see the vignette 'A quick introduction to iTOP' for more information). In addition, the configuration matrices can be centered and prepared for use with
#' the modified RV coefficient, both of which we will briefly explain here.
#'
#' The RV coefficient often results in values very close to one when both datasets are not centered around zero, even for orthogonal data.
#' For inner product similarity and Jaccard similarity, we recommend using centering. However, for some other similarity measures, centering
#' may not be beneficial (for example, because the measure itself is already centered, such as in the case of Pearson correlation). For more information on
#' centering of binary (and other non-continuous) data, for which we used kernel centering of the configuration matrix, we refer to our manuscript: Aben et al., 2018, doi.org/10.1101/293993.
#'
#' The modified RV coefficient was proposed for high-dimensional data, as the regular RV coefficient would result in values close to one even for
#' orthogonal data. We recommend always using the modified RV coefficient.
#'
#' @param data List of datasets.
#' @param similarity_fun Either a function pointer to the similarity function to be used for all datasets; or a list of function pointers,
#' if different similarity functions need to be used for different datasets (default=inner.product).
#' @param center Either a boolean indicating whether centering should be used for all datasets; or a list of booleans,
#' if centering should be used for some datasets but not all of them (default=TRUE).
#' @param mod.rv Either a boolean indicating whether the modified RV coefficient should be used for all datasets; or a list of booleans,
#' if the modified RV should be used for some datasets but not all of them (default=TRUE).
#' @return A list of n configuration matrices, where n is the number of datasets.
#' @examples
#' set.seed(2)
#' n = 100
#' p = 100
#' x1 = matrix(rnorm(n*p), n, p)
#' x2 = x1 + matrix(rnorm(n*p), n, p)
#' x3 = x2 + matrix(rnorm(n*p), n, p)
#' data = list(x1=x1, x2=x2, x3=x3)
#' config_matrices = compute.config.matrices(data)
#' cors = rv.cor.matrix(config_matrices)
#' @export
compute.config.matrices = function(data, similarity_fun=inner.product, center=TRUE, mod.rv=TRUE) {
  if(!is.list(data)) {
    stop("data should be a list.")
  }
  if(!length(similarity_fun) %in% c(1,length(data))) {
    stop("The length of similarity_fun should be equal to either 1 or to the number of datasets.")
  }
  if(!length(similarity_fun) %in% c(1,length(data))) {
    stop("The length of similarity_fun should be equal to either 1 or to the number of datasets.")
  }
  if(!length(similarity_fun) %in% c(1,length(data))) {
    stop("The length of similarity_fun should be equal to either 1 or to the number of datasets.")
  }
  if(!is.logical(center)) {
    stop("center should be a boolean.")
  }
  if(!is.logical(mod.rv)) {
    stop("mod.rv should be a boolean.")
  }

  if(length(similarity_fun)==1) {
    similarity_fun = replicate(length(data), similarity_fun)
  }
  if(length(center)==1) {
    center = rep(center,length(data))
  }
  if(length(mod.rv)==1) {
    mod.rv = rep(mod.rv,length(data))
  }

  config_matrices = list()
  for(i in 1:length(data)) {
    config_matrices[[i]] = compute.config.matrix(data[[i]], similarity_fun[[i]], center=center[i], mod.rv=mod.rv[i])
  }
  names(config_matrices) = names(data)
  return(config_matrices)
}

#' Computes the RV coefficient
#'
#' Computes the RV coefficient between dataset 1 and dataset 2. You'll typically won't need to call this function directly,
#' but should use rv.cor.matrix() instead, as it will make determining partial RV coefficients, p-values and confidence
#' intervals easier later on.
#'
#' @param S1 Configuration matrix corresponding to dataset 1
#' @param S2 Configuration matrix corresponding to dataset 2
#' @return The RV coefficient between dataset 1 and dataset 2
#' @examples
#' set.seed(2)
#' n = 100
#' p = 100
#' x1 = matrix(rnorm(n*p), n, p)
#' x2 = x1 + matrix(rnorm(n*p), n, p)
#' S1 = compute.config.matrix(x1)
#' S2 = compute.config.matrix(x1)
#' rv.coef(S1, S2)
#' @export
rv.coef = function(S1, S2) {
  if(!is.matrix(S1) | !is.numeric(S1)) {
    stop("S1 should be a numeric matrix.")
  }
  if(nrow(S1) != ncol(S1)) {
    stop("S1 should be square matrix.")
  }
  if(!is.matrix(S2) | !is.numeric(S2)) {
    stop("S2 should be a numeric matrix.")
  }
  if(nrow(S2) != ncol(S2)) {
    stop("S2 should be square matrix.")
  }
  if(length(S1) != length(S2)) {
    stop("S1 and S2 should be of the same size.")
  }

  rv = c(S1)%*%c(S2) / sqrt(c(S1)%*%c(S1) * c(S2)%*%c(S2))
  return(rv)
}

#' A correlation matrix of RV coefficients
#'
#' Given a list of n configuration matrices (corresponding to n datasets), this function computes an n x n matrix of pairwise RV coefficients.
#'
#' @param config_matrices The result from compute.config.matrices().
#' @return An n x n matrix of pairwise RV coefficients, where n is the number of datasets.
#' @examples
#' set.seed(2)
#' n = 100
#' p = 100
#' x1 = matrix(rnorm(n*p), n, p)
#' x2 = x1 + matrix(rnorm(n*p), n, p)
#' x3 = x2 + matrix(rnorm(n*p), n, p)
#' data = list(x1=x1, x2=x2, x3=x3)
#' config_matrices = compute.config.matrices(data)
#' cors = rv.cor.matrix(config_matrices)
#' @export
rv.cor.matrix = function(config_matrices) {
  if(!is.list(config_matrices)) {
    stop("config_matrices should be a list.")
  }
  n = length(config_matrices)
  cors = matrix(NA,n,n)
  colnames(cors) = names(config_matrices)
  rownames(cors) = names(config_matrices)
  for(i in 1:n) {
    for(j in i:n) {
      if(i==j) {
        cors[i,j] = 1
      } else {
        cors[i,j] = rv.coef(config_matrices[[i]], config_matrices[[j]])
      }
    }
  }
  cors = Matrix::triu(cors) + Matrix::t(Matrix::triu(cors, k=1))
  cors = as.matrix(cors)
  return(cors)
}

#' Performing a permutation
#'
#' Helper function for run.permutations(). It's unlikely you'll ever need to run this function directly.
#'
#' @param config_matrices The result from compute.config.matrices().
#' @return An n x n matrix of RV coefficients for the permutated data, where n is the number of datasets.
permute.config.matrices = function(config_matrices) {
  if(!is.list(config_matrices)) {
    stop("config_matrices should be a list.")
  }
  data_perm = config_matrices
  for(i in 1:length(data_perm)) {
    n = nrow(data_perm[[i]])
    ind = sample(1:n, n, replace=F)
    data_perm[[i]] = data_perm[[i]][ind,ind]
  }
  return(data_perm)
}

#' Permutations for significance testing
#'
#' Performs a permutations for significance testing. The result from this function can be used with
#' rv.pval() to determine a p-value. By decoupling this into two functions,
#' you don't have to redo the permutations for every p-value, hence increasing the runtime speed.
#'
#' @param config_matrices The result from compute.config.matrices().
#' @param nperm The number of permutations to perform (default=1000).
#' @return An n x n x nperms array of RV coefficients for the permutated data, where n is the number of datasets.
#' @examples
#' set.seed(2)
#' n = 100
#' p = 100
#' x1 = matrix(rnorm(n*p), n, p)
#' x2 = x1 + matrix(rnorm(n*p), n, p)
#' x3 = x2 + matrix(rnorm(n*p), n, p)
#' data = list(x1=x1, x2=x2, x3=x3)
#' config_matrices = compute.config.matrices(data)
#' cors = rv.cor.matrix(config_matrices)
#' cors_perm = run.permutations(config_matrices, nperm=1000)
#' rv.pval(cors, cors_perm, "x1", "x3", "x2")
#' @export
run.permutations = function(config_matrices, nperm=1000) {
  if(!is.list(config_matrices)) {
    stop("config_matrices should be a list.")
  }
  if(!is.numeric(nperm) | nperm%%1!=0 | nperm <0) {
    stop("nperm should be a positive integer.")
  }

  n = length(config_matrices)
  cors_perm = array(NA, dim=c(n,n,nperm))
  for(i in 1:nperm) {
    data_perm = permute.config.matrices(config_matrices)
    cors_perm[,,i] = rv.cor.matrix(data_perm)
  }
  colnames(cors_perm) = names(config_matrices)
  rownames(cors_perm) = names(config_matrices)
  return(cors_perm)
}

#' Performing a single bootstrap
#'
#' Helper function for run.bootstraps(). It's unlikely you'll ever need to run this function directly.
#'
#' @param config_matrices The result from compute.config.matrices().
#' @return An n x n matrix of RV coefficients for the bootstrapped data, where n is the number of datasets.
bootstrap.config.matrices = function(config_matrices) {
  data_boot = config_matrices
  n = nrow(data_boot[[1]])
  ind = sample(1:n, n, replace=TRUE)
  for(i in 1:length(data_boot)) {
    data_boot[[i]] = data_boot[[i]][ind,ind]
  }
  return(data_boot)
}

#' Bootstrapping procedure
#'
#' Performs a bootstrapping procedure. The result from this function can be used with
#' rv.conf.interval() to determine confidence intervals. By decoupling this into two functions,
#' you don't have to redo the bootstrapping for every confidence interval, hence increasing the runtime speed.
#'
#' @param config_matrices The result from compute.config.matrices().
#' @param nboots The number of bootstraps to perform (default=1000).
#' @return An n x n x nboots array of RV coefficients for the bootstrapped data, where n is the number of datasets.
#' @examples
#' set.seed(2)
#' n = 100
#' p = 100
#' x1 = matrix(rnorm(n*p), n, p)
#' x2 = x1 + matrix(rnorm(n*p), n, p)
#' x3 = x2 + matrix(rnorm(n*p), n, p)
#' data = list(x1=x1, x2=x2, x3=x3)
#' config_matrices = compute.config.matrices(data)
#' cors_boot = run.bootstraps(config_matrices, nboots=1000)
#' rv.conf.interval(cors_boot, "x1", "x3", "x2")
#' @export
run.bootstraps = function(config_matrices, nboots=1000) {
  if(!is.list(config_matrices)) {
    stop("config_matrices should be a list.")
  }
  if(!is.numeric(nboots) | nboots%%1!=0 | nboots <0) {
    stop("nboots should be a positive integer.")
  }
  n = length(config_matrices)
  cors_boot = array(NA, dim=c(n,n,nboots))
  for(i in 1:nboots) {
    data_boot = bootstrap.config.matrices(config_matrices)
    cors_boot[,,i] = rv.cor.matrix(data_boot)
  }
  colnames(cors_boot) = names(config_matrices)
  rownames(cors_boot) = names(config_matrices)
  return(cors_boot)
}

#' Determining a (partial) RV coefficient
#'
#' Determines the RV coefficient RV(a, b) or the partial RV coefficient
#' RV(a, b | set).
#'
#' @param cors The result from rv.cor.matrix().
#' @param a Either an index or a string to identify dataset a.
#' @param b Either an index or a string to identify dataset b.
#' @param set Optional parameter to define the datasets that need to be partialized for.
#' If set consists of one dataset, then provide an index or a string to identify set.
#' If set consists of multiple datasets, then provide a vector of indices or a vector of strings.
#' @return The (partial) RV coefficient.
#' @examples
#' set.seed(2)
#' n = 100
#' p = 100
#' x1 = matrix(rnorm(n*p), n, p)
#' x2 = x1 + matrix(rnorm(n*p), n, p)
#' x3 = x2 + matrix(rnorm(n*p), n, p)
#' data = list(x1=x1, x2=x2, x3=x3)
#' config_matrices = compute.config.matrices(data)
#' cors = rv.cor.matrix(config_matrices)
#' rv.pcor(cors, "x1", "x3", "x2")
#' @export
rv.pcor = function(cors, a, b, set=NULL) {
  if(!(is.character(a) | is.numeric(a)) | length(a)!=1) {
    stop("a should be either a numeric or a string.")
  }
  if(!(is.character(b) | is.numeric(b)) | length(b)!=1) {
    stop("b should be either a numeric or a string.")
  }
  if(!(is.character(set) | is.numeric(set) | is.null(set))) {
    stop("set should be either a numeric, a string, a vector of strings of numerics, or NULL.")
  }
  if(!is.matrix(cors) | nrow(cors)!=ncol(cors)) {
    stop("cors_boot should be an n x n matrix, with n the number of datasets.")
  }

  ind = c(a, b, set)
  pcors = corpcor::cor2pcor(cors[ind,ind])
  colnames(pcors) = rownames(pcors) = rownames(cors[ind,ind])
  return(pcors[1,2])
}

#' Wrapper function to determine significance in the PC algorithm
#'
#' This function is a wrapper function around rv.pval(), such that it
#' can easily be used with pc() from the pcalg package. If you have trouble installing the pcalg package, have a look at our vignette 'A quick start to iTOP'.
#'
#' @param a Either an index or a string to identify dataset a.
#' @param b Either an index or a string to identify dataset b.
#' @param set Datasets that need to be partialized for. Set to NULL if there are none (i.e. if you're computing a regular, non-partial RV).
#' If set consists of one dataset, then provide an index or a string to identify set.
#' If set consists of multiple datasets, then provide a vector of indices or a vector of strings.
#' @param suffStat A named list with two items: cors, which is the result from rv.cor.matrix(); and cors_perm, which is the result from run.permutations().
#' @return The p-value.
#' @examples
#' set.seed(2)
#' n = 100
#' p = 100
#' x1 = matrix(rnorm(n*p), n, p)
#' x2 = x1 + matrix(rnorm(n*p), n, p)
#' x3 = x2 + matrix(rnorm(n*p), n, p)
#' data = list(x1=x1, x2=x2, x3=x3)
#' config_matrices = compute.config.matrices(data)
#' cors = rv.cor.matrix(config_matrices)
#' cors_perm = run.permutations(config_matrices, nperm=1000)
#'
#' \dontrun{
#' library(pcalg)
#' suffStat = list(cors=cors, cors_perm=cors_perm)
#' pc.fit = pc(suffStat=suffStat, indepTest=rv.link.significance, labels=names(data),
#'             alpha=0.05, conservative=TRUE, solve.confl=TRUE)
#' plot(pc.fit, main="")}
#' @export
rv.link.significance = function(a, b, set, suffStat) {
  p = rv.pval(suffStat$cors, suffStat$cors_perm, a, b, set)
  return(p)
}

#' Determining a p-value the (partial) RV coefficient
#'
#' This function uses a permutation test to determine a p-value
#' for the RV coefficient RV(a, b) or the partial RV coefficient
#' RV(a, b | set).
#'
#' @param cors The result from rv.cor.matrix().
#' @param cors_perm The result from run.permutations().
#' @param a Either an index or a string to identify dataset a.
#' @param b Either an index or a string to identify dataset b.
#' @param set Optional parameter to define the datasets that need to be partialized for.
#' If set consists of one dataset, then provide an index or a string to identify set.
#' If set consists of multiple datasets, then provide a vector of indices or a vector of strings.
#' @return The p-value.
#' @examples
#' set.seed(2)
#' n = 100
#' p = 100
#' x1 = matrix(rnorm(n*p), n, p)
#' x2 = x1 + matrix(rnorm(n*p), n, p)
#' x3 = x2 + matrix(rnorm(n*p), n, p)
#' data = list(x1=x1, x2=x2, x3=x3)
#' config_matrices = compute.config.matrices(data)
#' cors = rv.cor.matrix(config_matrices)
#' cors_perm = run.permutations(config_matrices, nperm=1000)
#' rv.pval(cors, cors_perm, "x1", "x3", "x2")
#' @export
rv.pval = function(cors, cors_perm, a, b, set=NULL) {
  if(!(is.character(a) | is.numeric(a)) | length(a)!=1) {
    stop("a should be either a numeric or a string.")
  }
  if(!(is.character(b) | is.numeric(b)) | length(b)!=1) {
    stop("b should be either a numeric or a string.")
  }
  if(!(is.character(set) | is.numeric(set) | is.null(set))) {
    stop("set should be either a numeric, a string, a vector of strings of numerics, or NULL.")
  }
  if(!is.array(cors_perm) | length(dim(cors_perm))!=3 | nrow(cors_perm)!=ncol(cors_perm)) {
    stop("cors_perm should be an n x n x nboots array, with n the number of datasets and nperm the number of permutations.")
  }
  if(!is.matrix(cors) | nrow(cors)!=ncol(cors)) {
    stop("cors_boot should be an n x n matrix, with n the number of datasets.")
  }
  if(nrow(cors)!=ncol(cors_perm)) {
    stop("cors and cors_boot should have the same number of rows and columns.")
  }

  nperm = dim(cors_perm)[3]
  pcor = rv.pcor(cors, a, b, set)
  pcor_perm = c()
  for(i in 1:nperm) {
    pcor_perm[i] = rv.pcor(cors_perm[,,i], a, b, set)
  }
  if(pcor>0) {
    p = mean(pcor < pcor_perm)
  } else {
    p = mean(pcor > pcor_perm)
  }
  return(p)
}

#' Determining a confidence interval for the (partial) RV coefficient
#'
#' This function uses a bootstrapping procedure to determine a confidence
#' interval for the RV coefficient RV(a, b) or the partial RV coefficient
#' RV(a, b | set).
#'
#' @param cors_boot The result from run.bootstraps().
#' @param a Either an index or a string to identify dataset a.
#' @param b Either an index or a string to identify dataset b.
#' @param set Optional parameter to define the datasets that need to be partialized for.
#' If set consists of one dataset, then provide an index or a string to identify set.
#' If set consists of multiple datasets, then provide a vector of indices or a vector of strings.
#' @param conf The size of the confidence interval (default=0.95).
#' @return The confidence interval.
#' @examples
#' set.seed(2)
#' n = 100
#' p = 100
#' x1 = matrix(rnorm(n*p), n, p)
#' x2 = x1 + matrix(rnorm(n*p), n, p)
#' x3 = x2 + matrix(rnorm(n*p), n, p)
#' data = list(x1=x1, x2=x2, x3=x3)
#' config_matrices = compute.config.matrices(data)
#' cors_boot = run.bootstraps(config_matrices, nboots=1000)
#' rv.conf.interval(cors_boot, "x1", "x3", "x2")
#' @export
rv.conf.interval = function(cors_boot, a, b, set=NULL, conf=.95) {
  if(!(is.character(a) | is.numeric(a)) | length(a)!=1) {
    stop("a should be either a numeric or a string.")
  }
  if(!(is.character(b) | is.numeric(b)) | length(b)!=1) {
    stop("b should be either a numeric or a string.")
  }
  if(!(is.character(set) | is.numeric(set) | is.null(set))) {
    stop("set should be either a numeric, a string, a vector of strings of numerics, or NULL.")
  }
  if(!is.numeric(conf) | conf<0 | conf>1) {
    stop("conf should be a number between 0 and 1.")
  }
  if(!is.array(cors_boot) | length(dim(cors_boot))!=3 | nrow(cors_boot)!=ncol(cors_boot)) {
    stop("cors_boot should be an n x n x nboots array, with n the number of datasets and nboots the number of bootstraps.")
  }

  nperm = dim(cors_boot)[3]
  pcor_boot = c()
  for(i in 1:nperm) {
    pcor_boot[i] = rv.pcor(cors_boot[,,i], a, b, set)
  }
  conf = (1+conf)/2
  conf_interval = stats::quantile(pcor_boot, probs=c(1-conf, conf))
  return(conf_interval)
}
