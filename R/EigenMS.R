# This software is released unde The GNU General Public License:
# http://www.gnu.org/licenses/gpl.html
 
# EigenMS normalization
# Ref: "Normalization of peak intensities in bottom-up MS-based proteomics using
#       singular value decomposition" Karpievitch YV, Taverner T, Adkins JN,
#       Callister SJ, Anderson GA, Smith RD, Dabney AR. Bioinformatics 2009
# and
# Ref: "Metabolomics data normalization with EigenMS"
#       Karpievitch YK, Nikolic SB, Wilson R, Sharman JE, Edwards LM
#       Submitted to PLoS ONE
#
# Here we allow multiple factors to be preserved in ANOVA
# model before identifying bais. This requires a bit of
# thought on the side of the researchers, only few factors
# should be 'preserved'. For example 'treatment group' is
# important to preserve but 'age' may or may not be important
# to preserve. Here we do not utilize peptides with 1+ grp
# missing completely, this is a separate problem
# addressed by Wang et al, 2011
#
# Written by Yuliya Karpievitch, Tom Taverner, and Shelley Herbrich
# email: yuliya.k@gmail.com
#
# Version has been split into 2 main functions:
# eig_norm1 finds significant bias trends, then the
# user can decide if he/she wants to use that number.
# For Metabolomics data we suggest using 20% of the
# number of samples examined as the number of bias
# trends to eliminated. By setting
# ints_eig1$h.c = ceil(.2 * num_samples)
# if the number of automatically identified bias trends
# are close to that number we suggest
# using the estimated number of bias trends.
#
# eig_norm2 function normalizes the data
# Note that rescaling has been abandoned in the latest
# version due to discussions abot the fact
# that systematic bias is what has been inadvertently
# added, thus we need to remove the bias
# but nto add any additional variation.
# Metabolomics in particular has a lot of variation
# that still remains that we concluded
# there is not need for rescaling.
#
# EigenMS estimates and preserves fixed effects.
# Contributions to using Random effects are welcome.

# HOWTO RUN:
# source('~/EigenMS/EigenMS.R')
# dd = # read in the data
# grps = # read in group information file
# logInts = # subset the dd matrix to only the portion that contains intensities
#           # replace 0's with NAs, do a log-transformation to
#           # approximate Normality.
# prot.info = cbind(data.frame(rownames(dd)), data.frame(rownames(dd)))
#             # peptideIDs, no PR IDs here, just duplicate column
#             # in case od protein IDs, those are not importnat for
#             # normalization and will be ignored.
# ints_eig1 = eig_norm1(m=logInts, treatment=grps, prot.info=prot.info)
# ints_norm = eig_norm2(rv=ints_eig1)
#
# eig_norm1 = function(m, treatment, prot.info, write_to_file='')
#  m - m x n matrix of log-transformed intensities, num peptides x num samples
#  treatment - either a single factor or a data.frame of factors
#  (actual R factors, else code will fail)
#  eg:  bl = gl(3,14,168) # 3 factors, repeat every 14 times,
#  total of 168 samples
#  it = gl(2,7,168)  # 2 factors, repeat every 7 samples, 168 total
#  Pi = gl(2,42,168) # 2 factors, repeat twice: 42 of PI+, PI-, PI+, PI-
#  grpFactors = data.frame(bl,it,Pi)# factors we would
#  like to preserve in EigenMS
#  prot.info - 2 column data frame with peptide and protein IDs in that order.
#  for metabolites both columns should contain metabolite IDs.
#  write_to_file -  if a string is passed in, 'complete' peptides will be
#  written to that file name
#  Some peptides could be eliminated due to too many missing values
#  (if not able to do ANOVA)

# SUPPLEMENTARY FUNCTION USED BY EIGENMS
# plot top 3 eigentrends with a line at 0.
# in cse of a single treatment group u and v matrices are switched..
plot.eigentrends = function(svdr, title1){
  v = svdr$v
  d = svdr$d
  ss = d^2
  Tk = signif(ss/sum(ss)* 100, 2)

  titles = paste("Trend ", seq(1,3), " (", Tk[seq(1,3)], "%)", sep = "")
  do.text = function(j) graphics::mtext(titles[j], cex=0.7, padj=-0.7, adj=1)
  range.y = range(as.numeric(v[,seq(1,3)]), na.rm=TRUE)

  toplot1_1 = as.numeric(v[,1])
  toplot1_2 = as.numeric(v[,2])
  toplot1_3 = as.numeric(v[,3])

  plot(seq(1,length(toplot1_1) ), toplot1_1, type='b', ann=FALSE, ylim=range.y)
  do.text(1)
  graphics::abline(h=0, lty=3)
  graphics::title(title1, cex.main = 1.2, font.main= 1, col.main= "purple",
                  ylab=NULL)
  plot(seq(1,length(toplot1_2) ), toplot1_2, type='b', ann=FALSE, ylim=range.y)
  do.text(2)
  graphics::abline(h=0, lty=3)
  plot(seq(1,length(toplot1_3) ), toplot1_3, type='b', ann=FALSE, ylim=range.y)
  do.text(3)
  graphics::abline(h=0, lty=3)
  return(Tk)
}

# plot any 3 eigentrends with a line at 0. Strting at the trend number passed in
# pos1 parameter provide the starting index for teh 3 trends
# in cse of a single treatment group u and v matrices are switched..
plot.eigentrends.start = function(svdr, title1, pos1=1){
  # No check for valid range of pos1 is performed!!!
  v = svdr$v
  d = svdr$d
  ss = d^2
  Tk = signif(ss/sum(ss)* 100, 2)
  titles = paste("Trend ", pos1:(pos1+3),
                 " (", Tk[pos1:(pos1+3)], "%)", sep = "")
  do.text = function(j) graphics::mtext(titles[j], cex=0.7, padj=-0.7, adj=1)
  range.y = range(as.numeric(v[,pos1:(pos1+3)]), na.rm=TRUE)

  toplot1_1 = as.numeric(v[,pos1])
  toplot1_2 = as.numeric(v[,(pos1+1)])
  toplot1_3 = as.numeric(v[,(pos1+2)])

  plot(seq(1,length(toplot1_1) ), toplot1_1, type='b', ann=FALSE, ylim=range.y)
  do.text(1)
  graphics::abline(h=0, lty=3)
  graphics::title(title1, cex.main = 1.2, font.main= 1, col.main= "purple",
                  ylab=NULL)
  plot(seq(1,length(toplot1_2) ), toplot1_2, type='b', ann=FALSE, ylim=range.y)
  do.text(2)
  graphics::abline(h=0, lty=3)
  plot(seq(1,length(toplot1_3) ), toplot1_3, type='b', ann=FALSE, ylim=range.y)
  do.text(3)
  graphics::abline(h=0, lty=3)
  return(Tk)
}



#' String linear model formula suitable
#'
#' Makes a string linear model formula suitable for the right hand side of the
#' equasion passed into lm()
#'
#' eig_norm1 and eig_norm2
#' Here we incorporate the model matrix from EigenMS
#' normalization to find the significant
#' trends in the matrix of residuals.
#'
#' @param eff treatment group ordering for all samples being anlysed.
#'            Single factor with 2+ teatment groups.
#'            Used to generate formula and contrasts for lm().
#' @param var_name string variable name to use in the formula
#' @return data structure with linea model formula and contrasts
#'  \describe{
#'   \item{lm.formula}{Lienar model formula suitable for right hand side of '
#'         ~' in lm(), ~ is not included int eh formula}
#'   \item{lm.params}{contrasts for lm(), here sum-to-zero constraint only}
#' }
#' @examples
#' grps = as.factor(c('CG', 'CG', 'CG', 'mCG', 'mCG', 'mCG'))
#' makeLMFormula(grps, 'TREATS')
#' @export
makeLMFormula = function(eff, var_name='') {
  # eff - effects used in contrasts
  # var_name - for singe factor use var-name that is passed in as variable names
  if(is.factor(eff))
  {
    ndims = 1
    cols1 = var_name # ftemp in EigenMS
  }
  else
  {
    ndims = dim(eff)[2]
    cols1 = colnames(eff)
  }
  lhs = cols1[1]
  lm.fm = NULL
  # check if can have a list if only have 1 factor...

  params = paste('contrasts=list(', cols1[1], '=contr.sum', sep=)

  if (ndims > 1) { # removed ndims[2] here, now ndims holds only 1 dimention...
    for (ii in 2:length(cols1))
    {
      lhs = paste(lhs, "+", cols1[ii])  # bl="contr.sum",
	    params = paste(params, ',', cols1[ii], '=contr.sum', sep='')
	  }
  }
  params = paste(params,")")
  lm.formula = stats::as.formula(paste('~', lhs))
  lm.fm$lm.formula = lm.formula
  lm.fm$lm.params = params
  return(lm.fm)
}



# EigenMS normalization
#
# Ref: "Normalization of peak intensities in bottom-up MS-based proteomics using
#       singular value decomposition" Karpievitch YV, Taverner T, Adkins JN,
#       Callister SJ, Anderson GA, Smith RD, Dabney AR. Bioinformatics 2009
# and
# Ref: "Metabolomics data normalization with EigenMS"
#       Karpievitch YK, Nikolic SB, Wilson R, Sharman JE, Edwards LM
#       Submitted to PLoS ONE
#
# Here we allow multiple facotrs to be preserved in ANOVA model before
# identifying bais.
# This requires a bit of thought on the side of the researchers, only few
# factors should be 'preserved'.
# For example 'treatment group' is important to preserve but 'age' may or
# may not be important to preserve.
# Here we do not utilize peptides with 1+ grp missing completely,
# this is a separate problem  addressed by Wang et al, 2011.
#
# Written by Yuliya Karpievitch, Tom Taverner, and Shelley Herbrich
# email: yuliya.k@gmail.com
#
# Version has been split into 2 main functions: eig_norm1 finds significant
# bias trends, then the user can decide if he/she wants to use that number.
# For Metabolomics data we suggest useing 20% of the number of samples examined
# as the number of bias trends to eliminated. By setting
# ints_eig1$h.c = ceil(.2 * num_samples)
# if the number of automatically identified bias trends are close to that
# number we suggest
# using the estimated number of bias trends.

#' Identify bias trends
#'
#' First portion of EigenMS: Identify eigentrends attributable to bias, allow
#' the user to adjust the number (with causion! if desired) before normalizing
#' with eig_norm2.
#' Ref: "Normalization of peak intensities in bottom-up MS-based proteomics
#'       using singular value decomposition" Karpievitch YV, Taverner T, et al.
#'       2009, Bioinformatics
#' Ref:  "Metabolomics data normalization with EigenMS"
#'       Karpievitch YK, Nikolic SB, Wilson R, Sharman JE, Edwards LM 2014,
#'       PLoS ONE
#' @param m number of peptides x number of samples matrix of log-transformed
#'          expression data, metadata not included in this matrix
#' @param treatment either a single factor indicating the treatment group of
#'           each sample i.e. [1 1 1 1 2 2 2 2...]
#'           or a data frame of factors, eg:
#'            treatment= data.frame(cbind(data.frame(Group), data.frame(Time))
#' @param prot.info 2+ colum data frame, pepID, prID columns IN THAT ORDER.
#'              IMPORTANT: pepIDs must be unique identifiers and will be used
#'              as Row Names
#'              If normalizing non-proteomics data, create a column such as:
#'              paste('ID_',seq(1:num_rows), sep='')
#'              Same can be dome for ProtIDs, these are not used for
#'              normalization but are kept for future analyses
#' @param write_to_file if a string is passed in, 'complete' peptides
#'              (peptides with NO missing observations)
#'              will be written to that file name
#' @return A structure with multiple components
#' \describe{
#'   \item{m, treatment, prot.info, grp}{initial parameters passed into the
#'        function, returned for future reference}
#'   \item{my.svd}{matrices produced by SVD}
#'   \item{pres}{matrix of peptides that can be normalized, i.e. have enough
#'               observations for ANOVA}
#'   \item{n.treatment}{number of factors passed in}
#'   \item{n.u.treatment}{number of unique treatment facotr combinations, eg:
#'                   Factor A: a a a a c c c c
#'                   Factor B: 1 1 2 2 1 1 2 2
#'                   then:  n.treatment = 2; n.u.treatment = 4}
#'   \item{h.c}{number of bias trends identified}
#'   \item{present}{names/IDs of peptides in variable 'pres'}
#'   \item{complete}{complete peptides with no missing values,
#'                   these were used to compute SVD}
#'   \item{toplot1}{trends automatically identified in raw data,
#'                  can be plotted at a later time}
#'   \item{Tk}{scores for each bias trend, eigenvalues}
#'   \item{ncompl}{number of complete peptides with no missing observations}
#'}
#' @export
#' @examples
#' data(mm_peptides)
#' head(mm_peptides)
#' # different from parameter names as R uses outer name spaces
#' # if variable is undefined
#' intsCols = 8:13
#' metaCols = 1:7 # reusing this variable
#' m_logInts = make_intencities(mm_peptides, intsCols)  # will reuse the name
#' m_prot.info = make_meta(mm_peptides, metaCols)
#' m_logInts = convert_log2(m_logInts)
#' # 3 samples for CG and 3 for mCG
#' grps = as.factor(c('CG','CG','CG', 'mCG','mCG','mCG'))
#' mm_m_ints_eig1 = eig_norm1(m=m_logInts,treatment=grps,prot.info=m_prot.info)
#' mm_m_ints_eig1$h.c # check the number of bias trends detected
#' mm_m_ints_norm = eig_norm2(rv=mm_m_ints_eig1)
#' @export
eig_norm1 = function(m, treatment, prot.info, write_to_file=''){
  print("Data dimentions: ")
  print(dim(m))
  # check if treatment is a 'factor' vs data.frame',
  # i.e. single vs multiple factors
  if(class(treatment) == "factor") { # TRUE if one factor
     n.treatment = 1 # length(treatment)
     n.u.treatment = length(unique(treatment))[1]
  } else { # data.frame
    n.treatment = dim(treatment)[2]
    # all possible treatment combinations
    n.u.treatment = dim(unique(treatment))[1]
  }
  # convert m to a matrix from data.frame
  m = as.matrix(m) # no loss of information

  # filter out min.missing, here just counting missing values
  # if 1+ treatment completely missing, cannot do ANOVA,
  # thus cannot preserve grp diff.
  # IMPORTANT: we create a composite grp = number of unique combinations
  # of all groups, only for 'nested' groups for single layer
  # group is left as it is
  grpFactors = treatment # temporary var, leftover from old times...

  # length(colnames(treatment)) # not good: dim(grpFactors)
  nGrpFactors = n.treatment
  if(nGrpFactors > 1) { # got nested factors
    ugrps = unique(grpFactors)
    udims = dim(ugrps)
    grp = NULL
    for(ii in seq(1,udims[1]) ) {
      pos = grpFactors[,1] == ugrps[ii,1] # set to initial value
      for(jj in seq(2, udims[2]) ) {
	    pos = pos & grpFactors[,jj] == ugrps[ii,jj]
      }
	  grp[pos] = rep(ii, sum(pos))
    }
	grp = as.factor(grp)
  } else {
    grp = treatment
  }
  # nobs = number of observations
  nobs = array(NA, c(nrow(m), length(unique(grp))))
  print('Treatment groups:')
  print(grp)

  for(ii in seq_len(nrow(m) ) ) {
    for(jj in seq_len(length(unique(grp))) ) {
        # total number of groups num(g1) * num(g2) * ...
      nobs[ii,jj] = sum(!is.na(m[ii, grp==unique(grp)[jj]]))
    }
  }
  # now 'remove' peptides with missing groups
  present.min = apply(nobs, 1, min) # number present in each group
  ii = present.min == 0   # 1+ obs present in ALL of the groups
  pmiss = rbind(m[ii,]) # these have 1+ grp missing
  # rownames must be UNIQUE, if have possible duplicates: use 'ii' ?
  rownames(pmiss) = prot.info[ii,1]  # set rownames,

  # create matrix for peptides with enough observations for ANOVA
  # 'present' are names of the peptides (pepID) and 'pres' are abundances
  # NOTE: ! negates the proteins, so we get ones that have 1+ obs in each group
  present = prot.info[which(!prot.info[,1] %in% rownames(pmiss)), ] #rownames OK
  pres = m[which(!prot.info[,1] %in% rownames(pmiss)), ]
  rownames(pres) = prot.info[which(!prot.info[,1] %in% rownames(pmiss)),1]

  print('Selecting complete peptides')
  # Should issue an error message if we have NO complete peptides (unlikely)
  # select only 'complete' peptides, no missing values
  nobs = array(NA, nrow(pres)) # reassign noobs to dims of 'present'
  numiter = nrow(pres)
  for (ii in seq(1, numiter) ) {
    nobs[ii] = sum(!is.na(pres[ii,]))
  }

  iii = nobs == ncol(pres)
  complete = rbind(pres[iii,])

  #  write out a file of complete peptides if file name is passed in
  if(write_to_file != '') {
    utils::write.table(complete, file = write_to_file, append = FALSE,
             quote = FALSE, sep = "\t",
             eol = "\n", na = "NaN", dec = ".", row.names = TRUE,
             col.names = TRUE, qmethod = c("escape", "double"))
  }

  # compute bias with 'complete' matrix and residuals from 'present'
  # calculate eigenpeptides for 'complete' data only
  # if have only 1 group, we do not need to preserve group differernces,
  # everything is the same group, ex: QC samples
  # contrasts will fail if have only 1 group, thus have else
  if(n.u.treatment > 1) {
    print('Got 2+ treatment grps')
	  # using general function that can accomodate for 1+ number of factors
	  lm.fm = makeLMFormula(treatment, 'TREATS')
    TREATS = treatment
    # temp var to work if we got only 1 treatment vector.
    TREATS = data.frame(treatment)
    if(class(treatment) == "factor") {
      colnames(TREATS) = "TREATS"
     } else {
      colnames(TREATS) = colnames(treatment)
    }
    attach(TREATS)

    mod.c = stats::model.matrix(lm.fm$lm.formula, data=TREATS,
                         eval(parse(text=lm.fm$lm.params)))
    Y.c = as.matrix(complete) # used in eval() !!
	  options(warn = -1)

    # use stats::lm() to get residuals
    formula1 = paste('t(Y.c)~', as.character(lm.fm$lm.formula)[2], sep = '')
    TREATS = treatment
    fit_lmAll = stats::lm(eval(parse(text=formula1)))
    R.c = stats::residuals(fit_lmAll)  # Oct 2 messing with residuals...
  } else {  # 1 group only, set residuals to original matrix
    print('Got 1 treatment grp')
      mod.c = as.numeric(t(treatment))
	  # needs to be transposed to match the matrix returned from lm
	  R.c = t(as.matrix(complete))
      TREATS = treatment
  }

  print('Computing SVD, estimating Eigentrends...')
  # residuals are centered around 0, here center samples not
  # peptides/metabolites centering is basic normalization

  R.c_center = scale(R.c, center = TRUE, scale = FALSE)
  my.svd = svd(R.c_center)
  # can use wrapper below to chek if SVD has a problem...
  temp = my.svd$u
  my.svd$u = my.svd$v
  my.svd$v = temp

  #identify number of eigenvalues that account for
  # a significant amount of residual variation
  # save to return to the user as part of the return list
  numcompletepep = dim(complete)[1]
  # this is important info for publications
  # tell users how many peptides/metabolites the trends are based on
  # can also be determined by doing dim(return_value_fromEIg_norm1$pres)

  print(paste('Number of treatments: ', n.u.treatment))

  h.c = sva.id(complete, n.u.treatment, lm.fm=lm.fm, seed=1234)$n.sv
  print(paste("Number of significant eigenpeptides/trends", h.c) )

  # show RAW trends
  # center each peptide around zero (subtract its mean across samples)
  complete_center = scale(t(complete), center = TRUE, scale = FALSE)
  print('Preparing to plot...')

  n.u.treatment
  toplot1 = svd(complete_center) # scales above
  temp = toplot1$u
  toplot1$u = toplot1$v
  toplot1$v = temp

  graphics::par(mfcol=c(3,2))
  graphics::par(mar = c(2,2,2,2))
  plot.eigentrends(toplot1, "Raw Data")
  plot.eigentrends(my.svd, "Residual Data")

  d = my.svd$d;  ss = d^2;
  Tk = signif(ss/sum(ss)* 100, 2)

  retval = list(m=m, treatment=treatment, my.svd=my.svd,
            pres=pres, n.treatment=n.treatment, n.u.treatment=n.u.treatment,
            h.c=h.c, present=present, prot.info=prot.info,
            complete=complete, toplot1=toplot1, Tk=Tk, ncompl=numcompletepep,
            grp=grp)
  return(retval)
}

#' EigenMS normalization
#'
#' Eliminate the effects of systematic bias identified in eig_norm1()
#' Ref: "Normalization of peak intensities in bottom-up MS-based proteomics
#'       using singular value decomposition" Karpievitch YV, Taverner T et al.
#'       2009, Bioinformatics
#' Ref:  "Metabolomics data normalization with EigenMS"
#'       Karpievitch YK, Nikolic SB, Wilson R, Sharman JE, Edwards LM
#'       Submitted to PLoS ONE.
#'
#' @param rv return value from the eig_norm1 if user wants to change
#'           the number of bias trends that will be eliminated h.c in
#'           rv should be updates to the desired number
#' @return A structure with multiple components
#' \describe{
#'   \item{normalized}{matrix of normalized abundances with 2 columns
#'                     of protein and peptdie names}
#'   \item{norm_m}{matrix of normalized abundances, no extra columns}
#'   \item{eigentrends}{trends found in raw data, bias trends up to h.c}
#'   \item{norm.svd}{trends in normalized data, if one
#'                   wanted to plot at later time}
#'   \item{exPeps}{peptides excluded due to not enough peptides or
#'                 exception in fitting a linear model}
#'}
#' @examples
#' data(mm_peptides)
#' head(mm_peptides)
#' # different from parameter names as R uses outer name
#' # spaces if variable is undefined
#' intsCols = 8:13
#' metaCols = 1:7 # reusing this variable
#' m_logInts = make_intencities(mm_peptides, intsCols)
#' m_prot.info = make_meta(mm_peptides, metaCols)
#' m_logInts = convert_log2(m_logInts)
#' grps = as.factor(c('CG','CG','CG', 'mCG','mCG','mCG'))
#' mm_m_ints_eig1 = eig_norm1(m=m_logInts,treatment=grps,prot.info=m_prot.info)
#' mm_m_ints_eig1$h.c # check the number of bias trends detected
#' mm_m_ints_norm = eig_norm2(rv=mm_m_ints_eig1)
#' @export
eig_norm2 = function(rv) {
  # yuliya: use pres matrix, as we cannot deal with m anyways,
  # need to narrow it down to 'complete' peptides
  treatment = rv$treatment
  my.svd = rv$my.svd
  pres = rv$pres
  n.u.treatment = rv$n.u.treatment
  print(paste('Unique number of treatment combinations:', n.u.treatment) )
  h.c = rv$h.c
  present = rv$present
  toplot1 = rv$toplot1
  # vector of indicators of peptides that threw exeptions
  exPeps = vector(mode = "numeric", length = nrow(pres))

  print("Normalizing...")
  treatment = data.frame(treatment) # does this need to be done?
  if(n.u.treatment > 1) {
    lm.fm = makeLMFormula(treatment, 'ftemp')
    mtmp = stats::model.matrix(lm.fm$lm.formula, data=treatment,
                        eval(parse(text=lm.fm$lm.params)))
  } else {  # have 1 treatment group
    mtmp = treatment # as.numeric(t(treatment))
  }
  # above needed to know how many values will get back for some matrices
  # create some variables:
  betahat = matrix(NA,nrow=dim(mtmp)[2],ncol=nrow(pres))
  newR = array(NA, c(nrow(pres), ncol(pres))) #, n.treatment))
  norm_m = array(NA, c(nrow(pres), ncol(pres))) # , n.treatment))
  numsamp = dim(pres)[2]
  betahat_n = matrix(NA,nrow=dim(mtmp)[2],ncol=nrow(pres))
  rm(mtmp)

  V0 = my.svd$v[,seq(1,h.c),drop=FALSE]   # residual eigenpeptides

  if(n.u.treatment == 1) { # got 1 treatment group
    for (ii in seq(1,nrow(pres)) ) {
      if(ii%%250 == 0) { print(paste('Processing peptide ',ii))  }
      pep = pres[ii, ]
      pos = !is.na(pep)
      peptemp = as.matrix(pep[pos]) # take only the observed values
      resm = rep(NA, numsamp)
      resm[pos] = as.numeric(pep[pos])
      bias = array(NA, numsamp)
      bias[pos] = resm[pos] %*% V0[pos,] %*% t(V0[pos,])
      norm_m[ii, ] = as.numeric(pep - bias)
    }

  } else { # got 2+ treatment groups
    for (ii in seq(1, nrow(pres)) ) {
      if(ii %% 100 == 0) { print(paste('Processing peptide ',ii))  }
      pep = pres[ii, ]
      pos = !is.na(pep)
      # take only the observed values, may not be needed in R? but this works
      peptemp = as.matrix(pep[pos])
      ftemp = treatment[pos,]
      ftemp = data.frame(ftemp)
      #### use try, not entirely sure if need for modt, need it for solve lm?!
      options(warn = -1)
      # using general function that can accomodate for 1+ number of factors
      lm.fm = makeLMFormula(ftemp, 'ftemp')
      modt = try(stats::model.matrix(lm.fm$lm.formula, data=ftemp,
                              eval(parse(text=lm.fm$lm.params))), silent=TRUE)
      options(warn = 0)

      # do nothing if could not make model matrix
      if(!inherits(modt, "try-error")) {
        options(warn = -1)
        # if we are able to solve this, we are able to estimate bias
        bhat =  try(solve(t(modt) %*% modt) %*% t(modt) %*% peptemp)
        options(warn = 0)
        if(!inherits(bhat, "try-error")) {
          betahat[,ii] = bhat
          # these are the group effects, from estimated coefficients betahat
          ceffects = modt %*% bhat

          resm = rep(NA, numsamp) # really a vector only, not m
          resm[pos] = as.numeric(pep[pos] - ceffects)
          bias = array(NA, numsamp)
          bias[pos] = resm[pos] %*% V0[pos,] %*% t(V0[pos,])
          norm_m[ii, ] = as.numeric(pep - bias)

          # yuliya: newR should be computed on Normalized data
          resm_n = rep(NA, numsamp)
          bhat_n =  solve(t(modt) %*% modt) %*% t(modt) %*% norm_m[ii, pos]
          betahat_n[,ii] = bhat_n
          resm_n[pos] = norm_m[ii,pos] - ceffects
          newR[ii, ] = resm_n
        } else {
          print(paste('got exception 2 at peptide:', ii,
                      'should not get here...'))
          exPeps[ii] = 2
        }
      } else {
        print(paste('got exception at peptide:', ii))
        # keep track of peptides that threw exeptions, check why...
        exPeps[ii] = 1
      }
    }
  } # end else - got 2+ treatment groups

  ############################################################################
  # rescaling has been eliminated form the code after discussion that bias
  # adds variation and we remove it, so no need to rescale after as we removed
  # what was introduced
  y_rescaled = norm_m # for 1 group normalization only, we do not rescale
  # add column names to y-rescaled, now X1, X2,...
  colnames(y_rescaled) = colnames(pres) # these have same number of cols
  rownames(y_rescaled) = rownames(pres)
  y_resc = data.frame(present, y_rescaled)
  rownames(y_resc) = rownames(pres)  # rownames(rv$normalized)
  final = y_resc # row names are assumed to be UNIQUE, peptide IDs are unique

  # rows with all observations present
  complete_all = y_rescaled[rowSums(is.na(y_rescaled))==0,,drop=FALSE]

  #  x11() # make R open new figure window
  graphics::par(mfcol=c(3,2))
  graphics::par(mar = c(2,2,2,2))
  # center each peptide around zero (subtract its mean across samples)
  # note: we are not changing matrix itself, only centerig what we pass to svd
  complete_all_center = t(scale(t(complete_all), center = TRUE, scale = FALSE))
  toplot3 = svd(complete_all_center)
  plot.eigentrends(toplot1, "Raw Data")
  plot.eigentrends(toplot3, "Normalized Data")

  print("Done with normalization!!!")
  colnames(V0) =  paste("Trend", seq(1,ncol(V0)), sep="_")

  return(list(normalized=final, norm_m=y_rescaled, eigentrends=V0,
              norm.svd=toplot3, exPeps=exPeps))
} # end function eig_norm2



# EigenMS helper functions
# Tom had Sig set to 0.1, I and Storey's paper has 0.05
# sva.id = function(dat, mod, n.u.treatment, B=500, sv.sig=0.05, seed=NULL) {
# treatment, #' @param treatment treatment group information included in the
# analysis


#' Surrogate Variable Analysis
#'
#' Surrogate Variable Analysis function used internatlly by
#' eig_norm1 and eig_norm2
#' Here we incorporate the model matrix from EigenMS
#' normalization to find the significant
#' trends in the matrix of residuals.
#'
#' @param dat number of peptides/genes x number of samples
#'            matrix of expression data with no missing values
#' @param n.u.treatment number of treatment groups
#' @param lm.fm formular for treatment to be use on the right side of the call
#'        to stats::lm() as generated by makeLMFormula()
#' @param B The number of null iterations to perform
#' @param sv.sig The significance cutoff for the surrogate variables
#' @param seed A seed value for reproducible results
#' @return A data structure with the following values:
#' \describe{
#'   \item{n.sv}{Number of significant surrogate variables}
#'   \item{p.sv}{Significance for the returned surrogate variables}
#' }
sva.id = function(dat, n.u.treatment, lm.fm, B=500, sv.sig=0.05, seed=NULL)
{
  print("Number of complete peptides (and samples) used in SVD")
  print(dim(dat))


  if(!is.null(seed))  { set.seed(seed) }
  n = ncol(dat)

  ncomp = n.u.treatment # JULY 2013: as.numeric(n.u.treatment)
  print(paste("Number of treatment groups (in svd.id): ", ncomp))
  # should be true for either case and can be used later

  if(ncomp > 1) { #
    formula1 = paste('t(dat)~', as.character(lm.fm$lm.formula)[2], sep = '')
    fit_lmAll = stats::lm(eval(parse(text=formula1)))
    res = t(stats::residuals(fit_lmAll))
  } else {
    res = dat
  }
  # center each peptide around zero (subtract its mean across samples)
  res_center = t(scale(t(res), center = TRUE, scale = FALSE))

  uu = svd(t(res_center)) # NEED a WRAPPER for t(). the diag is min(n, m)
  temp = uu$u
  uu$u = uu$v
  uu$v = temp

  s0 = uu$d
  s0 = s0^2
  dstat = s0/sum(s0)
  ndf = length(dstat)
  dstat0 = matrix(0,nrow=B,ncol=ndf) # num samples (?)

print("Starting Bootstrap.....")
# Bootstrap procedure that determines the number of significant eigertrends...
  for(ii in seq(1, B)) {
    if(ii %% 50 == 0) { print(paste('Iteration ', ii)) }
    res0 = t(apply(res, 1, sample, replace=FALSE)) # regression
    # center each peptide around zero (subtract its mean across samples)
    # note: we are not changing matrix itself, only centerig what we pass to svd
    res0_center = t(scale(t(res0), center = TRUE, scale = FALSE))
    uu0 = svd(res0_center)
  	temp = uu0$u  # why did tom do this??
    uu0$u = uu0$v
    uu0$v = temp

    ss0 = uu0$d  # no need for diag....
    ss0 = ss0^2
    dstat0[ii,] = ss0 / sum(ss0) # Tk0 in Matlab
  }

# yuliya: check p-values here, Tom had mean value...
  psv = rep(1,n)
  for(ii in seq(1,ndf)) {
    # should this be compared to a MEAN?  Should this be dstat0[ii,] ?
    posGreater = dstat0[,ii] > dstat[ii]
    psv[ii] = sum(posGreater) / B
  }

  # p-values for trends have to be in the monotonically increasing order,
  # set equal to previous one if not the case
  for(ii in 2:ndf){
    if(psv[(ii-1)] > psv[ii]) {
      psv[ii] = psv[(ii-1)]
    }
  }
  nsv = sum(psv <= sv.sig)
  return(list(n.sv = nsv,p.sv=psv))
}
# end sva.id


mmul = function(A, B){
  # multiply square matrix by rectangle with NA's (???)
  X = A
  Y = B
  X[is.na(A)] = Y[is.na(B)] = 0
  R = X %*% Y
  R[is.na(B)] = NA
  R
}
# end mmul




#######################################################################
#######################################################################
my.Psi = function(x, my.pi){
  # calculates Psi
  exp(log(1-my.pi) + stats::dnorm(x, 0, 1, log=TRUE) -
        log(my.pi + (1-my.pi) * pnorm(x, 0, 1)))
}
# end my.Psi

my.Psi.dash = function(x, my.pi){
# calculates the derivative of Psi
-my.Psi(x, my.pi) * (x + my.Psi(x, my.pi))
}
# end my.Psi.dash

phi = function(x){ stats::dnorm(x) }

