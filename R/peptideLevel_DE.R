# Differetial Expression Nalsyis for follow after Model-based imputaion
#
# Ref:  "A statistical framework for protein quantitation in bottom-up MS-based
#        proteomics. Karpievitch Y, Stanley J, Taverner T, Huang J, Adkins JN,
#        Ansong C, Heffron F, Metz TO, Qian WJ, Yoon H, Smith RD, Dabney AR.
#        Bioinformatics 2009
#
# Written by Yuliya Karpievitch


#' Model-Based differential expression analysis
#'
#' Model-Based differential expression analysis is performed on peptide
#' level as desribed in
#' Karpievitch et al. 2009 "A statistical framework for protein quantitation in
#' bottom-up MS-based proteomics" Bioinformatics.
#' @param mm m x n matrix of intensities, number of peptides x number of samples
#' @param treatment vector indicating the treatment group
#'                  of each sample ie [1 1 1 1 2 2 2 2...]
#' @param prot.info 2+ colum data frame of peptide ID, protein ID, etc. columns
#' @param pr_ppos - column index for protein ID in prot.info.
#'                  Can restrict to be #2...
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{ProtID}{protein identification information taken from prot.info,
#'                 1 column used to identify proteins}
#'   \item{FC}{fold change}
#'   \item{p-value}{p-value for the comparison between 2 groups
#'                 (2 groups only here)}
#'   \item{BH-adjusted p-value}{Benjamini-Hochberg adjusted p-values}
#'}
#' @examples
#' data(mm_peptides)
#' head(mm_peptides)
#' # different from parameter names as R uses outer
#' # name spaces if variable is undefined
#' intsCols = 8:13
#' metaCols = 1:7 # reusing this variable
#' m_logInts = make_intencities(mm_peptides, intsCols)
#' m_prot.info = make_meta(mm_peptides, metaCols)
#' m_logInts = convert_log2(m_logInts)
#' grps = as.factor(c('CG','CG','CG', 'mCG','mCG','mCG'))
#' mm_m_ints_eig1 = eig_norm1(m=m_logInts,treatment=grps,prot.info=m_prot.info)
#' mm_m_ints_eig1$h.c # check the number of bias trends detected
#' mm_m_ints_norm = eig_norm2(rv=mm_m_ints_eig1)
#' mm_prot.info = mm_m_ints_norm$normalized[,1:7]
#' mm_norm_m =  mm_m_ints_norm$normalized[,8:13]
#' imp_mm = MBimpute(mm_norm_m, grps, prot.info=mm_prot.info,
#'                   pr_ppos=2, my.pi=0.05,
#'                   compute_pi=FALSE, sseed=131)
#' DE_res = peptideLevel_DE(imp_mm$y_imputed,
#'                          grps, mm_m_ints_norm$normalized[,metaCols],
#'                          pr_ppos=2)
#' @export
peptideLevel_DE = function(mm, treatment, prot.info, pr_ppos=2)
{
  # Match to protein
  all.proteins = unique(prot.info[,pr_ppos])
  numProts = length(all.proteins) # 1569
  y_out = matrix(NA, numProts, 5)
  de_ret = NULL
  u_prot_info = NULL
  ll = length(all.proteins)
  for (kk in seq(1,ll)) {
    prot = all.proteins[kk]
    pmid.matches = prot.info[prot.info[,pr_ppos]==prot,1]
    curr_prot.info = prot.info[prot.info[,pr_ppos]==prot,]
    idx.prot = which(prot.info[,1] %in% pmid.matches)
    # need to return unique prot.info, make it as we go
    ttt = prot.info[idx.prot,]
    if(!is.null(dim(ttt))) {
      u_prot_info = rbind(u_prot_info, ttt[1,])
    } else {
      u_prot_info = rbind(u_prot_info, ttt)
    }
    y_raw = mm[idx.prot,,drop=FALSE]
    n.peptide = nrow(y_raw)
    yy = as.vector(t(y_raw))
    ## n = length(yy)
    peptide = as.factor(rep(seq(1,n.peptide),
                            each=dim(data.frame(treatment))[1])) # used below
    # keep track of prIDs here
    curr_prot.info = curr_prot.info[kk,] # possibly a subset

    # good to know how many peptides were in a given protein
    y_out[kk,5] = n.peptide

    if (n.peptide != 1){
      # replicate treatment for # peptides
      treatment_hold = treatment
      treatment = rep(treatment, times=n.peptide)
      res = stats::lm(yy ~ treatment + peptide,
                      contrasts = list(peptide="contr.sum"))
      res1 = summary(res) # above line only good for 2 groups?
      # col 3 reserved for BH p-values
      y_out[kk,c(1,2,4)] = stats::coefficients(res1)[2,c(1,4,3)]
      treatment = treatment_hold
    } else {
      res =  stats::lm(yy~treatment) # only good for 2 groups?
      res1 = summary(res)
      # col 3 reserved for BH p-values
      y_out[kk,c(1,2,4)] = stats::coefficients(res1)[2,c(1,4,3)]
    }                      # coefs, p-value, t-value
  } # end for each protein

  colnames(y_out) = c('FC', 'P_val', 'BH_P_val', 'statistic', 'num_peptides')
  # BH adjustment
  y_out[,3] = stats::p.adjust(y_out[,2],"BH")

  # add protein names as 1st col in a data frame
  DE_res = data.frame(all.proteins, y_out, stringsAsFactors=FALSE)
  de_ret$DE_res = DE_res
  de_ret$prot.info = u_prot_info

  cols1 = colnames(de_ret$DE_res)
  cols1[1] = "ProtIDused"
  cols2 = colnames(de_ret$prot.info)
  de_ret = data.frame(de_ret, stringsAsFactors = FALSE)
  colnames(de_ret) = c(cols1, cols2)
  return(de_ret)
}


#' Plot peptide trends
#'
#' Plot Raw, Normalized and Normalized & Imputed peptide trends for a protein
#' @param mm matrix of raw intensities
#' @param prot.info metadata for the intensities in mm
#' @param sorted_norm_m normalized intensities, possibly fewer than in mm due to
#'        filtering out peptides with fewer than one
#'        observation per treatment group
#' @param sorted_prot.info metadata for the intensities in sorted_norm_m
#' @param imp_mm imputed intensities (post normalization)
#' @param imp_prot.info metadata for the imputed intensities in imp_mm
#' @param prot_to_plot protein ID to plot
#' @param prot_to_plot_col protein ID column index
#' @param gene_name gene ID to plot
#' @param gene_name_col gene ID to plot column index
#' @param mylabs sample labels to be plotted on x-axis
#'
#' @examples
#'
#' data("hs_peptides") # loads variable hs_peptides
#' intsCols = 8:13 # column indeces that contain intensities
#' m_logInts = make_intencities(hs_peptides, intsCols)
#' # replace 0's with NA's as NA's are more appropriate
#' # for anlysis and log2 transform
#' m_logInts = convert_log2(m_logInts)
#' # column indices that contain metadata such as protein IDs and sequences
#' metaCols = 1:7
#' m_prot.info = make_meta(hs_peptides, metaCols)
#' grps = as.factor(c('CG','CG','CG', 'mCG','mCG','mCG'))
#' hs_m_ints_eig1 = eig_norm1(m=m_logInts,treatment=grps,prot.info=m_prot.info)
#' hs_m_ints_eig1$h.c = 2 # visually looks like there are 2 bias trends at least
#' hs_m_ints_norm = eig_norm2(rv=hs_m_ints_eig1)
#' hs_prot.info = hs_m_ints_norm$normalized[,metaCols]
#' hs_norm_m =  hs_m_ints_norm$normalized[,intsCols]
#' imp_hs = MBimpute(hs_norm_m, grps, prot.info=hs_prot.info,
#'                   pr_ppos=3, my.pi=0.05, compute_pi=FALSE, sseed=171717)
#' mylabs = c( 'CG','CG','CG', 'mCG','mCG','mCG')
#' prot_to_plot = 'Prot32' # 43
#' gene_to_plot = 'Gene32'
#' plot_3_pep_trends_NOfile(as.matrix(hs_m_ints_eig1$m),
#'                          hs_m_ints_eig1$prot.info,
#'                          as.matrix(hs_norm_m),
#'                          hs_prot.info,
#'                          imp_hs$y_imputed,
#'                          imp_hs$imp_prot.info,
#'                          prot_to_plot, 3,
#'                          gene_to_plot, 4, mylabs)
#'
#' @return Nil
#' @export
plot_3_pep_trends_NOfile = function(mm, prot.info, sorted_norm_m,
                                    sorted_prot.info,
                                    imp_mm, imp_prot.info, prot_to_plot,
                                    prot_to_plot_col,
                                    gene_name, gene_name_col, mylabs) {
  graphics::par(mfcol=c(1,3))
  ppos = imp_prot.info[,prot_to_plot_col] == prot_to_plot
  tmp = imp_mm[ppos,]
  if(sum(ppos) == 1) { # only 1 row, duplicate 1 row...
    tmp = rbind(tmp,tmp)
  }
  # take a Raw data, need these for y-limits
  # unsorted but OK, peptides will be in different order
  ppos2 = prot.info[,prot_to_plot_col] == prot_to_plot
  tmp2 = mm[ppos2,]
  if(sum(ppos2) == 1) { # only 1 row, duplicate 1 row...
    tmp2 = rbind(tmp2,tmp2)
  }
  # could not normalize as all obs missing in 1+ group...
  if(dim(tmp)[1] == 0 ) { tmp = tmp2 }
  ylim_min = min(min(tmp, na.rm=TRUE), min(tmp2, na.rm=TRUE))
  ylim_max = max(max(tmp, na.rm=TRUE), max(tmp2, na.rm=TRUE))
  myylim = c(ylim_min-.2, ylim_max+.2) # give a bit of margin

  # print(paste(myylim, sep=' '))
  main_title = paste(gene_name, ' (',
                     prot_to_plot, ') Normalized & Imputed', sep='')
  graphics::par(mfcol=c(1,3))
  graphics::par(mar=c(10,3,3,3))
  # may need to transpose ylim will be dynamically set here,
  # and kept for the other plots
  graphics::matplot(t(tmp), type="l", main=main_title, xaxt = "n", ylim=myylim)
  graphics::axis(1, at=seq(1,length(mylabs)), labels=mylabs, las=3)
  graphics::matpoints(t(tmp), type="p", pch='*') # , ylim=c(15,35))
  graphics::lines(colMeans(tmp,na.rm = TRUE), lwd=3, col='black')
  limits = graphics::par("usr")  # upper and low ylim and upper and lower xlim
  myxlim = c(limits[1], limits[2])

  ppos = sorted_prot.info[,prot_to_plot_col] == prot_to_plot  #
  tmp = sorted_norm_m[ppos,]
  if(sum(ppos) == 1) { # only 1 row, duplicate 1 row...
    tmp = rbind(tmp,tmp)
  }
  # could not normalize as all obs missing in 1+ group...
  if(dim(tmp)[1] == 0 ) { tmp = tmp2 }
  main_title = paste(gene_name, ' (', prot_to_plot, ') Normalized', sep='')
  graphics::matplot(t(tmp), type="l", pch='19',  main=main_title,
          ylim=myylim, xlim=myxlim, xaxt = "n") # may need to transpose
  graphics::axis(1, at=seq(1,length(mylabs)), labels=mylabs, las=3)
  graphics::matpoints(t(tmp), type="p", pch='*')
  graphics::lines(colMeans(tmp,na.rm = TRUE), lwd=3, col='black')


  main_title = paste(gene_name, ' (', prot_to_plot, ') Raw', sep='')
  # may need to transpose
  graphics::matplot(t(tmp2), type="l", main=main_title, ylim=myylim, xaxt = "n")
  graphics::axis(1, at=seq(1,length(mylabs)), labels=mylabs, las=3)
  graphics::matpoints(t(tmp2), type="p", pch='*') # , ylim=c(15,35))
  graphics::lines(colMeans(tmp2,na.rm = TRUE), lwd=3, col='black')
  graphics::par(mar=c(3,3,3,3))
}


#' Plot trends for a single protien
#'
#' Plot peptide trends for a protein
#' @param mm matrix of raw intensities
#' @param prot.info metadata for the intensities in mm
#' @param prot_to_plot protein ID to plot
#' @param prot_to_plot_col protein ID column index
#' @param gene_name gene ID to plot
#' @param gene_name_col gene ID to plot column index
#' @param colors what colors to plot peptide abundances as,
#'        most commonly should be
#'        treatment groups
#' @param mylabs sample labels to be plotted on x-axis
#' @return Nil
plot_1prot = function(mm, prot.info, prot_to_plot, prot_to_plot_col,
                           gene_name, gene_name_col, colors, mylabs) {
  graphics::par(mfcol=c(1,1))
  graphics::par(mar=c(10,3,3,3))
  ppos = prot.info[,prot_to_plot_col] == prot_to_plot
  tmp = mm[ppos,]
  if(sum(ppos) == 1) {
      # only 1 row, duplicate 1 row to be able to sued matplot function
    tmp = rbind(tmp,tmp)
  }
  # may need to transpose ylim will be dynamically set here,
  # and kept for the other plots
  graphics::matplot(t(tmp), type="l", main=prot_to_plot, xaxt = "n", col=colors)
  graphics::axis(1, at=seq(1,length(mylabs)), labels=mylabs, las=3)
  graphics::matpoints(t(tmp), type="p", pch='*', col=colors)
  graphics::lines(colMeans(tmp,na.rm = TRUE), lwd=3, col='black')
  graphics::par(mar=c(3,3,3,3))
}




plot_3_peptide_trends = function(mm, prot.info, sorted_norm_m,
                                 sorted_prot.info, imp_mm,
                                 imp_prot.info, prot_to_plot,
                                 gene_name, mylabs) {

  graphics::par(mfcol=c(1,3))
  ppos = imp_prot.info[,2] == prot_to_plot
  tmp = imp_mm[ppos,]
  if(sum(ppos) == 1) { # only 1 row, duplicate 1 row...
    tmp = rbind(tmp,tmp)
  }
  # take a Raw data, need these for y-limits
  # unsorted but OK, peptides will be in different order
  ppos2 = prot.info[,2] == prot_to_plot
  tmp2 = mm[ppos2,]
  if(sum(ppos2) == 1) { # only 1 row, duplicate 1 row...
    tmp2 = rbind(tmp2,tmp2)
  }
  if(dim(tmp)[1] == 0 ) { tmp = tmp2 }
  # could not normalize as all obs missing in 1+ group...
  ylim_min = min(min(tmp, na.rm=TRUE), min(tmp2, na.rm=TRUE))
  ylim_max = max(max(tmp, na.rm=TRUE), max(tmp2, na.rm=TRUE))
  myylim = c(ylim_min-.2, ylim_max+.2) # give a bit of margin
  main_title = paste(gene_name, ' (', prot_to_plot,
                     ') Normalized & Imputed', sep='')

  outfnames_png = paste(gene_name, '_', prot_to_plot, '_3pepTrends.png',sep='')
  grDevices::png(outfnames_png, width = 20, height = 6.4, units='in', res = 400)
  # R cannot figue out how & when res is specifed... (?)
  graphics::par(mfcol=c(1,3))
  graphics::matplot(t(tmp), type="l", main=main_title, xaxt = "n", ylim=myylim)
  # may need to transpose ylim will be dynamically set here,
  # and kept for the other plots
  graphics::axis(1, at=seq(1,length(mylabs)), labels=mylabs)
  graphics::matpoints(t(tmp), type="p", pch='*') # , ylim=c(15,35))
  graphics::lines(colMeans(tmp,na.rm = TRUE), lwd=3, col='black')
  limits = graphics::par("usr")  # upper and low ylim and upper and lower xlim
  myxlim = c(limits[1], limits[2])

  ppos = sorted_prot.info[,2] == prot_to_plot  #
  tmp = sorted_norm_m[ppos,]
  if(sum(ppos) == 1) { # only 1 row, duplicate 1 row...
    tmp = rbind(tmp,tmp)
  }
  # could not normalize as all obs missing in 1+ group...
  if(dim(tmp)[1] == 0 ) { tmp = tmp2 }
  main_title = paste(gene_name, ' (', prot_to_plot, ') Normalized', sep='')
  graphics::matplot(t(tmp), type="l", pch='19',  main=main_title,
          ylim=myylim, xlim=myxlim, xaxt = "n") # may need to transpose
  graphics::axis(1, at=seq(1,length(mylabs)), labels=mylabs)
  graphics::matpoints(t(tmp), type="p", pch='*')
  graphics::lines(colMeans(tmp,na.rm = TRUE), lwd=3, col='black')


  main_title = paste(gene_name, ' (', prot_to_plot, ') Raw', sep='')
  # may need to transpose
  graphics::matplot(t(tmp2), type="l", main=main_title, ylim=myylim, xaxt = "n")
  graphics::axis(1, at=seq(1,length(mylabs)), labels=mylabs)
  graphics::matpoints(t(tmp2), type="p", pch='*') # , ylim=c(15,35))
  graphics::lines(colMeans(tmp2,na.rm = TRUE), lwd=3, col='black')
  grDevices::dev.off()
}


# function takes a vector of string gene IDs possibly separated by a ';'
# and returns a vector of the same lenfth with only the first gene ID
get1sttoken = function(ids)
{
  ll = length(ids) # 22780
  prots = NULL
  for(ii in seq(1,ll)) {
    tt_spl = strsplit(ids[ii], ';')
    prots = c(prots, c(tt_spl[[1]][1]))
  }
  return(prots)
}


# function takes a vector of string gene IDs possibly separated by a ';'
# and returns a vector of the same lenfth with the first gene ID
# with a '+' if more than 1 gene ID was present
get1sttokenPlus = function(ids)
{
  ll = length(ids)
  prots = NULL
  for(ii in seq(1,ll)) {
    tt_spl = strsplit(ids[ii], ';')
    tt = tt_spl[[1]][1]
    if(length(tt_spl[[1]]) > 1) {
      tt = paste(tt, '+', sep='')
    }
    prots = c(prots, tt)
  }
  return(prots)
}


# function takes a vector of string Protein IDs
# and removes any -2 from the end of the protein name
# input is one protein per row
clip_protID_end = function(ids){
  ll = length(ids)
  prots = vector('character', length = ll)
  for(ii in seq(1,ll)) {
    #if(ii %/% 1000 != 0) { print(ii) }
    ppos = regexpr("-", ids[ii], fixed=TRUE)[1]
    if(ppos != -1) {
      prots[ii] = substr(ids[ii], 1, (ppos-1) )
    } else {
      prots[ii] = ids[ii]
    }
  }
  return(prots)
}


## function
## remove_contaminats(dd, colsCheck = c('Reverse',
##  'Potential.contaminant'), colsCheckSymbol=c('+', '+'))
## check for the presence of '+' sign in the
## specificed columns but column number
remove_contaminats_symbol = function(dd, colsCheck = c('Reverse',
                                                       'Potential.contaminant'),
                                     colsCheckSymbol=c('+', '+') )
{
  # get column indeces for the specified column names,
  # eg: 'Reverse' and 'Potential.contaminant'
  # access 'dd' by the column index to check for '+'
  cols1 = colnames(dd)
  ll = length(colsCheck)
  for(ii in seq(1,ll)) {
    ppos = cols1 == colsCheck[ii]
    # position of the column we need to check
    colToCheck = seq(1,length(cols1))[ppos]
    tt1 = dd[,colToCheck] == '+'
    print(paste('Removing ', sum(tt1), 'rows via ', colsCheck[ii], 'column') )
    dd = dd[!tt1,]
  }
  return(dd)
}


## function
## remove_contaminats(dd, colsCheck = c('Reverse', 'Potential.contaminant'),
## colsCheckSymbol=c('+', '+'))
## check for the presence of '+' sign in the
## specificed columns but column number
remove_contaminats_prefix = function(dd,
                                     colsCheck = c('Leading.razor.protein',
                                                   'Leading.razor.protein'),
                                     colsCheckPrefix=c('REV_', 'CON_') )
{
  # get column indeces for the specified column names,
  # eg: 'Reverse' and 'Potential.contaminant'
  # access 'dd' by the column index to check for matches
  # to the prefix in the protein names
  cols1 = colnames(dd)
  ll = length(colsCheck)
  for(ii in seq(1,ll)) {
    ppos = cols1 == colsCheck[ii]
    # position of the column we need to check
    colToCheck = seq(1,length(cols1))[ppos]
    ppos1 = startsWith(as.character(dd[,colToCheck]), colsCheckPrefix[ii])
    print(paste('Removing ', sum(ppos1), 'rows via ', colsCheck[ii], 'column'))
    dd = dd[!ppos1,]
  }
  return(dd)
}


# in case of GeneID column having multiple gene names separated by ';'
# take the first one and make it the GeneID column,
# store original GeneID column in GeneIDLong
generate_single_geneID = function(dd, colCheck)
{
  cols1 = colnames(dd)
  ppos = cols1 == colCheck
  colToCheck = seq(1,length(cols1))[ppos]
  dd$GeneIDLong = dd[,colToCheck]
  dd$GeneID = get1sttoken(as.character(dd[,colToCheck]) )
  return(dd)
}


# in case of ProtID column having -2 appended to the protein name
# (or any other -#), eg: Q86U42-2
# take the portion prior to -# and make it the ProtID column,
# store original ProtID column in GeneIDLong
# If there is not -2 or other number following the protein name,
# the name will be unchanged
# Protein names with -2 may create problems when matching tot he
# Ensembl IDs - multiple organizms only!
generate_clipped_ProtID = function(dd, colCheck)
{
  cols1 = colnames(dd)
  ppos = cols1 == colCheck
  colToCheck = seq(1,length(cols1))[ppos]
  dd$ProtIDLong = dd[,colToCheck]
  dd$ProtID = clip_protID_end(as.character(dd[,colToCheck]) )
  return(dd)
}


# function to use after we have fitered duplicate
# Protein and Gene IDs based on Ens ID
# on rare occasion there exist 2 EnsIDs for the same
# GeneID but for different ProtID
# or valid protein ID in one row and missing protein ID
# int he 2nd row. Never seem 3 rows...
remove_dup_geneIDs = function(db)
{
  ppos_genid = db$GeneID !='' # 18351
  db_geneID = db[ppos_genid,]
  dim(db_geneID) # 18351
  db_nogeneID = db[!ppos_genid,]
  dim(db_nogeneID) # 424    => dim(db) = 18775, OK

  ppos = duplicated(db_geneID$GeneID )
  sum(ppos) # 53
  dup_geneIDs = db_geneID[ppos,1] # 53 IDs

  ppos_dup = db_geneID$GeneID %in% dup_geneIDs
  length(ppos_dup) # 18351
  sum(ppos_dup)    # 106

  db_geneID_nodup = db_geneID[!ppos_dup,]
  dim(db_geneID_nodup) # 18245
  db_geneID_dup = db_geneID[ppos_dup,]
  dim(db_geneID_dup) # 106

  # eliminate rows with no protID in the dup ones
  # and concatenate back the matrices
  u_dup_genids = unique(db_geneID_dup$GeneID) # 53
  nn = length(u_dup_genids)
  undup_db = NULL
  for(ii in seq(1,nn)) {
    pcur = db_geneID_dup$GeneID == u_dup_genids[ii]
    cur_gene = db_geneID_dup[pcur,]
    # keep row with ProtID
    ppos_protid = cur_gene$ProtID != ''
    if(sum(ppos_protid) > 1) {
      # keep the bottm row - in our dataset we get
      # that as Leading Razor Protein...
      undup_db = rbind(undup_db,cur_gene[2,])
      # http://www.enzolifesciences.com/ADI-SPP-502/hsp70-a1-mouse-recombinant/
    } else {
      undup_db = rbind(undup_db,cur_gene[ppos_protid,])
    }
  }
  # combine db back together without the duplicates
  # db_geneID_nodup  db_nogeneID  undup_db
  db_recomb = rbind(db_geneID_nodup,  db_nogeneID,  undup_db)
  return(db_recomb)
}


###############################################################################
## functions to help add Ensembl IDs to various organizms

# for each GeneID take all of the unique EnsIDs
# (and ProtIDs?) and make into a list?
# make a data structure:
# HS$GeneID - vector of unique gene IDs, strings
# HS$EnsID  - vector of lists of all EnsIDs that match to this gene
# HS$ProtID - vector of lists of all ProtIDs that match to this gene
# HS$numEnsList - vector of values of number of elelments in teh list of EnsIDs
# must have 3 columns: GeneID, ProtID and EnsID obtained from Biomart
make_ID_structure = function(ens_ids) {
  ugenes = unique(ens_ids$GeneID)
  ll = length(ugenes) # 19086
  HS = NULL
  HS$GeneID = NULL
  HS$EnsID = vector(mode="list", length=ll)
  HS$ProtID = vector(mode="list", length=ll)
  HS$NumEnsList = vector(mode="numeric", length=ll)
  for(ii in seq(1,ll)) {
    HS$GeneID[ii] = ugenes[ii]
    ppos = ens_ids$GeneID == HS$GeneID[ii]
    HS$EnsID[[ii]] = ens_ids$EnsID[ppos]
    # we already dropped rows with NO ProtID, so should not have blank ProtIDs
    HS$ProtID[[ii]] = NULL
    HS$ProtID[[ii]] = ens_ids$ProtID[ppos]
    HS$NumEnsList[ii] = length(HS$EnsID[[ii]])
  }
  return(HS)
}

get_EnsIDs = function(organizm='human') {
  lookup = matrix('', nrow = 2, ncol = 3)
  lookup[1,] = c('human', 'hsapiens_gene_ensembl', 'hgnc_symbol')
  lookup[2,] = c('mouse', 'mmusculus_gene_ensembl', 'mgi_symbol')
  # to be continued...

  lookup_pos = lookup[,1] == organizm

  ensembl = biomaRt::useMart("ensembl", dataset=lookup[lookup_pos,2])
  gene_symbol = lookup[lookup_pos,3]
  res = biomaRt::getBM(attributes=c(gene_symbol, 'uniprotswissprot',
                                    'ensembl_gene_id'),
                       mart=ensembl)
  dim(res)  # human: 37593
  colnames(res) = c('GeneID', 'ProtID', 'EnsID')

  # remove rows with blank ProtID, will not be an option
  # on the proteomics dataset (?)
  ppos = res$ProtID == ''
  sum(ppos) # 16509
  res = res[!ppos,]
  ens_ids = make_ID_structure(res)

  return(ens_ids)
}


# dd1 is the datastructure returned from make_ID_structure()
# must contain dd1$EnsID column
# linker - 2 column matrix with  EnsIDs that match to EnsIDs
#          in dd1 as well as the other organizm EnsIDs
#          which will be added as a list to dd1
match_linker_ids = function(dd1,  linker, l_col) {
  ll1 = length(dd1$EnsID) # 19086
  ll1
  match1 = vector(mode="list", length=ll1)
  for(ii in seq(1,ll1)) {
    if(ii%%100 == 0) { print(paste('Processing Gene', as.character(ii) ) ) }
    ppos = linker[,l_col] %in% dd1$EnsID[[ii]]
    match1[[ii]] = unique(linker[ppos,-l_col])
    linker = linker[!ppos,]
  }
  dd1$linkedIDs = match1
  return(dd1)
}




# mm_ens_ids_save = mm_ens_ids
# 2 structures "EnsID"      "ProtID"     "NumEnsList" "GeneID"     "linkedIDs"
# dd1 will be matched by EsnID column and dd2 by linkedIDs
match_ids = function(dd1, dd2) {
  ll1 = length(dd1$GeneID) # 19086
  ll2 = length(dd2$GeneID) # 15032
  # will NA's instead of default 0's  # vector('numeric', length=ll1)
  match1 = rep(NA, each=ll1)
  match2 = rep(NA, each=ll2)  # vector('numeric', length=ll2)
  num_match1 = vector('numeric', length=ll1)
  match_index = 1
  # add index such that we could track the original position of
  # the gene when elements are removed
  dd2$index = seq(1,ll2)
  print('Processing Gene 1')
  for(ii in seq(1,ll1)) {
    if(ii%%500 == 0) { print(paste('Processing Gene', as.character(ii) ) ) }
    ll2 = length(dd2$GeneID) # will shrink as we go through the matching process

    for(jj in seq(1,ll2)) {
      ppos = dd1$EnsID[[ii]] %in% dd2$linkedIDs[[jj]]
      if(sum(ppos) > 0) { # we got a match
        num_match1[ii] = sum(ppos)
        match1[ii] = match_index
        match_2_index = dd2$index[jj]
        match2[match_2_index] = match_index
        # remove elements that matched from the structure # 2
        dd2$GeneID = dd2$GeneID[-jj]
        dd2$EnsID = dd2$EnsID[-jj]
        dd2$ProtID =  dd2$ProtID[-jj]
        dd2$NumEnsList = dd2$NumEnsList[-jj]
        dd2$index = dd2$index[-jj]
        dd2$linkedIDs = dd2$linkedIDs[-jj]
        match_index = match_index + 1
        break # out of inner loop hopefully
      }
    }
  }
  ret = NULL
  ret$match1 = match1
  ret$match2 = match2
  return(ret)
}


# use hs_with_MM & mm_with_HS
# add matched IDs to the structures for different organisms
# function uses match_ids() to get most of the work done
# dd1 and dd2 must be structures returned from match_linker_ids
#
add_match_ids = function(dd1, dd2) {
  matched_ids = match_ids(dd1, dd2)
  dd1$matchedID = matched_ids$match1
  dd2$matchedID = matched_ids$match2
  ret = NULL
  ret$dd1 = dd1
  ret$dd2 = dd2
  return(ret)
}
