################################################################################
# Function to conduct multi-sample two-part test with
# permutation null distribution
# Initial implementation date: 22 July 2016
# Releae data: February 2018
################################################################################

#' hs_peptides - peptide-level intensities for human
#'
#' A dataset containing the protein and petide information and peptide-level
#' intensities for 6 samples: 3 CG and 3 mCG groups. There are 69 proteins.
#' The columns are as follows:
#'
#' \itemize{
#'   \item Sequence - peptide sequence - randomly chosen from a larger list of
#'         sequences
#'   \item MatchedID - numeric ID that links proteins in the two datasets,
#'         unnecessary if datasets are for the same species
#'   \item ProtID - protein ID, artificial protein ID, eg. Prot1, Prot2, ...
#'   \item GeneID - gene ID, artificial gene ID, eg. Gene1, Gene2, ...
#'   \item ProtName - artificial Protein Name
#'   \item ProtIDLong - long protein ID, full protein name, here artificially
#'         simulated
#'   \item GeneIDLong - long gene ID, full gene name, here artificially
#'         simulated
#'   \item CG1 - raw intensity column for sample 1 in CG group
#'   \item CG2 - raw intensity column for sample 2 in CG group
#'   \item CG3 - raw intensity column for sample 3 in CG group
#'   \item mCG1 - raw intensity column for sample 1 in mCG group
#'   \item mCG2 - raw intensity column for sample 2 in mCG group
#'   \item mCG3 - raw intensity column for sample 3 in mCG group
#' }
#'
#' @docType data
#' @keywords datasets
#' @name hs_peptides
#' @usage data(hs_peptides)
#' @format A data frame with 695 rows and 13 colummns, compiring 7 columns of
#'        metadata and 6 columns of peptide intensities. 69 proteins.
NULL

#' mm_peptides - peptide-level intensities for mouse
#'
#' A dataset containing the protein and petide information and peptide-level
#' intensities for 6 samples: 3 CG and 3 mCG groups. There are 69 proteins.
#' The columns are as follows:
#'
#' \itemize{
#'   \item Sequence - peptide sequence - randomly chosen from a larger list of
#'         sequences
#'   \item MatchedID - numeric ID that links proteins in the two datasets,
#'         unnecessary if datasets are for the same species
#'   \item ProtID - protein ID, artificial protein ID, eg. Prot1, Prot2, ...
#'   \item GeneID - gene ID, artificial gene ID, eg. Gene1, Gene2, ...
#'   \item ProtName - artificial Protein Name
#'   \item ProtIDLong - long protein ID, full protein name, here artificially
#'         simulated
#'   \item GeneIDLong - long gene ID, full gene name, here artificially
#'         simulated
#'   \item CG1 - raw intensity column for sample 1 in CG group
#'   \item CG2 - raw intensity column for sample 2 in CG group
#'   \item CG3 - raw intensity column for sample 3 in CG group
#'   \item mCG1 - raw intensity column for sample 1 in mCG group
#'   \item mCG2 - raw intensity column for sample 2 in mCG group
#'   \item mCG3 - raw intensity column for sample 3 in mCG group
#' }
#'
#' @docType data
#' @keywords datasets
#' @name mm_peptides
#' @usage data(mm_peptides)
#' @format A data frame with 1102 rows and 13 colummns, compiring 7 columns of
#'         metadata and 6 columns of peptide intensities. 69 proteins.
NULL



#' Subdivide data into intensities columns only
#'
#' Subdivide a data frame of protein intensities and
#' metadata into intensities only.
#' No row names will be provided.
#'
#' @param mm data frame of metadata and intensities as a single data frame
#' @param use_cols column numbers to subset and return,
#'                 no range checking no range
#'                 checking on the column indeces is performed
#' @return m_ints data frame of intensities only
#' @examples
#' data(mm_peptides)
#' head(mm_peptides)
#' intsCols = 8:13 # different from parameter names as R uses outer name
#'                 # spaces if variable is undefined
#' m_logInts = make_intencities(mm_peptides, intsCols)
#'
#' @export
make_intencities = function(mm, use_cols) {
  m_ints = mm[,use_cols] # Sequences should be unique
  return(m_ints)
}

#' Subdivide data into metadata columns only
#'
#' Subdivide a data frame of protein metadata and intensities
#' into a data frame of meta data only
#'
#' @param mm data frame of metadata and intensities as a signle data frame
#' @param use_cols column numbers to subset and return,
#'                 no range checking on the column
#'                 indeces is performed
#' @return m_ints data frame of intensities only
#' @examples
#' data(mm_peptides)
#' head(mm_peptides)
#' metaCols = 1:7 # reusing this variable
#' m_prot.info = make_meta(mm_peptides, metaCols)
#' @export
make_meta = function(mm, use_cols) {
  m_meta = mm[,use_cols] # Sequences should be unique
  return(m_meta)
}

#' Convert values in a matrix to log2 transfored values
#'
#' convert_log2 replaces 0's with NA's than does a log2 transformation
#' Replacing 0's with NA's is the correct approach to Proteomics data analysis
#' as 0's are not values that should be left in the data where no
#' observation was made, see citation below.
#' Karpievitch et al. 2009 "Normalization of peak intensities in
#'   bottom-up MS-based proteomics using singular value decomposition".
#'   PMID: 19602524
#' Karpievitch et al. 2009 "A statistical framework for protein
#'   quantitation in bottom-up MS-based proteomics".  PMID: 19535538
#'
#' @param mm a dataframe of raw intensities in format:
#'  (# peptides)x(# samples + possibly peptide & protein information (metadata))
#'
#' @param use_cols vector of column indecies that make up the intensities
#'              usually in sequential order but do not have to be
#'              user is responsible for making sure that specified columns are
#'              indeed numeric and correspond to intensities for each sample
#'
#' @return matrix of log2 transforemd intensities where 0's were
#'         replaced with NA's prior
#'         to transformation
#'
#' @examples
#' data(mm_peptides)
#' head(mm_peptides)
#' intsCols = 8:13
#' metaCols = 1:7
#' m_logInts = make_intencities(mm_peptides, intsCols)
#' m_prot.info = make_meta(mm_peptides, metaCols)
#' m_logInts = convert_log2(m_logInts) # 0's replaced with NAs and
#'                                     # log2 transnform applied
#'
#' @export
convert_log2 = function(mm, use_cols) {
  m_logInts = mm[,use_cols]
  # replace 0's with NA as more appropriate for analysis
  m_logInts[m_logInts==0] = NA
  m_logInts = log2(m_logInts)
  return(m_logInts)
}


# function compute_missing
# computes the nummber of missing and percent missing observations
# PARAMETER
#   mm is a matrix returned by convert_log2 or matrix with NA's representing
#   missing values
# PRINTS out % missing
# RETURN
#   dataframe with 3 columns (num_missing, num_total, perc_missing)
#   and 1 row of values
compute_missing = function(mm) {
  dims = dim(mm)
  nummiss = sum(is.na(mm))
  total = dims[1] * dims[2]
  perc_miss = nummiss / total
  print(paste('Percent missing observations:', perc_miss) )
  ret = data.frame(nummiss, total, perc_miss)
  colnames(ret) = c('num_missing', 'num_total', 'perc_missing')
}

#' Volcano plot
#'
#' Function plots fold changes and p-values as a volcano plot.
#' Two lines are plotted for the p-value cutoff at p = PV_cutoff (solid line)
#' and p = 0.1 (dashed line).
#' @param FC vector of fold changes
#' @param PV vctor of p-values, same lenght as FC
#' @param FC_cutoff fold change cutoff where to draw
#'                  vertical cutoff lines, default = 2
#' @param PV_cutoff p-value cutoff where to draw a
#'                  horisontal cutoff line, default ==.05
#' @param figtitle title to display at the top of the figure, default = ''
#' @return NULL
#' @examples
#' data(mm_peptides)
#' head(mm_peptides)
#' intsCols = 8:13 # different from parameter names as
#'                 # R uses outer name spaces if variable is undefined
#' metaCols = 1:7
#' m_logInts = make_intencities(mm_peptides, intsCols)
#' m_prot.info = make_meta(mm_peptides, metaCols)
#' m_logInts = convert_log2(m_logInts)
#'
#' # Normalize data
#' grps = as.factor(c('CG','CG','CG', 'mCG','mCG','mCG'))
#' mm_m_ints_eig1 = eig_norm1(m=m_logInts,treatment=grps,prot.info=m_prot.info)
#' mm_m_ints_eig1$h.c # check the number of bias trends detected
#'
#' # Impute missing values
#' mm_m_ints_norm = eig_norm2(rv=mm_m_ints_eig1)
#' mm_prot.info = mm_m_ints_norm$normalized[,1:7]
#' mm_norm_m =  mm_m_ints_norm$normalized[,8:13]
#' imp_mm = MBimpute(mm_norm_m, grps, prot.info=mm_prot.info,
#'                   pr_ppos=2, my.pi=0.05,
#'                   compute_pi=FALSE, sseed=131)
#' DE_res = peptideLevel_DE(imp_mm$y_imputed, grps, imp_mm$imp_prot.info,
#'                          pr_ppos=2)
#' plot_volcano(DE_res$FC, DE_res$BH_P_val, FC_cutoff=1.5,
#'              PV_cutoff=.05, figtitle='Mouse DE')
#' @return Nil
#' @export
plot_volcano = function(FC, PV, FC_cutoff=2, PV_cutoff=.05, figtitle='') {
  tmp_x = FC
  tt = PV
  # exact 0 values, replace with a small values as need to take a log
  ppos_rep = tt == 0
  tt[ppos_rep] = .001 # highly significant
  tmp_y = -log10(tt)
  ppos = tmp_y > 20
  num_replace = sum(ppos)
  print(paste("number to replace",num_replace) )
  if(num_replace) {
    tmp_y[ppos] = jitter(rep(19.5, times=num_replace) )
  }
  # graphics::par(mar=c(3,3,3,3))
  xlimits = max(abs(FC)) + .2
  ylimits = max(tmp_y) + .2
  graphics::par(mfcol=c(1,1))
  plot(tmp_x, tmp_y, pch=20, xlim=c(-xlimits,xlimits), ylim=c(0,ylimits),
       # add label to significantly different proteins
  xlab='FC', ylab='-log10(p-value)',  main=figtitle)
  graphics::lines(c(-(xlimits+2),(xlimits+2)), c(-log10(PV_cutoff),
                                       -log10(PV_cutoff)), col='blue')
  # also same as
  graphics::abline(h = -log10(.1), col='blue', lty=2) # permanent line for now
  graphics::lines(c(-FC_cutoff,-FC_cutoff), c(-10, 70), col='blue')
  graphics::lines(c(FC_cutoff,FC_cutoff), c(-10, 70), col='blue')
}


#' Volcano plot with labels for the differentially expressed proteins
#'
#' Function plots fold changes and p-values as a volcano plot.
#' Two lines are plotted for the p-value cutoff at p = PV_cutoff
#' (solid line) and p = 0.1 (dashed line).
#' @param FC vector of fold changes
#' @param PV vector of p-values, same lenght as FC
#' @param ProtID vector of protein IDs, can be gene IDs, same lenght as FC & PV.
#'        Namaes in this vector will be displayed in the volcano plot
#'        for differentially expressed proteins for this reason short names
#'        are preferred.
#' @param FC_cutoff fold change cutoff where to draw vertical cutoff
#'                lines, default = 2
#' @param PV_cutoff p-value cutoff where to draw a horisontal cutoff line,
#'                default ==.05
#' @param figtitle title to display at the top of the figure, default = ''
#' @return NULL
#' @examples
#' data(mm_peptides)
#' head(mm_peptides)
#' intsCols = 8:13 # different from parameter names as
#'                 # R uses outer name spaces if variable is undefined
#' metaCols = 1:7 # reusing this variable
#' m_logInts = make_intencities(mm_peptides, intsCols)  # will reuse the name
#' m_prot.info = make_meta(mm_peptides, metaCols)
#' m_logInts = convert_log2(m_logInts)
#'
#' # Normalize data
#' grps = as.factor(c('CG','CG','CG', 'mCG','mCG','mCG'))
#' mm_m_ints_eig1 = eig_norm1(m=m_logInts,treatment=grps,prot.info=m_prot.info)
#' mm_m_ints_eig1$h.c # check the number of bias trends detected
#'
#' # Impute missing values
#' mm_m_ints_norm = eig_norm2(rv=mm_m_ints_eig1)
#' mm_prot.info = mm_m_ints_norm$normalized[,1:7]
#' mm_norm_m =  mm_m_ints_norm$normalized[,8:13]
#' imp_mm = MBimpute(mm_norm_m, grps, prot.info=mm_prot.info,
#'                   pr_ppos=2, my.pi=0.05,
#'                   compute_pi=FALSE, sseed=131)
#' DE_res = peptideLevel_DE(imp_mm$y_imputed, grps, imp_mm$imp_prot.info,
#'                          pr_ppos=2)
#' plot_volcano_wLab(DE_res$FC, DE_res$BH_P_val, DE_res$ProtID, FC_cutoff=1.5,
#'                   PV_cutoff=.05, figtitle='Mouse DE')
#'
#' @import ggrepel ggplot2
#' @importFrom dplyr mutate
#' @return Nil
#' @export
plot_volcano_wLab = function(FC, PV, ProtID,
                             FC_cutoff=2,
                             PV_cutoff=.05, figtitle='') {
  plotdata = data.frame(FC, PV, ProtID) # combine into 2 data frame
  ppos_rep = plotdata$PV == 0
  plotdata$PV[ppos_rep] = .000000001
  plotdata = dplyr::mutate(plotdata, log_PV=-log10(PV))
  plotdata$threshold = as.factor(abs(plotdata$FC) >= FC_cutoff
                                 & plotdata$PV < PV_cutoff)
  dim(plotdata)

  ggplot() + geom_point(data=plotdata, aes(x=FC, y=log_PV), alpha=0.5, size=1) +
    theme(legend.position = "none") +
    xlab("log2 fold change") +
    ylab("-log10 p-value") +
    geom_text_repel(data=filter(plotdata, threshold==TRUE),
                    size = 3, alpha=.8, aes(x=FC, y=log_PV, label=ProtID) ) +
    theme_classic(base_size = 8) +
    geom_hline(yintercept=-log10(PV_cutoff), col='blue', alpha=.7) +
    geom_hline(yintercept=-log10(0.1),
               col='blue', linetype="dashed", alpha=.7) +
    geom_vline(xintercept=-FC_cutoff, col='blue', alpha=.7) +
      geom_vline(xintercept=FC_cutoff, col='blue', alpha=.7)
}

# DESCRIPTION
# This function calculates multi-sample two-part tests to test for a difference
# in means and proportions between two groups across multiple biological samples
# P-values are generated using a permutation null distribution.
#
# USAGE
# See TwoPartMultiMatrix_Example.r
#
# ARGUMENTS
# data - a list with each element one biological p x n matrix where p is the
#        the number of compounds and n is the number of subjects
# groups - a vector indicating group membership. Groups must be coded 0,1.
#          Only 2 group comparisons are supported with this code.
# n.perm - the number of permutations to generate for the null distribution
# set.seed - sets seed for permutations

# VALUE
# Returns a data.frame with the of T and T-Squared statistics for each compound
# and permutation p values.
# peptideLevel_DE = function(mm, treatment, prot.info, pr_ppos=2)
  # Imputes missing values based on information from multiple peptides
  # within a protein
  # Input:
  #   mm:    m x n matrix of intensities, numpeptides x numsamples
  #   treatment:  vector indicating the treatment group of
  #               each sample ie [1 1 1 1 2 2 2 2...]
  #   pr_ppos - column index for protein ID in prot.info.
  #     Can restrict to be #2...
  #
  # Output: matrix with 4 columns: protID, FC, p-value
  #   ProtID - taken from prot.info (col 2 usually)
  #   FC (fold change)
  #   p-value for comparison between 2 groups (2 groups only here)
  #      as interested in pariwise differences.
  #   BH-adjusted p-value, Benjamini-Hochberg


#' Multi-Matrix Differentia Expression Analysis
#'
#' Multi-Matrix Differential Expression Analysis computes Model-Based
#' statistics for each dataset, the sum of individual statistics is the
#' final statistic. The significance is determined via a permutation test
#' which computed the same statistics and sums them after permuting
#' the values across treatment groups. As is outlined in Karpievitch
#' et al. 2018.
#'
#' @param mm_list list of matrices for each experiment, length = number of
#'                datasets to compare internal dataset dimentions: numpeptides
#'                x numsamples for each dataset
#' @param treat list of data frames with treatment information to compute
#'              the statistic
#'              in same order as mm_list
#' @param prot.info list of protein and peptide mapping for each matrix
#'              in mm_list,
#'              in same order as mm_list
#' @param prot_col_name column name in prot.info that contains protein
#'              identifiers that link all datasets together. Not that
#'              Protein IDs will differ across different organizms and
#'              cannot be used as the linking identifier.
#'              Function match_linker_ids() produces numeric identifyers
#'              that link all datasets together
#' @param nperm number of permutations, default = 500,
#'              this will take a while, test code
#'              with fewer permutations
#' @param setseed random number generator seed, default = 12345
#' @param dataset_suffix vector of character strings that corresponds to the
#'        dataset being analysed. Same length as mm_list. Names will be appended
#'        to the columns names that will be generated for each analysed dataset.
#'        For example, if analysing mouse and human data this vector may be:
#'        c('Mouse', 'Human')
#' @return data frame with the following columns
#' \describe{
#'   \item{protIDused}{Column containing the protien IDs used to
#'                     link proteins across datasets}
#'   \item{FC}{Average fold change across all datasets}
#'   \item{P_val}{Permutation-based p-valu for the differences
#'                between the groups}
#'   \item{BH_P_val}{Multiple testing adjusted p-values}
#'   \item{statistic}{Statistic computed as a a sum of statistics
#'                    produced for each dataset}
#'   \item{Protein Information}{all columns passed into the function
#'                              for the 1st dataset
#'         in the list}
#'   \item{FCs}{Fold changes for individual datasets, these values
#'              should average to the
#'         FC above. As many columns as there are datasets being analyzed.}
#'   \item{PV}{p-values for individual datasets. As many
#'             columns as there are datasets
#'         being analyzed.}
#'   \item{BHPV}{Multiple testing adjusted p-values for
#'               individual datasets. As many
#'         columns as there are datasets being analyzed.}
#'   \item{NUMPEP}{Number of peptides presents in each protien
#'                 for each dataset. As many
#'         columns as there are datasets being analyzed.}
#'}
#' @examples
#' # Load mouse dataset
#' data(mm_peptides)
#' head(mm_peptides)
#' intsCols = 8:13 # different from parameter names as R uses
#'                 # outer name spaces if variable is undefined
#' metaCols = 1:7 # reusing this variable
#' m_logInts = make_intencities(mm_peptides, intsCols)  # will reuse the name
#' m_prot.info = make_meta(mm_peptides, metaCols)
#' m_logInts = convert_log2(m_logInts)
#' grps = as.factor(c('CG','CG','CG', 'mCG','mCG','mCG'))
#' mm_m_ints_eig1 = eig_norm1(m=m_logInts,treatment=grps,
#'                            prot.info=m_prot.info)
#' mm_m_ints_eig1$h.c # check the number of bias trends detected
#' mm_m_ints_norm = eig_norm2(rv=mm_m_ints_eig1)
#' mm_prot.info = mm_m_ints_norm$normalized[,1:7]
#' mm_norm_m =  mm_m_ints_norm$normalized[,8:13]
#' imp_mm = MBimpute(mm_norm_m, grps, prot.info=mm_prot.info,
#'                   pr_ppos=2, my.pi=0.05,
#'                   compute_pi=FALSE, sseed=131)
#'
#' # Load human dataset
#' data(hs_peptides)
#' head(hs_peptides)
#' intsCols = 8:13 # different from parameter names as R uses
#'                 # outer name spaces if variable is undefined
#' metaCols = 1:7 # reusing this variable
#' m_logInts = make_intencities(hs_peptides, intsCols)  # will reuse the name
#' m_prot.info = make_meta(hs_peptides, metaCols)
#' m_logInts = convert_log2(m_logInts)
#' grps = as.factor(c('CG','CG','CG', 'mCG','mCG','mCG'))
#' hs_m_ints_eig1 = eig_norm1(m=m_logInts,treatment=grps,prot.info=m_prot.info)
#' hs_m_ints_eig1$h.c # check the number of bias trends detected
#' hs_m_ints_norm = eig_norm2(rv=hs_m_ints_eig1)
#' hs_prot.info = hs_m_ints_norm$normalized[,1:7]
#' hs_norm_m =  hs_m_ints_norm$normalized[,8:13]
#' imp_hs = MBimpute(hs_norm_m, grps, prot.info=hs_prot.info,
#'                   pr_ppos=2, my.pi=0.05,
#'                   compute_pi=FALSE, sseed=131)
#'
#' # Multi-Matrix Model-based differential expression analysis
#' # Set up needed variables
#' mms = list()
#' treats = list()
#' protinfos = list()
#' mms[[1]] = imp_mm$y_imputed
#' mms[[2]] = imp_hs$y_imputed
#' treats[[1]] = grps
#' treats[[2]] = grps
#' protinfos[[1]] = imp_mm$imp_prot.info
#' protinfos[[2]] = imp_hs$imp_prot.info
#' nperm = 50
#' comb_MBDE = prot_level_multi_part(mm_list=mms, treat=treats,
#'                                   prot.info=protinfos,
#'                                   prot_col_name='ProtID', nperm=nperm,
#'                                   setseed=123, dataset_suffix=c('MM', 'HS'))
#'
#' # Analysis for proteins only present in mouse,
#' # there are no proteins suitable for
#' # Model-Based analysis in human dataset
#' subset_data = subset_proteins(mm_list=mms, prot.info=protinfos, 'MatchedID')
#' mm_dd_only = subset_data$sub_unique_mm_list[[1]]
#' hs_dd_only = subset_data$sub_unique_mm_list[[2]]
#' protinfos_mm_dd = subset_data$sub_unique_prot.info[[1]]
#' DE_mCG_CG_mm_dd = peptideLevel_DE(mm_dd_only, grps,
#'                                   prot.info=protinfos_mm_dd, pr_ppos=2)
#'
#' @export
prot_level_multi_part = function(mm_list, treat, prot.info,
                                  prot_col_name, nperm=500, setseed=12345,
                                  dataset_suffix){
  # select proteins that were detected in each experiment
  # make a list of unique protein IDs for each matrix in the list mm_list
  # grps will not change
  subset_data = subset_proteins(mm_list=mm_list, prot.info=prot.info,
                                prot_col_name)
  # names(subset_data)
  # "sub_mm_list"  "sub_prot.info"  "sub_unique_mm_list"
  # "sub_unique_prot.info"  "common_list"

  print('Computing statistics')
  sub_mm_list = subset_data$sub_mm_list
  sub_prot.info = subset_data$sub_prot.info
  nsets = length(mm_list)

  tt = colnames(sub_prot.info[[1]])
  ttt = tt == prot_col_name
  # should be position of the column that was passed in for ID
  pos_prot_id_col = seq(1,length(tt))[ttt]

  ## Tstat = CalcT(data, groups)  # my test here return t-stat
  # for each dataset loop through, compute, and add up t-stat values
  tmp = peptideLevel_DE(sub_mm_list[[1]], treat[[1]], sub_prot.info[[1]],
                        pr_ppos=pos_prot_id_col)
  tstat_all = list()
  tstat_all[[1]] = tmp
  # t_value si stored in col 5:  ProtID  FC  p-val
  # BH_p-val  t_value   num_peptides
  tstat = as.double(tmp[,5])
  FCs = as.double(tmp[,2])
  PV = as.double(tmp[,3])
  BHPV = as.double(tmp[,4])
  NUMPEP = as.numeric(tmp[,6])
  col_FC = paste('FC_', dataset_suffix[1], sep='')
  col_PV = paste('PV_', dataset_suffix[1], sep='')
  col_BHPV = paste('BHPV_', dataset_suffix[1], sep='')
  col_NUMPEP = paste('NUMPEP_', dataset_suffix[1], sep='')
  # prot names will be the same, will not combine them in the loop
  PROTIDS = tmp[,1]
  for(ii in 2:nsets){ # for second and more datasets
    tmp = peptideLevel_DE(sub_mm_list[[ii]], treat[[ii]],
                          sub_prot.info[[ii]], pr_ppos=pos_prot_id_col)
    tstat_all[[ii]] = tmp
    # yuliya: may need to subset here, tmp is complex var
    tstat = cbind(tstat, as.double(tmp[,5]))
    FCs = cbind(FCs, tmp[,2])
    PV =  cbind(PV, tmp[,3])
    BHPV = cbind(BHPV, tmp[,4])
    NUMPEP = cbind(NUMPEP, tmp[,6])

    # column headers
    col_FC = c(col_FC, paste('FC_', dataset_suffix[ii], sep=''))
    col_PV = c(col_PV, paste('PV_', dataset_suffix[ii], sep=''))
    col_BHPV = c(col_BHPV, paste('BHPV_', dataset_suffix[ii], sep=''))
    col_NUMPEP = c(col_NUMPEP, paste('NUMPEP_', dataset_suffix[ii], sep=''))
  }
  colnames(FCs) = col_FC
  colnames(PV) = col_PV
  colnames(BHPV) = col_BHPV
  colnames(NUMPEP) = col_NUMPEP

  sum_tstat = rowSums(tstat)

  print('Perfoming permutation test')

  tstat_perm = list()
  for(ii in seq(1,nsets)) {
    print(paste('Dataset', as.character(ii) ) )
    tstat_perm[[ii]] = NULL
    for(jj in seq(1,nperm)) {
      # get permuted labels for each iteration, then compute T_p
      perm_pos = sample(length(treat[[ii]]), length(treat[[ii]]) )
      tmp = peptideLevel_DE(sub_mm_list[[ii]],
                            treat[[ii]][perm_pos], sub_prot.info[[ii]],
                            pr_ppos=pos_prot_id_col)
      if(jj == 1) {
        tstat_perm[[ii]] =  as.matrix(as.double(tmp[,5]))
      } else {
        tstat_perm[[ii]] = cbind(tstat_perm[[ii]],as.matrix(as.double(tmp[,5])))
      }
    }
  }

  # sum the matrices
  T_perm = tstat_perm[[1]]
  for(ii in 2:nsets) {
    T_perm = T_perm + tstat_perm[[ii]]
  }

  num_prot = dim(tstat)[1]
  p_vals = vector('numeric', length=num_prot)
  pos_stat_pos = sum_tstat >= 0
  for(ii in seq(1,2)) { # positive and negative values separately
    if(ii == 1) {
      ppos = which(pos_stat_pos)
      for(kk in seq(1,length(ppos))) {
        p_vals[ppos[kk]] = (.5+sum(T_perm[ppos[kk],] >=
                                     sum_tstat[ppos[kk]])) / (nperm+1)
      }
    } else {
      ppos = which(!pos_stat_pos)
      for(kk in seq(1,length(ppos))) {
        p_vals[ppos[kk]] = (.5+ sum(T_perm[ppos[kk],] <
                                      sum_tstat[ppos[kk]])) / (nperm+1)
      }
    }
  }
 # It seesm that p-values produced by the permitation
 # test are [0, .6], thus standard
 # adjustments will not do very well.
 # I will use 'fdr' orption in p.adjust and then rescale the interval [0 1].
 # p-values look the best, according to te theoretical
 # distribution, after such adjustment
 p_vals_tmp =  p.adjust(p_vals, method="fdr")
 mmin = min(p_vals_tmp)
 mmax = max(p_vals_tmp)
 adj_PV = (p_vals_tmp - mmin) / (mmax-mmin)
 FC = rowMeans(FCs)

 # protein info is on peptide level, so convert to
 # the protein level info, no duplication
 # take prot IDs from dataset 1
 unik = !duplicated(sub_prot.info[[1]][,prot_col_name])
 ppos_u_prots = seq_along(sub_prot.info[[1]][,prot_col_name])[unik]  ## indices
 u_prot_info = sub_prot.info[[1]][ppos_u_prots,]

 res = data.frame(protIDused=PROTIDS, FC, P_val=p_vals,
                  BH_P_val=adj_PV, statistic=sum_tstat,
                  u_prot_info, FCs, PV, BHPV, NUMPEP)
 # column names in res are inherited from the structures that are combined
 return(res)
}


#' Presence/Absence peptide-level analysis
#'
#' Presence/Absence peptide-level analysis uses
#' all peptides for a protein as IID
#' to produce 1 p-value across multiple (2+) datasets.
#' Significance is estimated using a g-test which is suitable
#' for two treatment groups only.
#'
#' @param mm m x n matrix of intensities, number of peptides x number of samples
#' @param treatment vector indicating the treatment
#'                  group of each sample ie [1 1 1 1 2 2 2 2...]
#' @param prot.info 2+ colum data frame of peptide ID, protein ID, etc. columns
#' @param pr_ppos - column index for protein ID in
#'               prot.info. Can restrict to be #2...
#'
#' @return A list of length two items:
#' \describe{
#'   \item{ProtIDused}{protein identification information taken from prot.info,
#'         a column used to identify proteins}
#'   \item{FC}{Approximation of the fold change computed
#'        as percent missing observations
#'        group 1 munis in percent missing observations group 2}
#'   \item{P_val}{p-value for the comparison between
#'                2 groups (2 groups only here)}
#'   \item{BH_P_val}{Benjamini-Hochberg adjusted p-values}
#'   \item{statistic}{statistic returned by
#'          the g-test, not very useful as depends on
#'         the direction of the test and can produce all 0's}
#'   \item{num_peptides}{number of peptides within a protein}
#'   \item{metadata}{all columns of metadata from teh matrix that was passed in}
#'}
#' @examples
#' # Load mouse dataset
#' data(mm_peptides)
#' head(mm_peptides)
#' intsCols = 8:13 # different from parameter names as R uses
#'                 # outer name spaces if variable is undefined
#' metaCols = 1:7 # reusing this variable
#' m_logInts = make_intencities(mm_peptides, intsCols)  # will reuse the name
#' m_prot.info = make_meta(mm_peptides, metaCols)
#' m_logInts = convert_log2(m_logInts)
#' grps = as.factor(c('CG','CG','CG', 'mCG','mCG','mCG'))
#' mm_m_ints_eig1 = eig_norm1(m=m_logInts,treatment=grps,prot.info=m_prot.info)
#' mm_m_ints_eig1$h.c # check the number of bias trends detected
#' mm_m_ints_norm = eig_norm2(rv=mm_m_ints_eig1)
#'
#' # Load human dataset
#' data(hs_peptides)
#' head(hs_peptides)
#' intsCols = 8:13 # different from parameter names as R
#'                 # uses outer name spaces if variable is undefined
#' metaCols = 1:7 # reusing this variable
#' m_logInts = make_intencities(hs_peptides, intsCols)  # will reuse the name
#' m_prot.info = make_meta(hs_peptides, metaCols)
#' m_logInts = convert_log2(m_logInts)
#' grps = as.factor(c('CG','CG','CG', 'mCG','mCG','mCG'))
#' hs_m_ints_eig1 = eig_norm1(m=m_logInts,treatment=grps,prot.info=m_prot.info)
#' hs_m_ints_eig1$h.c # check the number of bias trends detected
#' hs_m_ints_norm = eig_norm2(rv=hs_m_ints_eig1)
#'
#' # Set up for presence/absence analysis
#' raw_list = list()
#' norm_imp_prot.info_list = list()
#' raw_list[[1]] = mm_m_ints_eig1$m
#' raw_list[[2]] = hs_m_ints_eig1$m
#' norm_imp_prot.info_list[[1]] = mm_m_ints_eig1$prot.info
#' norm_imp_prot.info_list[[2]] = hs_m_ints_eig1$prot.info
#'
#' protnames_norm_list = list()
#' protnames_norm_list[[1]] = unique(mm_m_ints_norm$normalized$MatchedID)
#' protnames_norm_list[[2]] = unique(hs_m_ints_norm$normalized$MatchedID)
#'
#' presAbs_dd = get_presAbs_prots(mm_list=raw_list,
#'                               prot.info=norm_imp_prot.info_list,
#'                               protnames_norm=protnames_norm_list,
#'                               prot_col_name=2)
#'
#' presAbs_de = peptideLevel_PresAbsDE(presAbs_dd[[1]][[1]],
#'                                     grps, presAbs_dd[[2]][[1]],
#'                                     pr_ppos=2)
#' @export
peptideLevel_PresAbsDE = function(mm, treatment, prot.info, pr_ppos=2){
  #   XsqDE (degrees of freedom) -- not used, place holder for where
  #   FC goes in abundace-based DE
  #   p-value for comparison between 2 groups (2 groups only here)
  #      as interested in pariwise differences.
  #   BH-adjusted p-value, Benjamini-Hochberg multiple testing adjustment

  # Match to protein
  all.proteins = unique(prot.info[,pr_ppos])
  numProts = length(all.proteins) # 1569
  y_out = data.frame(matrix(NA, numProts, 5))
  nummiss = data.frame(matrix(NA, numProts, 2))
  u_treat = unique(treatment)
  numgrps = length(u_treat)
  numeeachgroup =  vector('numeric', length=numgrps)
  for(ii in seq(1,numgrps)) {
    numeeachgroup[ii] = sum(treatment == u_treat[ii])
  } # needed for the FC estimation

  de_ret = NULL
  u_prot_info = NULL
  for (kk in seq(1,length(all.proteins))) {
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
    nn = length(yy)
    # keep track of prIDs here...
    curr_prot.info = curr_prot.info[kk,] # possibly a subset
    # good to know how many peptides were in a given protein
    y_out[kk,5] = n.peptide

    # replicate treatment for # peptides
    treatment_hold = treatment
    treatment = rep(treatment, times=n.peptide)
    # use g-test to compute stat significance in differences between groups
    xx = is.na(yy)

    treatsX = unique(treatment)
    if(sum(xx) < nn & sum(xx) > 0) {
      res = g.test(xx,treatment)
      y_out[kk,2] = res$p.value
      y_out[kk,1] = res$parameter
      y_out[kk,4] = res$statistic
    } else {
      y_out[kk,2] = 2
      # all values are OUT of RANGE of proper results, will be
      # used to flag (and remove) all preent or all absent values
      y_out[kk,1] = 2
      y_out[kk,4] = 2 # more likely to have all NAs
    }
    # count # of missing values in each tretment,
    # will vary depending on the number of peptides
    nummiss[kk,1] = sum(xx[treatment==treatsX[1]]==TRUE)
    nummiss[kk,2] = sum(xx[treatment==treatsX[2]]==TRUE)
    treatment = treatment_hold
  } # end for each protein

  colnames(y_out) = c('FC', 'P_val', 'BH_P_val', 'statistic', 'num_peptides')
  # BH adjustment - only values on interval [0 1]
  # we may have 2 if g.test cannot be performed
  ppos = y_out[,2] <= 1
  y_out[ppos,3] = stats::p.adjust(y_out[ppos,2],"BH")
  y_out[!ppos,3] = 1 # these would have been 2 in Raw p-values

  # XsqDF returned from the g-test is not useful,
  # and is the same as statistic here
  # so calculate my estimate of fold change as:
  # (% pep missing grp1) / (% pep missing grp1)

  # HERE make a dataframe to be returned -
  # add protein names as 1st col in a data frame
  DE_res = data.frame(all.proteins, y_out, stringsAsFactors=FALSE)
  de_ret$DE_res = DE_res
  de_ret$prot.info = u_prot_info

  num_obs =  matrix(0, length(all.proteins), numgrps)
  for(ii in seq(1,numgrps)) {
    num_obs[,ii] = de_ret$DE_res$num_peptides * numeeachgroup[ii]
  }
  percmiss = nummiss / num_obs
  de_ret$DE_res$FC = percmiss[,1] - percmiss[,2]
  cols1 = colnames(de_ret$DE_res)
  cols1[1] = "ProtIDused"
  cols2 = colnames(de_ret$prot.info)
  de_ret = data.frame(de_ret, stringsAsFactors = FALSE)
  colnames(de_ret) = c(cols1, cols2)
  return(de_ret)
}

#########################################################################
#' Multi-Matrix Presence Absence analysis
#'
#' Multi-Matrix Presence Absence Analysis computes Model-Based
#' statistics for each dataset and sums them up to produce the final
#' statistic. The significance is determined via a permutation
#' test which computes the same statistics and sums them
#' after permuting the values across treatment groups,
#' as is outlined in Karpievitch et al. 2018. Whenever possible
#' proteins should be analysed using the Model-Based
#' Differential Expression Analysis due to higher statistical
#' power over the Presence Absence analysis.
#'
#' @param mm_list list of matrices of intensities for each experiment,
#'        dimentions: numpeptides x numsamples
#' @param treat list of data frames with treatment information to
#'        compute the statistic,
#'        parallel to mm_list and prot.info
#' @param prot.info list of protein metadata for each matrix in
#'        mm_list, data.frame
#'        parallel to mm_list and treat
#' @param prot_col_name column names present in all datasets that
#'        identifies protein IDs
#'        across all datasets
#' @param nperm number of permutations
#' @param setseed random number generator seed, is used to
#'                permute the data, default 13457
#' @param dataset_suffix a list of strings that will be
#'         appended to the column names
#'         for FC, PV, BHPV and numebers of peptides
#'
#' @return a data frame with the following columns:
#' \describe{
#'   \item{protIDused}{protein metadata, peptide sequence if was
#'        passed in as one of the columns is the first peptide
#'        equence encountered in the data for that protein}
#'   \item{FCs}{Avegares across all datasets of the approximation
#'        of the fold change computed as percent missing observations
#'        group 1 munis in percent missing observations group 2 in
#'        peptideLevel_PresAbsDE() function}
#'   \item{P_val}{p-value for the comparison between 2 groups
#'               (2 groups only here) obtained from a permutation test}
#'   \item{BH_P_val}{Benjamini-Hochberg adjusted p-values}
#'   \item{statistic}{statistic returned by the g-test and
#'        summed across all datasets,
#'        not very useful as depends on the direction of the
#'        test and can produce all 0's}
#'   \item{u_prot_info}{column containing ptoein identifiers
#'                      across all datasets}
#'   \item{FCs}{Approximation of the fold change computed
#'        as percent missing observations
#'        group 1 munis in percent missing observations
#'        group 2 in peptideLevel_PresAbsDE() function}
#'   \item{PV}{p-values produced by g-test for individual datasets}
#'   \item{BHPV}{adjusted p-values produced by g-test for individual datasets}
#'   \item{NUMPEP}{number of peptides observed for
#'         each protein in each of the datasets}
#' }
#' @export
#' @examples
#' # Load mouse dataset
#' data(mm_peptides)
#' head(mm_peptides)
#' intsCols = 8:13
#' metaCols = 1:7
#' m_logInts = make_intencities(mm_peptides, intsCols)  # will reuse the name
#' m_prot.info = make_meta(mm_peptides, metaCols)
#' m_logInts = convert_log2(m_logInts)
#' grps = as.factor(c('CG','CG','CG', 'mCG','mCG','mCG'))
#' mm_m_ints_eig1 = eig_norm1(m=m_logInts,treatment=grps,prot.info=m_prot.info)
#' mm_m_ints_eig1$h.c # check the number of bias trends detected
#' mm_m_ints_norm = eig_norm2(rv=mm_m_ints_eig1)
#'
#' # Load human dataset
#' data(hs_peptides)
#' head(hs_peptides)
#' intsCols = 8:13
#' metaCols = 1:7
#' m_logInts = make_intencities(hs_peptides, intsCols)
#' m_prot.info = make_meta(hs_peptides, metaCols)
#' m_logInts = convert_log2(m_logInts)
#' grps = as.factor(c('CG','CG','CG', 'mCG','mCG','mCG'))
#' hs_m_ints_eig1 = eig_norm1(m=m_logInts,treatment=grps,prot.info=m_prot.info)
#' hs_m_ints_eig1$h.c # check the number of bias trends detected
#' hs_m_ints_norm = eig_norm2(rv=hs_m_ints_eig1)
#'
#' # Set up for presence/absence analysis
#' raw_list = list()
#' norm_imp_prot.info_list = list()
#' raw_list[[1]] = mm_m_ints_eig1$m
#' raw_list[[2]] = hs_m_ints_eig1$m
#' norm_imp_prot.info_list[[1]] = mm_m_ints_eig1$prot.info
#' norm_imp_prot.info_list[[2]] = hs_m_ints_eig1$prot.info
#'
#' protnames_norm_list = list()
#' protnames_norm_list[[1]] = unique(mm_m_ints_norm$normalized$MatchedID)
#' protnames_norm_list[[2]] = unique(hs_m_ints_norm$normalized$MatchedID)
#'
#' presAbs_dd = get_presAbs_prots(mm_list=raw_list,
#'                               prot.info=norm_imp_prot.info_list,
#'                               protnames_norm=protnames_norm_list,
#'                               prot_col_name=2)
#'
#' ints_presAbs = list()
#' protmeta_presAbs = list()
#' ints_presAbs[[1]] = presAbs_dd[[1]][[1]] # Mouse
#' ints_presAbs[[2]] = presAbs_dd[[1]][[2]] # HS
#' protmeta_presAbs[[1]] = presAbs_dd[[2]][[1]]
#' protmeta_presAbs[[2]] = presAbs_dd[[2]][[2]]
#'
#' treats = list()
#' treats[[1]] = grps
#' treats[[2]] = grps
#'
#' subset_presAbs = subset_proteins(mm_list=ints_presAbs,
#'                         prot.info=protmeta_presAbs, 'MatchedID')
#'
#' nperm = 50  # set to 500+ for publication
#' presAbs_comb = prot_level_multiMat_PresAbs(
#'                            mm_list=subset_presAbs$sub_mm_list,
#'                            treat=treats,
#'                            prot.info=subset_presAbs$sub_prot.info,
#'                            prot_col_name='MatchedID', nperm=nperm,
#'                            setseed=123372, dataset_suffix=c('MM', 'HS') )
#'
#' plot_volcano(presAbs_comb$FC, presAbs_comb$BH_P_val,
#'              FC_cutoff=.5, PV_cutoff=.05,
#'              'Combined Pres/Abs CG vs mCG')
#'
prot_level_multiMat_PresAbs = function(mm_list, treat, prot.info, prot_col_name,
                                       nperm=500,
                                       setseed=13457, dataset_suffix){
  # select proteins that were detected in each experiment
  # make a list of unique protein IDs for each matrix in the list mm_list
  subset_data = subset_proteins(mm_list=mm_list,
                                prot.info=prot.info,
                                prot_col_name)  # grps will not change
  # names(subset_data)
  # "sub_mm_list"  "sub_prot.info"
  # "sub_unique_mm_list"  "sub_unique_prot.info"  "common_list"

  print('Computing statistics')
  sub_mm_list = subset_data$sub_mm_list
  sub_prot.info = subset_data$sub_prot.info
  nsets = length(mm_list)

  tt = colnames(sub_prot.info[[1]])
  ttt = tt == prot_col_name
  # should be position of the column that was passed in for ID
  pos_prot_id_col = seq(1,length(tt))[ttt]
  tmp = peptideLevel_PresAbsDE(sub_mm_list[[1]],
                               treat[[1]], sub_prot.info[[1]],
                               pr_ppos=pos_prot_id_col)
  tstat_all = list()
  tstat_all[[1]] = tmp
  # t_value si stored in col 5:  ProtID  FC  p-val
  # BH_p-val  t_value   num_peptides
  tstat = as.double(tmp[,5])
  FCs = as.double(tmp[,2])
  PV = as.double(tmp[,3])
  BHPV = as.double(tmp[,4])
  NUMPEP = as.numeric(tmp[,6])

  col_FC = paste('FC_', dataset_suffix[1], sep='')
  col_PV = paste('PV_', dataset_suffix[1], sep='')
  col_BHPV = paste('BHPV_', dataset_suffix[1], sep='')
  col_NUMPEP = paste('NUMPEP_', dataset_suffix[1], sep='')
  # prot names will be the same, will not combine them in the loop
  PROTIDS = tmp[,1] # ??
  # prot names will be the same, will not combine them in the loop
  for(ii in 2:nsets){ # for second and more datasets
    tmp = peptideLevel_PresAbsDE(sub_mm_list[[ii]],
                                 treat[[ii]],
                                 sub_prot.info[[ii]],
                                 pr_ppos=pos_prot_id_col)
    tstat_all[[ii]] = tmp
    # yuliya: may need to subset here, tmp is complex var
    tstat = cbind(tstat, as.double(tmp[,5]))
    FCs = cbind(FCs, tmp[,2])
    PV =  cbind(PV, tmp[,3])
    BHPV = cbind(BHPV, tmp[,4])
    NUMPEP = cbind(NUMPEP, tmp[,6])

    col_FC = c(col_FC, paste('FC_', dataset_suffix[ii], sep=''))
    col_PV = c(col_PV, paste('PV_', dataset_suffix[ii], sep=''))
    col_BHPV = c(col_BHPV, paste('BHPV_', dataset_suffix[ii], sep=''))
    col_NUMPEP = c(col_NUMPEP, paste('NUMPEP_', dataset_suffix[ii], sep=''))
  }
  colnames(FCs) = col_FC
  colnames(PV) = col_PV
  colnames(BHPV) = col_BHPV
  colnames(NUMPEP) = col_NUMPEP
  sum_tstat = rowSums(tstat) # check that correctly summed over columns

  print('Perfoming permutation test')

  tstat_perm = list()
  for(ii in seq(1,nsets)) {
    print(paste('Dataset', as.character(ii) ) )
    tstat_perm[[ii]] = NULL
    for(jj in seq(1,nperm)) {
      # get permuted labels for each iteration, then compute T_p
      perm_pos = sample(length(treat[[ii]]), length(treat[[ii]]) )
      tmp = peptideLevel_PresAbsDE(sub_mm_list[[ii]],
                                   treat[[ii]][perm_pos],
                                   sub_prot.info[[ii]],
                                   pr_ppos=2)
      if(jj == 1) {
        tstat_perm[[ii]] =  as.matrix(as.double(tmp[,5]))
      } else {
        tstat_perm[[ii]] = cbind(tstat_perm[[ii]],
                                 as.matrix(as.double(tmp[,5])))
      }
    }
  }

  # sum the matrices
  T_perm = tstat_perm[[1]]
  for(ii in 2:nsets) {
    T_perm = T_perm + tstat_perm[[ii]]
  }

  num_prot = dim(tstat)[1]
  p_vals = vector('numeric', length=num_prot)
  pos_stat_pos = sum_tstat >= 0
  for(ii in seq(1,num_prot)) { # positive and negative values separately
    if(ii == 1) {
    p_vals[ii] = (.5+ sum(T_perm[ii,] >= sum_tstat[ii])) / (nperm+1)
      ppos = which(pos_stat_pos)
      for(kk in seq(1,length(ppos))) {
        p_vals[ppos[kk]] = (.5+ sum(T_perm[ppos[kk],] >=
                                      sum_tstat[ppos[kk]])) / (nperm+1)
      }
    } else {
      ppos = which(!pos_stat_pos)
      for(kk in seq(1,length(ppos))) {
        p_vals[ppos[kk]] = (.5+ sum(T_perm[ppos[kk],] <
                                      sum_tstat[ppos[kk]])) / (nperm+1)
      }
    }
  }

  # It seems that sometimes p-values produced by the
  # permitation test are [0, .5],
  # thus standard adjustments may not do very well.
  # One option is to rescale raw p-values to the
  # interval with max of 1.
  # P-values can look the best (in accordance with the
  # theoretical distribution, after such adjustment.
  adj_PV = stats::p.adjust(p_vals, method = 'BH')
  # I will use 'fdr' option in p.adjust and then rescale the interval [0 1].
  # p-values look the best, according to te theoretical
  # distribution, after such adjustment - below works well for ANOVA
  # p_vals_tmp =  stats::p.adjust(p_vals, method="fdr")
  # mmin = min(p_vals_tmp)
  # mmax = max(p_vals_tmp)
  # adj_PV = (p_vals_tmp - mmin) / (mmax-mmin)
  FC = rowMeans(FCs)

  # protein info is on peptide level, so convert to
  # the protein level info, no duplication
  # take prot IDs from dataset 1
  unik = !duplicated(sub_prot.info[[1]][,prot_col_name])
  ppos_u_prots = seq_along(sub_prot.info[[1]][,prot_col_name])[unik]  # indices
  u_prot_info = sub_prot.info[[1]][ppos_u_prots,]

  res = data.frame(protIDused=PROTIDS, FC, P_val=p_vals,
                   BH_P_val=adj_PV, statistic=sum_tstat,
                 u_prot_info, FCs, PV, BHPV, NUMPEP)
  # column names in res are inherited from the structures
  # that are combined into the data frame
  return(res)
}


#' Subset proteins
#'
#' Subset proteins into ones common to all datasets passed
#' into the function and unique to each dataset. Note: for 3+ datasets
#' no intermediate combinations of proteins are returned, only
#' proteins common to all datasets, the rest are
#' returned as unique to each dataset.
#'
#' @param mm_list list of matrices for each experiment,
#'                length = number of datasets to compare
#'                internal dataset dimentions:
#'                numpeptides x numsamples for each dataset
#' @param prot.info list of protein and peptide mapping
#'              for each matrix in mm_list,
#'              in same order as mm_list
#' @param prot_col_name column name in prot.info that contains
#'              protein identifiers that
#'              link all datasets together.
#'              Not that Protein IDs will differ across
#'              different organizms and cannot be used
#'              as the linking identifier.
#'              Function match_linker_ids() produces
#'              numeric identifyers that link all
#'              datasets together
#' @return data frame with the following columns
#' \describe{
#'   \item{sub_mm_list}{list of dataframes of intensities
#'         for each of the datasets
#'         passed in with proteins present in all datasets}
#'   \item{sub_prot.info}{list of dataframes of metadata
#'         for each of the datasets
#'         passed in with proteins present in all datasets.
#'         Same order as sub_mm_list}
#'   \item{sub_unique_mm_list}{list of dataframes of
#'         intensities not found in all
#'         datasets}
#'   \item{sub_unique_prot.info}{ist of dataframes of
#'         metadata not found in all
#'         datasets}
#'   \item{common_list}{list of protein IDs commnon to all datasets}
#' }
#' @examples
#' # Load mouse dataset
#' data(mm_peptides)
#' head(mm_peptides)
#' # different from parameter names as R uses
#' # outer name spaces if variable is undefined
#' intsCols = 8:13
#' metaCols = 1:7 # reusing this variable
#' m_logInts = make_intencities(mm_peptides, intsCols)  # will reuse the name
#' m_prot.info = make_meta(mm_peptides, metaCols)
#' m_logInts = convert_log2(m_logInts)
#' grps = as.factor(c('CG','CG','CG', 'mCG','mCG','mCG'))
#' mm_m_ints_eig1 = eig_norm1(m=m_logInts,treatment=grps,prot.info=m_prot.info)
#' mm_m_ints_eig1$h.c # check the number of bias trends detected
#' mm_m_ints_norm = eig_norm2(rv=mm_m_ints_eig1)
#' mm_prot.info = mm_m_ints_norm$normalized[,1:7]
#' mm_norm_m =  mm_m_ints_norm$normalized[,8:13]
#' imp_mm = MBimpute(mm_norm_m, grps,
#'                   prot.info=mm_prot.info, pr_ppos=2, my.pi=0.05,
#'                   compute_pi=FALSE, sseed=131)
#'
#' # Load human dataset
#' data(hs_peptides)
#' head(hs_peptides)
#' intsCols = 8:13
#' metaCols = 1:7 # reusing this variable
#' m_logInts = make_intencities(hs_peptides, intsCols)  # will reuse the name
#' m_prot.info = make_meta(hs_peptides, metaCols)
#' m_logInts = convert_log2(m_logInts)
#' grps = as.factor(c('CG','CG','CG', 'mCG','mCG','mCG'))
#' hs_m_ints_eig1 = eig_norm1(m=m_logInts,treatment=grps,prot.info=m_prot.info)
#' hs_m_ints_eig1$h.c # check the number of bias trends detected
#' hs_m_ints_norm = eig_norm2(rv=hs_m_ints_eig1)
#' hs_prot.info = hs_m_ints_norm$normalized[,1:7]
#' hs_norm_m =  hs_m_ints_norm$normalized[,8:13]
#' imp_hs = MBimpute(hs_norm_m, grps,
#'                   prot.info=hs_prot.info, pr_ppos=2,
#'                   my.pi=0.05,
#'                   compute_pi=FALSE, sseed=131)
#'
#' # Multi-Matrix Model-based differential expression analysis
#' # Set up needed variables
#' mms = list()
#' treats = list()
#' protinfos = list()
#' mms[[1]] = imp_mm$y_imputed
#' mms[[2]] = imp_hs$y_imputed
#' treats[[1]] = grps
#' treats[[2]] = grps
#' protinfos[[1]] = imp_mm$imp_prot.info
#' protinfos[[2]] = imp_hs$imp_prot.info
#'
#' subset_data = subset_proteins(mm_list=mms, prot.info=protinfos, 'MatchedID')
#' mms_mm_dd = subset_data$sub_unique_mm_list[[1]]
#' protinfos_mm_dd = subset_data$sub_unique_prot.info[[1]]
#' # DIfferential expression analysis for mouse specific protiens
#' DE_mCG_CG_mm_dd = peptideLevel_DE(mms_mm_dd, grps,
#'                                   prot.info=protinfos_mm_dd, pr_ppos=2)
#'
#' @export
subset_proteins = function(mm_list, prot.info, prot_col_name) {
  ll = length(mm_list)
  numuprots = vector('numeric', ll)
  common_list = ''
  uprots = list()
  for(ii in seq(1,ll)) {
    uprots[[ii]] = unique(prot.info[[ii]][,c(prot_col_name)])
    numuprots = length(uprots[[ii]])
    if(ii == 1) {
      common_list = uprots[[ii]]
    } else {
      # match protein names across multiple datasets, keep only overlapping ones
      common_list = intersect(common_list,uprots[[ii]])
    }
  }

  # subset each experiment matrix to the proteins
  # that are in ALL of the datasets/experiment
  # stored in common_list; need to have unique proteins
  # from each dataset, so not doing in the loop above
  sub_mm_list = list()
  sub_prot.info = list()
  sub_unique_mm_list = list()
  sub_unique_prot.info = list()

  for(ii in seq(1,ll)) {
    ppos = prot.info[[ii]][,c(prot_col_name)] %in% common_list
    sub_mm_list[[ii]] = mm_list[[ii]][ppos,]
    sub_prot.info[[ii]] = prot.info[[ii]][ppos,]
    sub_unique_mm_list[[ii]] = mm_list[[ii]][!ppos,]
    sub_unique_prot.info[[ii]] = prot.info[[ii]][!ppos,]
    # each list with proteins in common between
    # datasets needs to be in the same order
    # for future comparisons.
    # Sort in increasing order of ID being used by the values in ret$ix
    # ret = sort(sub_prot.info[[ii]][,c(prot_col_name)], index.return=TRUE)
    indx = mixedorder(as.character(sub_prot.info[[ii]][,c(prot_col_name)]))
    sub_mm_list[[ii]] = sub_mm_list[[ii]][indx,]
    sub_prot.info[[ii]] = sub_prot.info[[ii]][indx,]
  }
  ret = NULL
  ret$sub_mm_list = sub_mm_list
  ret$sub_prot.info = sub_prot.info
  ret$sub_unique_mm_list = sub_unique_mm_list
  ret$sub_unique_prot.info = sub_unique_prot.info
  ret$common_list = common_list
  return(ret)
}


#' Get Presence/Absence Proteins
#'
#' Function get_presAbs_prots() produces a
#' subset of protein meta data and intencities
#' for multiple datasets pass in as a list.
#' If a single dataset is passed in
#' (list of length one) it will be processed in the same way as longer lists.
#'
#' @param mm_list list of matrices of intensities for each experiment.
#'                Dimentions: numpeptides x numsamples
#'                different for each dataset.
#' @param prot.info list of protein and peptide metadata/mappings
#'                for each matrix in mm_list, data.frames "parallel"
#'                to matrices in mm_list.
#' @param protnames_norm list of protein pdentifies to be used
#'                to determine peptides that will be placed into
#'                Presence/Absence analysis category due to
#'                too many missing peptides. Taken from the
#'                return value from eig_norm2().
#' @param prot_col_name column name (string) that will be used to get
#'                ProteinIDs in the raw data matrices
#'
#' @return list of lists of length 2
#' \describe{
#'   \item{intensities}{list of intecities in the same
#'         order and of the same length as
#'         the number of datasets that were passed into the function}
#'   \item{protein metadata}{list of protein metadata in the
#'         same order and of the same length as
#'         the number of datasets that as were passed into the function}
#' }
#'@examples
#'#' # Load mouse dataset
#' data(mm_peptides)
#' head(mm_peptides)
#' intsCols = 8:13
#' metaCols = 1:7 # reusing this variable
#' m_logInts = make_intencities(mm_peptides, intsCols)  # will reuse the name
#' m_prot.info = make_meta(mm_peptides, metaCols)
#' m_logInts = convert_log2(m_logInts)
#' grps = as.factor(c('CG','CG','CG', 'mCG','mCG','mCG'))
#' mm_m_ints_eig1 = eig_norm1(m=m_logInts,treatment=grps,prot.info=m_prot.info)
#' mm_m_ints_eig1$h.c # check the number of bias trends detected
#' mm_m_ints_norm = eig_norm2(rv=mm_m_ints_eig1)
#'
#' # Load human dataset
#' data(hs_peptides)
#' head(hs_peptides)
#' intsCols = 8:13
#' metaCols = 1:7 # reusing this variable
#' m_logInts = make_intencities(hs_peptides, intsCols)  # will reuse the name
#' m_prot.info = make_meta(hs_peptides, metaCols)
#' m_logInts = convert_log2(m_logInts)
#' grps = as.factor(c('CG','CG','CG', 'mCG','mCG','mCG'))
#' hs_m_ints_eig1 = eig_norm1(m=m_logInts,treatment=grps,prot.info=m_prot.info)
#' hs_m_ints_eig1$h.c # check the number of bias trends detected
#' hs_m_ints_norm = eig_norm2(rv=hs_m_ints_eig1)
#'
#' # Set up for presence/absence analysis
#' raw_list = list()
#' norm_imp_prot.info_list = list()
#' raw_list[[1]] = mm_m_ints_eig1$m
#' raw_list[[2]] = hs_m_ints_eig1$m
#' norm_imp_prot.info_list[[1]] = mm_m_ints_eig1$prot.info
#' norm_imp_prot.info_list[[2]] = hs_m_ints_eig1$prot.info
#'
#' protnames_norm_list = list()
#' protnames_norm_list[[1]] = unique(mm_m_ints_norm$normalized$MatchedID)
#' protnames_norm_list[[2]] = unique(hs_m_ints_norm$normalized$MatchedID)
#'
#' presAbs_dd = get_presAbs_prots(mm_list=raw_list,
#'                               prot.info=norm_imp_prot.info_list,
#'                               protnames_norm=protnames_norm_list,
#'                               prot_col_name=2)
#' @export
get_presAbs_prots = function(mm_list, prot.info,
                             protnames_norm, prot_col_name) {
  # function get_presAbs_prots() produces a subset
  # of protein meta data and intencities
  # for multiple datasets pass in as a list, single dataset
  # will be processed in the same way
  # INPUT
  #   mm_list   - list of matrices of intensities for each
  #         experiment, dimentions: numpeptides x numsamples
  #   prot.info - list of protein and peptide mappings for each
  #          matrix in mm_list, data.frame parallel to mm-list
  #   protnames_norm - list of Protein Identifies to be used to
  #         determine peptides that will be
  #         placed into Presence/Absence analysis category due
  #         to too many missing peptides
  #   prot_col_name - column name (string) that will be
  #         used to get ProteinIDs in the raw data matrices
  # OUTPUT
  #   list of lists - position 1: list of intecities,
  #            position 2: list of protein metadata,
  #            in the order matrices were pass in
  ll = length(mm_list)
  presAbs_ints = list()
  presAbs_prot.info = list()
  for(ii in seq(1,ll)) {
    # negation of these are what we want...
    prots_removed_pos = prot.info[[ii]][,c(prot_col_name)] %in%
      protnames_norm[[ii]]
    # peptides kept
    print(paste('Number of peptides normalized:',sum(prots_removed_pos) ) )
    # peptides eliminated
    print(paste('Number of peptides Pres/Abs:',sum(!prots_removed_pos) ) )

    presAbs_prot.info[[ii]] = prot.info[[ii]][!prots_removed_pos,]
    presAbs_ints[[ii]] = mm_list[[ii]][!prots_removed_pos,]
  }
  return(list(presAbs_ints, presAbs_prot.info))
}
