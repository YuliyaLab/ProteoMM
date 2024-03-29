% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TwoPart_MultiMS.R
\name{prot_level_multi_part}
\alias{prot_level_multi_part}
\title{Multi-Matrix Differentia Expression Analysis}
\usage{
prot_level_multi_part(
  mm_list,
  treat,
  prot.info,
  prot_col_name,
  nperm = 500,
  dataset_suffix
)
}
\arguments{
\item{mm_list}{list of matrices for each experiment, length = number of
datasets to compare internal dataset dimentions: numpeptides
x numsamples for each dataset}

\item{treat}{list of data frames with treatment information to compute
the statistic
in same order as mm_list}

\item{prot.info}{list of protein and peptide mapping for each matrix
in mm_list,
in same order as mm_list}

\item{prot_col_name}{column name in prot.info that contains protein
identifiers that link all datasets together. Not that
Protein IDs will differ across different organizms and
cannot be used as the linking identifier.
Function match_linker_ids() produces numeric identifyers
that link all datasets together}

\item{nperm}{number of permutations, default = 500,
this will take a while, test code
with fewer permutations}

\item{dataset_suffix}{vector of character strings that corresponds to the
dataset being analysed. Same length as mm_list. Names will be appended
to the columns names that will be generated for each analysed dataset.
For example, if analysing mouse and human data this vector may be:
c('Mouse', 'Human')}
}
\value{
data frame with the following columns
\describe{
  \item{protIDused}{Column containing the protien IDs used to
                    link proteins across datasets}
  \item{FC}{Average fold change across all datasets}
  \item{P_val}{Permutation-based p-valu for the differences
               between the groups}
  \item{BH_P_val}{Multiple testing adjusted p-values}
  \item{statistic}{Statistic computed as a a sum of statistics
                   produced for each dataset}
  \item{Protein Information}{all columns passed into the function
                             for the 1st dataset
        in the list}
  \item{FCs}{Fold changes for individual datasets, these values
             should average to the
        FC above. As many columns as there are datasets being analyzed.}
  \item{PV}{p-values for individual datasets. As many
            columns as there are datasets
        being analyzed.}
  \item{BHPV}{Multiple testing adjusted p-values for
              individual datasets. As many
        columns as there are datasets being analyzed.}
  \item{NUMPEP}{Number of peptides presents in each protien
                for each dataset. As many
        columns as there are datasets being analyzed.}
}
}
\description{
Multi-Matrix Differential Expression Analysis computes Model-Based
statistics for each dataset, the sum of individual statistics is the
final statistic. The significance is determined via a permutation test
which computed the same statistics and sums them after permuting
the values across treatment groups. As is outlined in Karpievitch
et al. 2018. Important to set the random number generator seed for 
reprodusibility with set.seed() function.
}
\examples{
# Load mouse dataset
data(mm_peptides)
head(mm_peptides)
intsCols = 8:13 # different from parameter names as R uses
                # outer name spaces if variable is undefined
metaCols = 1:7 # reusing this variable
m_logInts = make_intencities(mm_peptides, intsCols)  # will reuse the name
m_prot.info = make_meta(mm_peptides, metaCols)
m_logInts = convert_log2(m_logInts)
grps = as.factor(c('CG','CG','CG', 'mCG','mCG','mCG'))
set.seed(135)
mm_m_ints_eig1 = eig_norm1(m=m_logInts,treatment=grps,
                           prot.info=m_prot.info)
mm_m_ints_eig1$h.c # check the number of bias trends detected
mm_m_ints_norm = eig_norm2(rv=mm_m_ints_eig1)
mm_prot.info = mm_m_ints_norm$normalized[,1:7]
mm_norm_m =  mm_m_ints_norm$normalized[,8:13]
set.seed(125) # Needed for reprodicibility of results
imp_mm = MBimpute(mm_norm_m, grps, prot.info=mm_prot.info,
                  pr_ppos=2, my.pi=0.05, compute_pi=FALSE)

# Load human dataset
data(hs_peptides)
head(hs_peptides)
intsCols = 8:13 # different from parameter names as R uses
                # outer name spaces if variable is undefined
metaCols = 1:7 # reusing this variable
m_logInts = make_intencities(hs_peptides, intsCols)  # will reuse the name
m_prot.info = make_meta(hs_peptides, metaCols)
m_logInts = convert_log2(m_logInts)
grps = as.factor(c('CG','CG','CG', 'mCG','mCG','mCG'))
set.seed(1237) # needed for reproducibility
hs_m_ints_eig1 = eig_norm1(m=m_logInts,treatment=grps,prot.info=m_prot.info)
hs_m_ints_eig1$h.c # check the number of bias trends detected
hs_m_ints_norm = eig_norm2(rv=hs_m_ints_eig1)
hs_prot.info = hs_m_ints_norm$normalized[,1:7]
hs_norm_m =  hs_m_ints_norm$normalized[,8:13]

set.seed(125) # or any value, ex: 12345
imp_hs = MBimpute(hs_norm_m, grps, prot.info=hs_prot.info,
                  pr_ppos=2, my.pi=0.05,
                  compute_pi=FALSE)

# Multi-Matrix Model-based differential expression analysis
# Set up needed variables
mms = list()
treats = list()
protinfos = list()
mms[[1]] = imp_mm$y_imputed
mms[[2]] = imp_hs$y_imputed
treats[[1]] = grps
treats[[2]] = grps
protinfos[[1]] = imp_mm$imp_prot.info
protinfos[[2]] = imp_hs$imp_prot.info
nperm = 50

# ATTENTION: SET RANDOM NUMBER GENERATOR SEED FOR REPRODUCIBILITY !!
set.seed(131) # needed for reproducibility

comb_MBDE = prot_level_multi_part(mm_list=mms, treat=treats,
                                  prot.info=protinfos,
                                  prot_col_name='ProtID', nperm=nperm,
                                  dataset_suffix=c('MM', 'HS'))

# Analysis for proteins only present in mouse,
# there are no proteins suitable for
# Model-Based analysis in human dataset
subset_data = subset_proteins(mm_list=mms, prot.info=protinfos, 'MatchedID')
mm_dd_only = subset_data$sub_unique_mm_list[[1]]
hs_dd_only = subset_data$sub_unique_mm_list[[2]]
protinfos_mm_dd = subset_data$sub_unique_prot.info[[1]]
DE_mCG_CG_mm_dd = peptideLevel_DE(mm_dd_only, grps,
                                  prot.info=protinfos_mm_dd, pr_ppos=2)

}
