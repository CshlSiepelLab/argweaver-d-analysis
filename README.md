# Scripts supporting ARGweaver-D paper

These directories contain scripts for reproducing the results and main
plots in the paper "Mapping gene flow between ancient hominins through
demography-aware inference of the ancestral recombination graph"
by M.J. Hubisz, A. L. Williams, and A. Siepel (2019)

The scripts are not published as finished software but for informational purposes
and transparency of our methods. Most of the analysis was done on a cluster, with
segments of the genome analyzed separately, and then results combined and plotted
at a later stage.

Not every detail is described within these scripts. Some details, such as where
data was downloaded, or how filtering files were generated, are described in
the methods section of the text.

The scripts assume the installation of several software/tools:

1. ARGweaver: http://github.com/CshlSiepelLab/argweaver.git
   Note that the package also includes an R package "argweaver" that is
   also used for plotting some results
2. R and packages: rphast, Hmisc, spatstat
3. bedops (http://bedops.readthedocs.io/en/latest)
4. msprime (https://msprime.readthedocs.io). We used version 0.7.0.
5. Seq-Gen (http://tree.bio.ed.ac.uk/software/seqgen/) We used version 1.1.3.


## Important files/directories:

models/: Contains all the demographic models (divergence times/migration paths)
  used in various ARGweaver-D runs. These models are passed to the
  --pop-tree-file argument of ARGweaver-D.

simulations/generate/: Contains all scripts for generating simulated data
  sets in the paper. More information given in the README.txt file therein.

simulations/analyze/: Contains all scripts for running ARGweaver on
  the simulated data sets, and plotting results.
  runArgweaverRecent.sh runs ARGweaver on a single simulated data set, with a
     single set of parameters, with some defaults used to look for recent
     introgression into non-African humans from Neanderthal
  runArgweaverDeep.sh runs ARGweaver on a single simulated data set, with
     some defaults used to look for older introgression events
  runCRF.sh runs the CRF method on a single simulated data set. Software
     was supplied by Sriram Sankararaman.
  The script ARG_vs_CRF.R produces Figure 2 of the main paper.
  The script plotTPRates.R produces Figure 5.

analysis/: Contains all scripts for running ARGweaver on data set
  consisting of some SGDP individuals, Vindija/Altai Neanderthals,
       Denisovan, and chimpanzee outgroup
  runOOA.sh runs ARGweaver on a single chromosomal segment looking for recent
     Nea or Den introgression into non-African humans
  runDeep.sh runs ARGweaver on a single segment looking for deeper (older)
     introgression events
  combineResults.sh combines all the results across many segments (using the same
      model) into a single set of predictions
  The script ooa_coverage.R produces Figure 4 of the main paper.
  The script deep_coverage.R produces Figure 6.

scripts/: Contains some miscellaneous scripts, including functions.R which
  is required by other R code in this repository


