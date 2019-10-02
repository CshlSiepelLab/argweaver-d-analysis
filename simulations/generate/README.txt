Main files in this directory:

genDeepSims.sh: simulation script to generate simulations with "deep"
 introgression (from Africans to Neanderthals, or super-archaics into
 Africans or Denisovans). Type "bash genDeepSims.sh -h" for options.

  When finished, produces simulated data sets. Each is in its own directory
    and contains files named:
    sim.sites.gz: simulated sites file (for ARGweaver input)
    aToN.txt: regions with humToNea introgression
    sToD.txt: regions with supToDen introgression
    sToA.txt: regions with supToAfr introgression
    recomb_map.bed.gz: recombination map that applies to simulated region
    mu_map.bed.gz: mutation map that applies to simulated region
    trueArg.bed.gz: local trees for simulated region

genRecentSims.sh: similar to genDeepSims.sh, but simulates Europeans and
  Africans (with realistic demographic parameters). By default it only
  has Nea->Eur introgression, but can optionally include Afr->Nea
  introgression as well. As well as the same types of files mentioned
  above, it also produces a directory "crf" which contains input files
  formatted to run the CRF software.

Other files in this directory:
  simDeep.py: helper script for genDeepSims.sh
  used_rates_to_mask.sh: helper script for genDeepSims.sh
  archaic_popsizes.txt: population sizes used by genDeepSims.sh;
    format is time (in generations), population size, population num
    (population nums: 1=Vindija; 2=Altai; 3=Denisovan)
  simFilters: show regions that are masked in real data; these files are
    used by genDeepSims.sh to replicate missing data pattern in simulated
    data
    sitesToCRF.R, par.caller.test, mcle: helper files for producing data
       formatted for the CRF software
    mu_rho.bed.gz: autosomal mutation/recombination rate maps that the
       scripts sample from
    
    
