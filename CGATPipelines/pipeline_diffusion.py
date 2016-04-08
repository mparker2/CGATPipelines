"""

===========================
Diffusion analysis pipeline
===========================

:Author: Katherine Fawcett
:Release: $Id$
:Date: |today|
:Tags: Python

The diffusion analysis pipeline imports gene lists (eg. candidate genes from exome sequencing
and ranks the genes according to how closely connected they are with seed gene (eg. known
disease genes) within a given network.  The pipeline performs the following steps:

1.  Converts the gene names in the gene lists into ensembl identifiers
2.  Matches seed genes to random genes for building a null distribution
2.  Removes duplicate entries in the gene list
3.  Selects only one ensembl identifier per gene (where there are >1 identifier)
4.  Performs diffusion analysis
5.  Corrects the p values for multiple testing
6.  Converts ensembl identifiers in the output back into gene names

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

Input
-----

Input gene lists are simply text files containing a list of gene names (one gene per line) places in
the working directory and labelled '<list_id>-candidates.txt' and a list of seed genes labelled
'<seed_id>-seeds.txt'.

Documentation
-------------

Requirements
------------

The following scripts must be available within the working directory:

match_gene_degree.py
pipeline_diffusion.py
prediction_network.py
correction_fdr.py
data.py
fdr.py
functions.py

The ini file must provide a path to:
1.  get_ensg.py
2.  the network file.

Pipeline output
===============

Three files:
1.  <list_id>__GenesIheritedVariants_prediction contains columns gene, score, pvalue
2.  <list_id>__GenesIheritedVariants_prediction_fdr contains columns gene, score, pvalue, fdr-corrected p value
3.  <list_id>__GenesIheritedVariants_prediction_fdr_ord is file 2 sorted by score

Code
====

"""

# load modules
from ruffus import *
from ruffus.combinatorics import *

import sys
import os
import CGAT.IOTools as IOTools
import CGATPipelines.Pipeline as P
import PipelineDiffusion
import CGAT.Experiment as E

USECLUSTER = True

# load options from the config file


P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini", "pipeline.ini"])
PARAMS = P.PARAMS

#########################################################################
#########################################################################
# Get  ensembl ids for candidate gene lists


@transform('*seed.tsv', suffix(".tsv"), "_clean.tsv")
def translateSeeds(infile, outfile):
    genes = [line.strip()
             for line in IOTools.openFile(infile).readlines()]
    seeds = PipelineDiffusion.symbol2Ensembl(
        genes, PARAMS["ensembl_ids"], PARAMS["ensembl_inds"])
    if PARAMS["diffusion_dedup"]:
        # get unique seeds?
        seeds = set(seeds)

    outf = IOTools.openFile(outfile, "w")
    outf.write('\n'.join(seeds))
    outf.close()

@follows(mkdir("clean_candidates.dir"))
@transform("candidates.dir/*.tsv", regex("candidates.dir/(.*)-target.tsv"),
           add_inputs(r'\1-seed_clean.tsv'),
           r"clean_candidates.dir/\1-target.tsv")
def cleanGeneLists(infiles, outfile):
    '''Get ensembl ids for candidate gene lists'''
    cand_file, seed_file = infiles
    genes = [line.strip()
             for line in IOTools.openFile(cand_file).readlines()]
    cands = PipelineDiffusion.symbol2Ensembl(
        genes, PARAMS["ensembl_ids"], PARAMS["ensembl_inds"],
        submit=False)
    seeds = [line.strip()
             for line in IOTools.openFile(seed_file).readlines()]
    cands = set(cands) - set(seeds)
    outf = IOTools.openFile(outfile, "w")
    outf.write('\n'.join(cands))
    outf.close()


@transform(translateSeeds, suffix("_clean.tsv"), "_matches.tsv")
def getMatches(infile, outfile):
    '''Find matching genes for seeds.
    These are used to generate null expectation'''

    nmatches = PARAMS["diffusion_nmatches"]
    network = PARAMS["diffusion_network"]

    PipelineDiffusion.getMatchedGenes(nmatches, network, infile, outfile,
                                      job_memory='2G', submit=True)


@follows(mkdir("diffusion.dir"))
@follows(getMatches)
@transform(cleanGeneLists, regex("clean_candidates.dir/(.*)-target.tsv"),
           add_inputs(r'\1-seed_matches.tsv'),
           r"diffusion.dir/\1-target.tsv")
def buildDiffusionNetwork(infiles, outfile):
    '''Perform diffusion analysis for each candidate gene list'''
    infile, seedfile = infiles

    network = PARAMS["diffusion_network"]
    niter = PARAMS["diffusion_niter"]
    alpha = PARAMS["diffusion_alpha"]
    nsims = PARAMS["diffusion_nsims"]

    PipelineDiffusion.predictionNetwork(infile, seedfile, network, niter,
                                        alpha, nsims, outfile,
                                        job_memory='2G', submit=True)


@transform(buildDiffusionNetwork, suffix(".tsv"), ["_fdr.tsv", "_fdr_ord.tsv"])
def correctFDR(infile, outfiles):
    PipelineDiffusion.correctFDR(infile, outfiles[0])
    statement = '''sort -grk 2 %s > %s''' % (outfiles[0], outfiles[1])
    P.run()


@transform(correctFDR, suffix(".tsv"), "_symbol.tsv")
def convertGeneNames(infiles, outfile):
    '''convert ensembl ids back to gene names in output files'''
    infile = infiles[0]
    genes1, scores1, scores2 = [], [], []
    for line in open(infile).readlines():
        line = line.strip().split("\t")
        genes1.append(line[0])
        scores1.append(line[1])
        scores2.append(line[2])

    e_inds_2 = ",".join(PARAMS['ensembl_inds'].split(",")[::-1])

    symbols1 = PipelineDiffusion.symbol2Ensembl(
        genes1, PARAMS["ensembl_ids"], e_inds_2)

    z = zip(symbols1, scores1, scores2)
    
    z.sort(key=lambda t: (t[1]), reverse=True)
    z.sort(key=lambda t: (t[2]))

    outf = open(outfile, "w")
    for row in z:
        outf.write("%s\n" % ("\t".join(row)))
    outf.close()
# ###########################################################################


@follows(convertGeneNames)
def full():
    pass

###########################################################################


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
