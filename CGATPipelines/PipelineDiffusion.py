import CGAT.IOTools as IOTools
from CGATPipelines.Pipeline import cluster_runnable
import random
import operator
import copy


@cluster_runnable
def symbol2Ensembl(genelist, transfile, inds):
    inds = [int(x) for x in inds.split(",")]
    trans = [line.strip().split("\t")
             for line in IOTools.openFile(transfile).readlines()]
    ens = []
    for line in trans:
        if len(line) > max(inds[1], inds[0]):
            if line[inds[0]] in genelist:
                ens.append(line[inds[1]])
    return ens


@cluster_runnable
def getMatchedGenes(nmatches, network, infile, outfile, submit=True):
    ScorePair = getLinks(network)
    ScoreGene = ScoreGeneTot(ScorePair)
    genelist = [line.strip()
                for line in IOTools.openFile(infile).readlines()]
    matchedgenes = Match(ScoreGene, nmatches, genelist)
    reportMatch(outfile, matchedgenes)


@cluster_runnable
def predictionNetwork(genefile, matched, network,
                      niter, alpha, nsims, outfile):
    EvalGene = dict()
    Match = dict()
    for gene in IOTools.openFile(genefile).readlines():
        EvalGene[gene.strip()] = 1
    for line in IOTools.openFile(matched).readlines():
        line = line.strip().split("\t")
        Match.setdefault(line[0], [])
        Match[line[0]] += line[1:]
    ScorePair = getLinks(network)
    U, f0 = Create_U_f0(ScorePair, Match)
    U = NormAdjcencyMatrix(U)
    f = IterativeRanking_Lmt(U, f0, niter, alpha, EvalGene)
    fSim = {}
    f0temp = {}
    for Gene in f0.iterkeys():
        f0temp[Gene] = 0
    for i in range(nsims):
        f0sim = MakefoRandom(f0, f0temp, Match)
        fSim[i] = IterativeRanking_Lmt(U, f0sim, niter, alpha, EvalGene)
        print i
    reportPValObs(outfile, f, fSim, Match)


@cluster_runnable
def correctFDR(infile, outfile):
    score = dict()
    score_pval = dict()
    for line in IOTools.openFile(infile).readlines():
        line = line.strip().split("\t")
        score[line[0]] = float(line[1])
        score_pval[line[0]] = float(line[2])
    score_pval_fdr = compute_fdr(score_pval)
    outf = open(outfile, "w")
    for gene in score.iterkeys():
        outf.write("%s\t%s\t%s\t%s\n" %
                   (gene, score[gene], score_pval[gene], score_pval_fdr[gene]))
    outf.close()


def getLinks(infile):
    '''
    Read the network file - this consists of gene pairs and their
    linkage score as gene1\tgene2\tscore

    Returns a dictionary of dictionaries
    Each gene is the name of a subdictionary with keys corresponding
    to all other genes and values the linkage score with that gene.
    '''
    f = open(infile, 'r')
    ScorePair = {}
    for line in f:
        line = line.rstrip().split("\t")
        Gene1 = line[0]
        Gene2 = line[1]
        Score = float(line[2])
        ScorePair.setdefault(Gene1, {})
        ScorePair[Gene1][Gene2] = Score
    f.close()
    return ScorePair


def ScoreGeneTot(Score):
    '''
    Takes the total score of each gene vs every other gene in the network.
    '''
    ScoreGene = {}
    for Gene1 in Score.iterkeys():
        for Gene2 in Score[Gene1].iterkeys():
            ScoreGene.setdefault(Gene1, 0)
            ScoreGene[Gene1] += Score[Gene1][Gene2]
            ScoreGene.setdefault(Gene2, 0)
            ScoreGene[Gene2] += Score[Gene1][Gene2]
    return ScoreGene


def Match(ScoreGene, NbGene, ListGene):
    '''
    Takes the list of seed genes and extracts them from the ScoreGene
    dictionary.
    Calculates the difference between the  total linkage score of seed genes
    and all other genes in the ScoreGene dict.
    Sorts these by size of difference and takes the NbGene genes with scores
    most similar to the seed gene.
    '''
    InfoGene = {}
    for Gene1 in ScoreGene.iterkeys():
        if Gene1 not in ListGene:
            continue
        Score = {}
        InfoGene.setdefault(Gene1, [])
        ListAttribut = []
        for Gene2 in ScoreGene.iterkeys():
            if Gene1 == Gene2:
                continue
            Score = abs(ScoreGene[Gene2] - ScoreGene[Gene1])
            ListAttribut.append({'Score': Score, 'Gene': Gene2})

        SortedListAttribut = sorted(ListAttribut,
                                    key=lambda elem: elem['Score'])

        if NbGene > len(SortedListAttribut):
            NbGeneTemp = len(SortedListAttribut)
        else:
            NbGeneTemp = NbGene

        for i in range(0, NbGeneTemp):
            Gene = SortedListAttribut[i]['Gene']
            InfoGene[Gene1].append(Gene)

    return InfoGene


def reportMatch(outfile, Match):
    '''
    Outputs the matched genes for each seed gene.
    '''
    f = open(outfile, 'w')
    for Gene in Match.iterkeys():
        f.write("%s" % (Gene))
        for Gene2 in Match[Gene]:
            f.write("\t%s" % (Gene2))
        f.write("\n")
    f.close()


def Create_U_f0(ScorePair, ListGene):
    '''
    f0 = keys are all possible gene ids,
    values are 1 if the gene is in the candidate list, 0 otherwise
    U = dictionary of dictionaries.  Dictionary values are initially 1 for
    all pairs.
    '''
    U, f0 = {}, {}
    for Gene1 in ScorePair.iterkeys():
        if Gene1 in ListGene:
            f0[Gene1] = 1
        else:
            f0[Gene1] = 0
        for Gene2 in ScorePair[Gene1].iterkeys():
            if Gene2 in ListGene:
                f0[Gene2] = 1
            else:
                f0[Gene2] = 0

            U.setdefault(Gene1, {})
            U.setdefault(Gene2, {})
            U[Gene1][Gene2] = 1
            U[Gene2][Gene1] = 1

    # fill in diagonals
    for Gene in f0.iterkeys():
        U.setdefault(Gene, {})
        U[Gene][Gene] = 1

    return U, f0


def NormAdjcencyMatrix(U):
    '''
    Normalise adjacency matrix - divide each score by the total score.
    '''
    for g1 in U.iterkeys():
        total = sum([U[g1][g2] for g2 in U[g1].iterkeys()])
        for g2 in U[g1].iterkeys():
            U[g1][g2] = float(U[g1][g2])/total
    return U


def IterativeRanking_Lmt(U, f0, NbIter, alpha, GenesEval):
    '''
    Predicts the probability that each candidate gene
    is involved in the disease.
    '''
    p = copy.deepcopy(f0)

    f = {}
    geneset = set(GenesEval.keys()) & set(f0.keys())
    for i in range(NbIter):
        for g in geneset:
            if g not in U:
                f[g] = 0
            else:
                f[g] = ComputeProduct_U_ft(g, p, U)
                f[g] = f0[g]*alpha+(1-alpha)*f[g]
                p[g] = f[g]
    return f


def ComputeProduct_U_ft(g, f, U):
    product = 0
    for g2 in U[g].iterkeys():
        product += U[g][g2]*f[g2]
    return product


def MakefoRandom(f0, f0temp, Match):
    f0New = copy.deepcopy(f0temp)
    for Gene in f0.iterkeys():
        if f0[Gene] == 1:
            GeneNew = RandomGene(Gene, f0, Match)
            f0New[GeneNew] = 1

    return f0New


def RandomGene(Gene, f0, Match):
    return random.choice(Match[Gene])


def reportPValObs(outfile, f, fSim, Match):
    out = open(outfile, 'w')

    for Gene in f.iterkeys():
        nbtime, nbtotal = 0, 0
        for sim in fSim.iterkeys():
            if f[Gene] <= fSim[sim][Gene]:
                nbtime += 1
            nbtotal += 1
        pVal = float(nbtime)/float(nbtotal)
        out.write("%s\t%s\t%s\n" % (Gene, f[Gene], pVal))
    out.close()


def compute_fdr(score_pathway):
    sorted_score = sorted(score_pathway.iteritems(),
                          key=operator.itemgetter(1))
    pval = []
    for i in range(0, len(sorted_score)):
        pval.append(sorted_score[i][1])

    score_pathway_fdr = dict()
    pval_adj = fdr_adjusted(pval)
    for i in range(0, len(sorted_score)):
        term = sorted_score[i][0]
        score_pathway_fdr[term] = pval_adj[i]
    return score_pathway_fdr


def fdr_adjusted(pval):
    '''
    the p-adjusted for the current p-value is the smallest slope among
    all the slopes of each of the p-values larger than the current one
    Yekutieli & Benjamini (1999), equation
    '''

    pval.sort()
    V = len(pval)
    cV = 1
    padj = [0 for i in range(0, V)]
    prev = 1

    for i in range(V, 0, -1):
        padj[i-1] = min(prev, pval[i-1]*V*cV/i)
        prev = padj[i-1]

    return padj
