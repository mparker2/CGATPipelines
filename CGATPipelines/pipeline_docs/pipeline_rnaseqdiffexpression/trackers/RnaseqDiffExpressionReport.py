import re
import glob
import os

from CGATReport.Tracker import TrackerSQL

from CGATReport.Utils import PARAMS as P

from CGATPipelines.PipelineGeneset import mapUCSCToEnsembl
import CGATPipelines.PipelineTracks as PipelineTracks

# get from config file
UCSC_DATABASE = P["genome"]

###################################################################
###################################################################
# parameterization

DATABASE_ANNOTATIONS = P['annotations_database']
CUFFDIFF_LEVELS = ("gene", "isoform", "cds", "tss")


def splitLocus(locus):
    if ".." in locus:
        contig, start, end = re.match("(\S+):(\d+)\.\.(\d+)", locus).groups()
    elif "-" in locus:
        contig, start, end = re.match("(\S+):(\d+)\-(\d+)", locus).groups()

    return contig, int(start), int(end)

###########################################################################


def linkToUCSC(contig, start, end):
    '''build URL for UCSC.'''

    ucsc_database = UCSC_DATABASE
    link = "`%(contig)s:%(start)i-%(end)i <http://genome.ucsc.edu/cgi-bin/hgTracks?db=%(ucsc_database)s&position=%(contig)s:%(start)i..%(end)i>`_" \
        % locals()
    return link


def linkToEnsembl(id):
    ensembl_info = mapUCSCToEnsembl(UCSC_DATABASE)
    ensembl_species = ensembl_info.species

    if id.startswith(ensembl_info.gene_prefix):
        link = "`%(id)s <http://www.ensembl.org/%(ensembl_species)s/Gene/Summary?g=%(id)s>`_" \
            % locals()
    elif id.startswith(ensembl_info.transcript_prefix):
        link = "`%(id)s <http://www.ensembl.org/%(ensembl_species)s/Transcript/Summary?t=%(id)s>`_" \
            % locals()
    else:
        link = id
    return link


class ProjectTracker(TrackerSQL):

    '''Define convenience tracks for plots'''

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self,
                            *args,
                            attach=[(DATABASE_ANNOTATIONS, 'annotations')],
                            **kwargs)

        designs = glob.glob(os.path.join(self.datadir, "design*.tsv"))
        self.designs = sorted([os.path.splitext(os.path.basename(x))[0]
                               for x in designs])

        genesets = glob.glob(os.path.join(self.datadir, "*.gtf.gz"))
        self.genesets = sorted([os.path.splitext(os.path.basename(x))[:-7]
                                for x in genesets])

        self.samples = PipelineTracks.Tracks(PipelineTracks.Sample).loadFromDirectory(
            glob.glob("%s/*.bam" % self.datadir), "(\S+).bam")
        self.experiments = PipelineTracks.Aggregate(
            self.samples, labels=("condition", "tissue"))
        self.conditions = PipelineTracks.Aggregate(
            self.samples, labels=("condition", ))
        self.tissues = PipelineTracks.Aggregate(
            self.samples, labels=("tissue", ))
