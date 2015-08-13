import itertools
import glob

from CGATReport.Tracker import TrackerSQL
from CGATReport.Utils import PARAMS as P
import CGATPipelines.Pipeline as Pipeline
import CGATPipelines.PipelineTracks as PipelineTracks

# PARAMS_PIPELINE = Pipeline.peekParameters(".",
#                                           "pipeline_chipseq.py")


# Sample = PipelineTracks.Sample3

# suffixes = ["export.txt.gz",
#             "sra",
#             "fastq.gz",
#             "fastq.1.gz",
#             "csfasta.gz"]

# TRACKS = sum(itertools.chain([PipelineTracks.Tracks(Sample).loadFromDirectory(
#     [x for x in glob.glob("%s/*.%s" % (DATADIR, s)) if "input" not in x],
#     "%s/(\S+).%s" % (DATADIR, s)) for s in suffixes]),
#     PipelineTracks.Tracks(Sample))

# Sample.setDefault("asTable")

# ALL = PipelineTracks.Aggregate(TRACKS)
# EXPERIMENTS = PipelineTracks.Aggregate(TRACKS, labels=("condition", "tissue"))
# CONDITIONS = PipelineTracks.Aggregate(TRACKS, labels=("condition", ))
# TISSUES = PipelineTracks.Aggregate(TRACKS, labels=("tissue", ))

# ############################################################################
# # The folllowing need to be parameterized in a config file
# # TISSUES=["GM00855", "GM00861" ]
# # CONDITIONS=["D3", "unstim" ]
# # REPLICATES=["R1", "R2" ]
# TAG_UNSTIM = PARAMS_PIPELINE["tracks_unstimulated"]
# UCSC_GENOME = PARAMS_PIPELINE["genome"]

# if "motifs_plot" in P and P["motifs_plot"]:
#     MOTIFS = [x.strip() for x in P["motifs_plot"].split(",")]
# else:
#     MOTIFS = None

# ###########################################################################
# # shorthand
# # use list to convert trackers to strings
# MAP_TRACKS = {
#     'master': map(str, list(EXPERIMENTS) + list(CONDITIONS)),
#     'replicates': map(str, list(TRACKS)),
#     'default': map(str, list(EXPERIMENTS)),
#     'experiments': map(str, list(EXPERIMENTS)),
#     'conditions': map(str, list(CONDITIONS)),
#     'tissues': map(str, list(TISSUES)),
#     'merged': map(str, list(EXPERIMENTS)),
# }


def getReplicates(track):
    '''return replicates for a track.'''


class CallingTracker(TrackerSQL):
    '''Define convenience tracks for plots'''


class DefaultTracker(CallingTracker):
    '''Define convenience tracks for plots'''
