# Configuration Files
meansAndStdDevFile = "../data/loweststd.txt"
differenceInAminoFile = "../data/differenceInAmino.txt"
hydrophobicFile = "../data/hydrophobic.txt"

# For the Indexator, choosing the best values

# For the analyzePdbVolumes analysis
# The number between which we'll filter out the actual volumes
MIN_ACTUAL_VOL = 5
MAX_ACTUAL_VOL = 17.5

# When using the Lowest Std method
# The top number of Amino Acids motifs to include when deciding
# which motifs to include in the trackers
NUMBER_OF_TOP_MOTIFS = 1000000

# When using the Lowest Vol method
# The threshold of the std. of the index values such that
# the motif amino acid are far enough part in the protein
MEAN_PAIRWISE = 3
MIN_PAIRWISE = 2
# The N_MULTIPLE * sizeOfProtein number of motifs chosen
N_MULTIPLE = 3

# Determines how the motifs for the least sqaures are chosen
USE_LOWEST_STD = 0
USE_LOWEST_VOL = 1

WHICH_CHOOSE_METHOD = USE_LOWEST_VOL

# If we're using 1D constraints, determines if we'll use a tub method
# or polynomial method.
USE_TUB = True
MAX_LENGTH = 10

# Determines which Least Square method to use
USE_1D_CONSTRAINTS = 0

WHICH_LSQ_METHOD = USE_1D_CONSTRAINTS
