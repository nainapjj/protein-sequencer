meansAndStdDevFile = "data/loweststd.txt"
differenceInAminoFile = "data/differenceInAmino.txt"

# The top number of Amino Acids motifs to include when deciding
# which motifs to include in the trackers
NUMBER_OF_TOP_MOTIFS = 1000000

# If we're using 1D constraints, determines if we'll use a tub method
# or polynomial method.
USE_TUB = True
MAX_LENGTH = 10

# Determines whether we use our own gradient function or not
USE_1D_CONSTRAINTS = 0

WHICH_METHOD = USE_1D_CONSTRAINTS
