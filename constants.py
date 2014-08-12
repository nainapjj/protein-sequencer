meansAndStdDevFile = "data/loweststd.txt"
differenceInAminoFile = "data/differenceInAmino.txt"
pdbFile = "data/3ZOB.pdb"
aminoFile = "data/ListOfAminoNames.txt"

# Our current Amino Acid Count
sizeOfProtein = 67

# The top number of Amino Acids motifs to include when deciding
# which motifs to include in the trackers
NUMBER_OF_TOP_MOTIFS = 1000000

# Limit for how many times an amino acid must be included in a
# for it's amino acid to be tracked. (0 is default)
MIN_NUMBER_OF_MOTIFS_PER_AA = 0

# Determines whether we use our own gradient function or not
USE_OUR_GRADIENT = 0
USE_NO_GRADIENT = 1
USE_MINIMIZE = 2
USE_1D_CONSTRAINTS = 3

WHICH_METHOD = USE_1D_CONSTRAINTS