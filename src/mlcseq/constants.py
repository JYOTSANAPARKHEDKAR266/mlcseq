"""Mass constants (Da) used by MLC-Seq.

All values are copied unchanged from the original repository so results match.
NOTE: the original code uses three slightly different methyl masses in different
places (14.01, 14.0106, 14.0156). The chemically correct monoisotopic CH2 shift is
14.01565. They are kept as-is here for reproducibility and named separately so
the inconsistency is visible.
"""

# Monoisotopic nucleotide residue masses
A = 329.0525
C = 305.0413
G = 345.0474
U = 306.0253

H2O = 18.0106
P = 79.9663          # phosphate
NA = 21.9819         # sodium adduct
K = 37.9559          # potassium adduct
CO = 27.9949

METHYL_UTILS = 14.01      # used for adduct/shift search in the total-tRNA workflow
METHYL_BASECALL = 14.0106 # used by homology search ("Me")
METHYL_GAPFILL = 14.0156  # used by GapFill composition search (bases_methyl.csv)

# MassSum
MASS_SUM_P = 79.96633
MASS_SUM_TAG = 826.3184

# GapFill ladder-end tags
TAG_5P = 79.9663
TAG_3OH = 18.0106

# Rough average residue mass used to convert a fragment mass to a position.
# (Original behaviour; can drift by >1 position for heavily modified tRNAs.)
POSITION_MASS = 320

# Intact tRNA lower mass bound used by homology search
INTACT_MIN_MASS = 23000

# Base-calling dictionaries (original values)
BASES_DEFAULT = {
    "C": 305.0413, "A": 329.0525, "G": 345.0474, "U": 306.0253,
    "mG": 359.0631, "D": 308.041, "mC": 319.0569, "mU": 320.041,
    "mA": 343.0682, "m22C": 333.0625, "mcm5U": 378.0464,
}
BASES_HOMOLOGY = {"C": 305.0413, "A": 329.0525, "Ox": 16.0, "Me": 14.0106, "2*Me": 28.0212}
BASES_ACID_LABILE = {"Y": 358.1599}
BASES_METHYL = {"Me": 14.0106}

# Building blocks for GapFill composition search, in the original column order
# (A, C, G, U, D, y, Methyl) from statics/bases_methyl.csv.
GAPFILL_BLOCKS = {
    "A": 329.0525, "C": 305.0413, "G": 345.0474, "U": 306.0253,
    "D": 308.0410, "y": 212.0106, "Methyl": 14.0156,
}
