
# MAIN FUNCTION TO CALL THE L1B MODULE

from l1b.src.l1b import l1b

# Directory - this is the common directory for the execution of the E2E, all modules
auxdir = r'C:/Users/Marcos Porto/Desktop/EODP/EODP_MP/auxiliary'
indir = r"C:/Users/Marcos Porto/Desktop/EODP/EODP-TS-ISM/myismoutput"
outdir = r"C:/Users/Marcos Porto/Desktop/EODP/EODP-TS-L1B/myl1boutput_noeq"

# Initialise the ISM
myL1b = l1b(auxdir, indir, outdir)
myL1b.processModule()
