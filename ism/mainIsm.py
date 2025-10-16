
# MAIN FUNCTION TO CALL THE ISM MODULE

from ism.src.ism import ism

# Directory - this is the common directory for the execution of the E2E, all modules
auxdir = r"C:\\Users\\Marcos Porto\\Desktop\\EODP\\EODP_MP\\auxiliary"
indir = r"C:\\Users\\Marcos Porto\\Desktop\\EODP\\EODP-TS-E2E\sgm_out" # small scene
outdir = r"C:\\Users\\Marcos Porto\\Desktop\\EODP\\EODP-TS-ISM\\myismoutput"

# Initialise the ISM
myIsm = ism(auxdir, indir, outdir)
myIsm.processModule()
