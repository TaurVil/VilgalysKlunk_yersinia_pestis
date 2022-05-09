import numpy as np
import sys

# can be applied per-row
def get_quantile(target, pool):
	return np.mean(pool >= target)

file1, file2 = sys.argv[1:]

candidates = np.loadtxt(file1, comments="#", delimiter="\t", skiprows=5, usecols=(0,1,6,7))
pool = np.loadtxt(file2, comments="#", delimiter="\t", skiprows=5, usecols=(0,1,6,7))

numSites = candidates.shape[0]
perc = np.apply_along_axis(get_quantile, axis=1, arr=candidates[:,3].reshape((numSites,1)), pool=pool[:,3])

report = np.hstack([candidates,perc.reshape((numSites,1))])
pool_distn = '\t'.join(['%.3g' % x for x in np.quantile(pool[:,3], [0, 0.5, 0.75, 0.975, 0.999,1])])

np.savetxt(sys.stdout.buffer, report, fmt="%.9g", delimiter="\t", 
	header="chr\tposition\ts2_hat\tMLR\tperc", comments=f"## neutral MLR quantiles: (0, 0.5, 0.75, 0.975, 0.999,1)\n##\t = {pool_distn}\n#")