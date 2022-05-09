'''
Read two input files of london and denmark, and compute lls.

usage: python Compute_LLs_est_shat.py <infile_suffix> <outfile_prefix>

'''
import os, sys, time, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

maindir = os.getcwd()

# local imports
import diplo_locus.likelihood as likelihood
# import diplo_locus.utility as utility

# for interpolation
from scipy import interpolate, optimize

# helper function
def get_geom_grid(left, right, grid_num):
    '''generate geometric grid points that cover every 10s'''
    # assume symmetric. If not, make it so.
    if left < 0 and right > 0: 
        bound = max(right, -1*left)
        digits = round(np.log10(bound/1e-4))
        n = (grid_num/2 -1) // digits # the # of grid pts bewteen any 10^k and 10^(k-1)
        pos_half = np.geomspace(bound/(10**digits), bound, int(n*digits + 1))
        neg_half = np.geomspace(-1*bound, -1*bound/(10**digits), int(n*digits + 1))
        s_grid = np.hstack((neg_half, 0, pos_half))
    
    elif left * right > 0: # 0 will be automatically included
        digits = round(np.log10(right/left))
        n = (grid_num -1) // digits # the # of grid pts bewteen any 10^k and 10^(k-1)
        s_grid = np.sort(np.hstack((np.geomspace(left, right, int(n*digits + 1)) , 0)))

    elif left * right == 0:
        bound = max( abs(left), abs(right))
        sign = ( (left+right) > 0 )*2 -1
        s_grid = _get_geom_grid(sign * 1e-3 * bound, sign*bound, grid_num)

    else:
        print(f'Please double-check your input for geometric grid [{left}, {right}] grid_num = {grid_num}.')
        sys.exit()

    return(s_grid)



# helpful function to be applied along axis
def findMax(ll, s2_list):
		# make an interp function
		f = interpolate.interp1d(s2_list, ll, kind="cubic")
		# the right bounds
		bounds = (min(s2_list), max(s2_list))

		# and find the optimum (minus, cause we can only minimize)
		Opt = optimize.minimize_scalar(lambda x: - f(x), method='bounded', bounds=(bounds[0], bounds[1]))
		# return the results
		maxS2 = float(Opt.x)
		maxLL = - float(Opt.fun)

		# we might be slightly negative here, so make sure to only return 0s
		if (maxLL < 0):
				# that is pretty lenient, but I guess so
				if (maxLL > -1e-4):
						maxLL = 0
				else:
						assert (False), f"Numerics for minimum too imprecise: {maxLL}."
		return np.array([maxS2, maxLL])


def get_lon_den_LLs(samples, sampleSizes, parameters, mRate, s2_grid, h):
		s2_grid = np.array(s2_grid)
		s_pairs = [(s1, s2) for (s1, s2) in zip(h * s2_grid, s2_grid) if (abs(s1) < 1)]
		## for later convenience:
		sorted_s_pairs = sorted(list(set(s_pairs)))
		num_s_pairs = len(s_pairs)
		neutIdx = sorted_s_pairs.index((0, 0))

		llAlpha, llBeta = mRate

		LL_matrice = {}
		LR_matrice = {}
		MLRs = {}
		print(time.ctime(), f"Start computing log likelihoods on {num_s_pairs} pairs of (s1, s2)...")
		for group in ("lon13", "den13", "lon12"):
				# extract parameters
				(samplingTimes, selInterval, Ne) = parameters[group]

				# make a container
				numSites = samples[group].shape[0]
				LL_matrix = np.zeros((numSites, num_s_pairs))

				# going through the grid
				for selIdx, (s1, s2) in enumerate(sorted_s_pairs):
						if group == "lon12":
								HMMgen = likelihood.SelHmm(Ne, [0, -s1, 0], [0, -s2, 0], llAlpha, llBeta, "uniform",
																				   emissionType="fractional", transitionType="piecewise",
																				   selectionChangeTimes=selInterval)
						else:
								HMMgen = likelihood.SelHmm(Ne, [0, s1, 0], [0, s2, 0], llAlpha, llBeta, "uniform",
																				   emissionType="fractional", transitionType="piecewise",
																				   selectionChangeTimes=selInterval)

						LL_matrix[:, selIdx] = HMMgen.computeLogLikelihood(samplingTimes, sampleSizes[group], samples[group])
				LL_matrice[group] = LL_matrix

		# get joint LL
		LL_matrice['joint_lon'] = LL_matrice['lon12'] + LL_matrice['lon13']
		LL_matrice['joint_all'] = LL_matrice['lon12'] + LL_matrice['lon13'] + LL_matrice['den13']

		for group in LL_matrice:
				# get likelihood ratio
				LR_matrice[group] = 2 * (
								LL_matrice[group] - np.repeat(LL_matrice[group][:, neutIdx, np.newaxis], num_s_pairs, axis=1))
				# maximize LR
				MLRs[group] = np.apply_along_axis(findMax, axis=1, arr=LR_matrice[group], s2_list=s2_grid)

		return (LL_matrice, LR_matrice, MLRs)


def write_output(LL_matrice, LR_matrice, MLR_matrice, positions, sorted_s_pairs, outpre, comment):
		# write output
		outputLLs_name = outpre + "_all_per-site_LLs.txt"
		outputMLRs_name = outpre + "_all_maxLLRs.txt"

		print(time.ctime(), "writing computed log-likelihoods to", outputLLs_name, "\nand interpolated max-LL ratios to\n",
				  outputMLRs_name)

		outputLL = open(outputLLs_name, 'w')
		outputMLR = open(outputMLRs_name, 'w')
		outputLL.write(comment)
		outputMLR.write(comment)
		# header
		if type(positions[0]) is tuple:
				outputLL.write(
						'chr\tposition\ts1\ts2\tLL_lon13\tLR_lon13\tLL_lon12\tLR_lon12\tLL_den13\tLR_den13\tLL_joint_lon\tLR_joint_lon\tLL_joint_all\tLR_joint_all\n')
				MLR_header = 'chr\tposition\ts2hat_lon\tMLR_lon\ts2hat_den\tMLR_den\ts2hat_all\tMLR_all\n'
				outputMLR.write(MLR_header)
		elif type(positions[0]) is not tuple:
				outputLL.write(
						'position\ts1\ts2\tLL_lon13\tLR_lon13\tLL_lon12\tLR_lon12\tLL_den13\tLR_den13\tLL_joint_lon\tLR_joint_lon\tLL_joint_all\tLR_joint_all\n')
				MLR_header = 'position\ts2hat_lon\tMLR_lon\ts2hat_den\tMLR_den\ts2hat_all\tMLR_all\n'
				outputMLR.write(MLR_header)

		# start writing
		numSites = len(positions)
		header_written = False
		for pos_i in range(numSites):
				target_SNP = False
				if type(positions[0]) is tuple:
						pos = "\t".join([str(p) for p in positions[pos_i]])
						if int(positions[pos_i][1]) in {96244549, 114915460, 77287776, 204738938}:
								target_SNP = True
				else:
						pos = positions[pos_i]

				outputline = f'{pos}\t{MLR_matrice["joint_lon"][pos_i, 0]}\t{MLR_matrice["joint_lon"][pos_i, 1]}\t{MLR_matrice["den13"][pos_i, 0]}\t{MLR_matrice["den13"][pos_i, 1]}\t{MLR_matrice["joint_all"][pos_i, 0]}\t{MLR_matrice["joint_all"][pos_i, 1]}\n'
				outputMLR.write(outputline)
				if target_SNP:
						if not header_written:
								print(MLR_header[:-1])
								header_written = True
						print(outputline)

				for selIdx, (s1, s2) in enumerate(sorted_s_pairs):
						outputLL.write(
								f'{pos}\t{s1}\t{s2}\t{LL_matrice["lon13"][pos_i, selIdx]}\t{LR_matrice["lon13"][pos_i, selIdx]}\t{LL_matrice["lon12"][pos_i, selIdx]}\t{LR_matrice["lon12"][pos_i, selIdx]}\t{LL_matrice["den13"][pos_i, selIdx]}\t{LR_matrice["den13"][pos_i, selIdx]}\t{LL_matrice["joint_lon"][pos_i, selIdx]}\t{LR_matrice["joint_lon"][pos_i, selIdx]}\t{LL_matrice["joint_all"][pos_i, selIdx]}\t{LR_matrice["joint_all"][pos_i, selIdx]}\n')

		# done writing, close file
		outputLL.close()
		outputMLR.close()


def main():
		insuffix, outpre = sys.argv[1], sys.argv[2]

		if "--Ne" in sys.argv:
				Nes = sys.argv[(sys.argv.index("--Ne")) + 1].split(",")
				Ne_lon, Ne_den = map(lambda x: int(float(x)), Nes)
				# change H if specified
				if len(sys.argv) > 5:
					H = float(sys.argv[3])
				else:
					H = 0.5
		else:
				Ne_lon, Ne_den = 5000, 3000
				# change H if specified
				if len(sys.argv) > 3:
					H = float(sys.argv[3])
				else:
					H = 0.5
		# infile names
		lon_infile = os.path.join(maindir, f"london{insuffix}")
		den_infile = os.path.join(maindir, f"denmark{insuffix}")

		# prepare input
		london = pd.read_csv(lon_infile, sep="\t", comment="#")
		denmark = pd.read_csv(den_infile, sep="\t", comment="#")

		## merge, only keep shared sites
		twoPops = pd.merge(london, denmark, how="inner", on=["chr", "position"], suffixes=["_lon", "_den"])
		# print(twoPops.head())
		numSites = twoPops.shape[0]
		del (london, denmark)

		## assemble sample matrices
		samples = {
				# group: sample matrix
				"lon13": np.array(twoPops.loc[:, ["x1_lon", "x3_lon"]]),
				"den13": np.array(twoPops.loc[:, ["x1_den", "x3_den"]]),
				"lon12": np.array(twoPops.loc[:, ["x1_lon", "x2_lon"]])
		}

		sampleSizes = {
				"lon13": np.array(twoPops.loc[:, ["n1_lon", "n3_lon"]]),
				"den13": np.array(twoPops.loc[:, ["n1_den", "n3_den"]]),
				"lon12": np.array(twoPops.loc[:, ["n1_lon", "n2_lon"]])
		}

		mAlpha, mBeta = 0, 0

		# parameters
		Parameters = {
				# group: (samplingTimes, selectionInterval, Ne)
				"lon13": ([0, 0, 15], [5, 8], Ne_lon),
				"den13": ([0, 0, 17], [6, 9], Ne_den),
				"lon12": ([0, 0, 8], [5, 8], Ne_lon),
		}

		# prep grid
		s2_grid = get_geom_grid(-0.9, 0.9, 50)
		# exclude grid points that are too small
		s2_grid = s2_grid[np.where((np.abs(s2_grid) >= 8.9e-4) | np.isclose(s2_grid, 0))]
		s_pairs = [(s1, s2) for (s1, s2) in zip(H * s2_grid, s2_grid) if (abs(s1) < 1)]
		## for later convenience:
		sorted_s_pairs = sorted(list(set(s_pairs)))

		# call the function
		LL_matrice, LR_matrice, MLRs = get_lon_den_LLs(samples, sampleSizes, Parameters, (mAlpha, mBeta), s2_grid, H)

		# write output
		comment = f'## Ne_lon = {Parameters["lon13"][2]}; Ne_den = {Parameters["den13"][2]}; u01 = {mAlpha}; u10 = {mBeta}; default numStates & deltaT\n## london pre/post samp.times = {Parameters["lon13"][0]}, selection interval = {Parameters["lon13"][1]}\n## london pre/during samp.times = {Parameters["lon12"][0]}, selection interval = {Parameters["lon12"][1]}\n## denmark pre/post samp.times = {Parameters["den13"][0]}, selection interval = {Parameters["den13"][1]}\n'

		write_output(LL_matrice, LR_matrice, MLRs, [(c, p) for (c, p) in zip(twoPops.chr, twoPops.position)],
								 sorted_s_pairs, outpre, comment)

		print(time.ctime(), "Done.")


if __name__ == '__main__':
		main()




