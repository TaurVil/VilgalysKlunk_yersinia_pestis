'''
This script extract characteristic summary stats from empirical data and simulate reps accordingly

usage: python Simulate_CIreps_from_SNPfreq_distn.py <numRepliates> [--seed <sd> --Ne <N1,N2> --write_txt --trueS2 <list of s2s to simulate> ]

'''
import sys, os, time, re, pickle
import numpy as np
import pandas as pd

import diplo_locus.simulate as simulate
import diplo_locus.likelihood as likelihood

maindir = os.getcwd()

from GTprobs_to_counts_wFilter import readGenolik_asMatrice


def filter_missing_and_extract_SNP_GTs(pop, panels, chrs, positions):
	'''# get WG GT matrix, filter (Report and remove individuals with >50% GT calls genome-wide)
	output dict with filtered GT matrix
	'''
	GTs_by_chr = {};
	ids_to_keep = {}

	for panel in set(panels):
		filename = f'genolik.{panel}_{pop}_pre.genolik'
		filepath = os.path.join(os.getcwd(), "genoliks", filename)
		GTs_by_chr[panel] = readGenolik_asMatrice(filepath)

		sampNum = GTs_by_chr[panel][0][1].shape[1]

		total_loci = 0
		missing_loci = np.empty((0, sampNum))
		for c in range(22):
			# chrom = c+1
			N_Matrix = GTs_by_chr[panel][c][1]
			# sanity check
			assert N_Matrix.shape[1] == sampNum
			total_loci += N_Matrix.shape[0]  # nrow
			missing_loci = np.vstack((missing_loci, np.sum((N_Matrix == 0), axis=0)))  # sum up all rows by each col
		print(f'Total number of locus on panel {panel}: {total_loci}')
		# sum up all chromosome
		missing = np.sum(missing_loci, axis=0)
		# get the ids
		filtered_num = np.sum(missing >= 0.5 * total_loci)
		# filtered_ids = np.array(np.where(missing >= 0.5*total_loci)).reshape((filtered_num))
		# generate filtered Freqs
		ids_to_keep[panel] = np.array(np.where(missing < 0.5 * total_loci)).reshape((sampNum - filtered_num))
	# print(np.shape(ids_to_keep))

	Filtered_samples = {}
	for i in range(len(positions)):
		chrom, pos, panel = chrs[i], positions[i], panels[i]
		X, N, pos_list = GTs_by_chr[panel][chrom - 1]
		# extract row
		try:
			row_idx = pos_list.index(pos)
		except:
			print(pos, pos_list[:5], (pos in pos_list))
			sys.exit()
		snpX = X[row_idx, ids_to_keep[panel]]
		snpN = N[row_idx, ids_to_keep[panel]]
		Filtered_samples[f'snp{i + 1}'] = (snpX, snpN)

	return (Filtered_samples)


def bootstrap_freq_distn(sample, repNum=1e6, histResolution=1001, seed=None):
	'''Bootstrap sampling & return dict of arrays of sample frequency & their distn'''
	# SNP_distns = {}
	if seed is not None:
		np.random.seed(seed)
	# for snp in samples:
	x, n = sample
	numIndv = len(x)
	assert len(n) == numIndv

	# generate all the idx:
	samp_idxs = np.random.choice(range(numIndv), (int(repNum), numIndv))

	# get x's & n's, each row is a rep
	boot_x = x[samp_idxs]
	boot_n = n[samp_idxs]
	print(boot_x.shape, boot_n.shape)
	freqs = np.sum(boot_x, axis=1) / np.sum(boot_n, axis=1)

	# get distn
	hist = np.histogram(freqs, bins=np.linspace(0, 1, histResolution))
	xCoord = hist[1][:-1] + np.diff(hist[1]) / 2
	# distn from histogram; a tuple of np.array
	# make sure they sum to one
	distn = (xCoord, hist[0] / np.sum(hist[0]))

	return (distn)

#from Simulate_reps_from_empr_stats import simulate_reps_emprInit
import scipy.special


# parameter for the inverse frequency change
# need to make sure it works with multiple reps
# logitLambda = 1
def inverse_delta_freq(startFreq, endFreq, invSampSize, logitLambda=1):
	preLogit = scipy.special.logit(startFreq)
	postLogit = scipy.special.logit(endFreq)
	# sanity checks
	if type(preLogit) is np.array : len(preLogit) == len(postLogit)

	deltaLogit = postLogit - preLogit
	invLogit = preLogit - logitLambda * deltaLogit
	postInvFreq = scipy.special.expit(invLogit)
	# print (postInvFreq)

	if type(preLogit) is np.ndarray :
		assert type(postInvFreq) is np.ndarray
		assert len(postInvFreq) == len(startFreq)
		## make sure it makes at least some sense
		assert all( (0 <= postInvFreq) & (postInvFreq <= 1)) , postInvFreq
		## change in opposing directions
		assert all(0 >= (endFreq - startFreq) * (postInvFreq - startFreq)), (endFreq - startFreq) * (postInvFreq - startFreq)
		invSample = np.random.binomial(invSampSize, postInvFreq.reshape((len(postInvFreq),1)) )
	else:
		#assert len(postInvFreq) == 1 #TypeError: object of type 'numpy.float64' has no len()
		assert (0 >= (endFreq - startFreq) * (postInvFreq - startFreq)), (endFreq - startFreq) * (postInvFreq - startFreq)
		invSample = np.random.binomial(invSampSize, postInvFreq, size = 1 )


	return (invSample)


# init_Distns = {'lon': (bins, distn), 'den': (bins, distn)}
def simulate_reps_emprInit(repNum, init_Distns, scenarioParameters, popGenParameters, s2s, h=0.5, deltaT=1,
						   initList=("SFS", 'snp1', 'snp2', 'snp3', 'snp4')):
	# extract values
	mAlpha, mBeta = popGenParameters
	# initialize
	samples = {}
	trajs = {}
	thrownOut = {}

	# loop through all s2 values
	for s2 in s2s:
		# make containers for thess keys
		for init in initList:
			thrownOut[(s2, init)] = {}
			samples[(s2, init)] = {}
			trajs[(s2, init)] = {}

		# for group, parameters in scenarioParameters:
		for group in ('lon13', 'den13', 'lon12'):
			# unpack parameters
			sampTimes, sampSizes, selInterval, Ne = scenarioParameters[group]
			# retrieve init
			pop = group[:3]
			# get init dis
			for init in initList:
				initcond = f'{pop}_{init}'
				xcoord, distn = init_Distns[initcond]

				# initialize LL container
				## dimension: repNum x K (#sampling time points)
				samples[(s2, init)][group] = np.empty((repNum, len(sampSizes)))

				if group[3:] == "13":
					# print((selInterval[-1] - selInterval[0])/deltaT)
					trajs[(s2, init)][group] = np.empty((repNum, int((sampTimes[-1] - sampTimes[0]) / deltaT) + 2))

				# simulate reps
				r = 0
				while r < repNum:
					# repSamples = {}
					# repTrajs = {}
					oneZero = False
					oneMax = False
					# simulate
					if group[3:] == "13":
						thisSample, thisTraj = simulate.simulateSamples(Ne, [0, h * float(s2), 0], [0, float(s2), 0],
																		mAlpha, mBeta, sampTimes, sampSizes,
																		initCond='choice', numReplicates=1,
																		initProbs=distn, initValues=xcoord,
																		selectionChangeTimes=selInterval, deltaT=deltaT)
						# check this here
						oneZero |= any(thisSample[0, :] == np.zeros((thisSample.shape[1])))
						oneMax |= any(thisSample[0, :] == sampSizes)
						if (oneZero or oneMax):
							if group not in thrownOut[(s2, init)]: thrownOut[(s2, init)][group] = 0
							thrownOut[(s2, init)][group] += 1
							continue
						else:
							# store it
							samples[(s2, init)][group][r, :] = thisSample
							if group[3:] == "13": trajs[(s2, init)][group][r, :] = thisTraj
							# more replicates
							r += 1

					# when it's lon12's turn, lon13 would've been generated already
					elif group == "lon12":

						# get index
						selStartIdx = int((selInterval[0] - sampTimes[0]) / deltaT)
						selEndIdx = int((selInterval[1] - sampTimes[0]) / deltaT)
						# get delta freq
						thatTraj = trajs[(s2, init)]['lon13'][r, :]
						preSelFreq = thatTraj[selStartIdx]
						postSelFreq = thatTraj[selEndIdx]

						duringSampSize = sampSizes[1]
						assert sampSizes[1] > 10

						duringSample = inverse_delta_freq(preSelFreq, postSelFreq, duringSampSize)

						oneZero |= (duringSample <= 0)
						oneMax |= (duringSample >= duringSampSize)

						# and check if good
						# if (allZero or allMax):
						if (oneZero or oneMax):
							if group not in thrownOut[(s2, init)]: thrownOut[(s2, init)][group] = 0
							thrownOut[(s2, init)][group] += 1
							continue
						else:
							thisSample = [int(samples[(s2, init)]['lon13'][r, 0]), duringSample[0]]
							thisSample = np.array(thisSample)
							# thisSample = np.array([int(samples[(s2, init)]['lon13'][0,0]), duringSample[0]])
							samples[(s2, init)][group][r, :] = thisSample
							# don't need to save traj for this one
							# more replicates
							r += 1
					# end of while loop for rep
		# end of for loop for group/path
	# end of for loop for s2
	return (samples, trajs, thrownOut)


def main():
	NUM = sys.argv[1]
	# mutation rates for diffusion simulation
	#simAlpha = 1.25e-8
	simAlpha = 0
	simBeta = simAlpha
	# mutation rates for likelihood computation
	llAlpha = 0
	llBeta = llAlpha

	# simulation paramters
	#numReplicates = 801
	numReplicates = int(float(NUM))
	H = 0.5

	# we need to make this somewhat explicit
	# this only applies to simulation
	deltaT = 0.1

	if "--Ne" in sys.argv:
		Nes = sys.argv[(sys.argv.index("--Ne")) + 1].split(",")
		print(Nes)
		strNe = f"_Ne{'-'.join(Nes)}"
		Ne_lon, Ne_den = map(lambda x: int(float(x)), Nes)
	else:
		Ne_lon, Ne_den = 5000, 3000
		strNe = ""

	# set up paramters for analyses
	# {'ln1': 71.39858903038059,
	# 'ln2': 74.29769406253362,
	# 'ln3': 111.13117536633305,
	# 'dn1': 51.08092692249311,
	# 'dn3': 63.27491903374267}
	scenarioParameters = { # group: (samplingTimes, sampleSizes, selectionInterval, Ne)
		"lon13": ([0, 0, 15], [72, 112], [5, 8], Ne_lon),#5000
		"den13": ([0, 0, 17], [52, 64], [6, 9], Ne_den),#3000
		"lon12": ([0, 0, 8], [72, 75], [5, 8], Ne_lon),#5000
	}

	# selection coefficients
	if "--trueS2" in sys.argv:
		trueS2_str = sys.argv[(sys.argv.index("--trueS2")+1)]
		print("Check out simulations under s = ",trueS2_str)
		trueS2_list_x = list(map(float, trueS2_str.split(",")))
		#print((0 not in trueS2_list), trueS2_list)
		if 0 not in trueS2_list_x:
			trueS2_list_x = [0] + trueS2_list_x
			print(trueS2_list_x)
	else:
		#trueS2_list = [0, 0.05, 0.25, 0.5, 0.75]
		trueS2_list_x = [-0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
		trueS2_str = "-".join(str(s) for s in trueS2_list_x)
	trueS2_list = [f"{x:.4f}" for x in trueS2_list_x]
	trueS2_list_2dg = [f"{x:g}" for x in trueS2_list_x]
	print(trueS2_list)

	# do we write out txts?
	if '--write_txt' in sys.argv:
		write_txt = True
	else:
		write_txt = False

	# get seed
	if "--seed" in sys.argv:
		seed = int(sys.argv[(sys.argv.index("--seed")) + 1])
	else:
		seed = np.random.randint(1e4)

	chrs = [5, 5, 18, 2]
	positions = [96244549, 114915460, 77287776, 204738938]
	panels = ['gwas','exons','exons','gwas']
	sMLEs = [0.38887195407942104, 0.27907186611219886, 0.442118998767634, -0.258994285919875]
	rsIDs = ['rs2549794', 'rs17473484', 'rs1052025', 'rs11571319']

	import pickle
	# read init pikl if provided
	if '--init' in sys.argv:
		init_file = sys.argv[(sys.argv.index("--init")+1)]
		print(f'Loading initial distributions from {init_file}')
		with open(init_file, "rb") as ini:
			init_Distns = pickle.load(ini)
		print(init_Distns.keys())
	else:
		# bootstrap and get initial frequency
		bootstrap_file = f'bootstrap{seed}_init_Distns_for_target_SNPs.pkl'
		if os.path.isfile(bootstrap_file):  # 'init_Distns_for_target_SNPs.pkl'
			with open(bootstrap_file, "rb") as ini:
				init_Distns = pickle.load(ini)
		else:
			# get initial distn
			print(f"Bootstrapping initial frequency distribution for target SNPs with seed {seed}...")
			init_Distns = {}; bootstrap_seed = seed
			for pop in ("london", "denmark"):
				# read input and retrieve sites
				popSamples = filter_missing_and_extract_SNP_GTs(pop, panels, chrs, positions)
				for i in range(4):
					key = f'{pop[:3]}_snp{i + 1}'
					init_Distns[key] = bootstrap_freq_distn(popSamples[f'snp{i + 1}'], seed=int(bootstrap_seed))
					# update the seed after each use (same # of digits)
					bootstrap_seed = np.random.randint(10 ** (len(str(bootstrap_seed)) - 1), 10 ** (len(str(bootstrap_seed))) - 1)
					print("seed update to", bootstrap_seed)
			# save stuff
			with open(bootstrap_file, "wb") as ini:
				pickle.dump(init_Distns, ini)


	# create simulations folder if not exist
	if os.path.isdir("simulations") is False:
		os.mkdir("simulations")
	# load if already simulated
	sim_path = os.path.join("simulations", f"sim_u{simAlpha:g}_{numReplicates}reps{strNe}_seed{seed}_SNPfDistn_s{trueS2_str}.pkl")
	print(sim_path)
	if os.path.isfile(sim_path):
		print('Loading simulations with seed', seed)
		with open(sim_path, "rb") as f:
			samples, trajs, thrownOut = pickle.load(f)
	else:
		print('Start simulating with seed', seed)
		np.random.seed(seed)
		samples, trajs, thrownOut = simulate_reps_emprInit(numReplicates, init_Distns, scenarioParameters, (simAlpha, simBeta), trueS2_list, h = H, deltaT = deltaT, initList = ('snp1','snp2','snp3','snp4'))
		print(thrownOut)
		# save stuff
		with  open(sim_path,"wb") as f:
			pickle.dump( (samples, trajs, thrownOut), f)

	# now compute likelihood
	init_list = ['snp1', 'snp2', 'snp3', 'snp4'] #"SFS",
	## load pkl if already exist
	LL_pklpath = f"LL_LR_MLR_uSim{simAlpha:g}_uLL{llAlpha:g}{strNe}_seed{seed}_{numReplicates}reps_SNPfDistn_s{trueS2_str}.pkl"

	if os.path.isfile(LL_pklpath):
		print("Loading LLs from", LL_pklpath)
		with open(LL_pklpath, "rb") as LL:
			(allLL_matrice, allLR_matrice, all_MLR) = pickle.load(LL)
	else:
		from Compute_LLs_est_shat import get_lon_den_LLs, write_output

		## prep sampleSizes
		sampleSizes = {
			"lon13": np.repeat([scenarioParameters["lon13"][1]], numReplicates, axis = 0),
			"lon12": np.repeat([scenarioParameters["lon12"][1]], numReplicates, axis = 0),
			"den13": np.repeat([scenarioParameters["den13"][1]], numReplicates, axis = 0)
		}


		# parameters
		compParameters = {
			# group: (samplingTimes, selectionInterval, Ne)
			"lon13": ([0, 0, 15], [5, 8], Ne_lon),
			"den13": ([0, 0, 17], [6, 9], Ne_den),
			"lon12": ([0, 0, 8], [5, 8], Ne_lon),
		}

		## prep grid
		import diplo_locus.utility as utility

		s2_grid = utility._get_geom_grid(-0.9, 0.9, 50)
		# exclude grid points that are too small
		s2_grid = s2_grid[np.where( (np.abs(s2_grid) >= 8.9e-4) | np.isclose(s2_grid, 0))]
		## for later convenience:
		s_pairs = [(s1, s2) for (s1, s2) in zip( H*s2_grid , s2_grid) if ( abs(s1) < 1) ]
		sorted_s_pairs = sorted(list(set(s_pairs)))

		import itertools

		allLL_matrice = {} ; allLR_matrice = {}
		all_MLR = {}
		for true_s2, init in itertools.product(trueS2_list, init_list):

			LL_matrice, LR_matrice, MLR_matrice = get_lon_den_LLs(samples[(true_s2, init)], sampleSizes, compParameters, (llAlpha, llBeta), s2_grid, H)

			allLL_matrice[(true_s2, init)] = LL_matrice

			allLR_matrice[(true_s2, init)] = LR_matrice

			all_MLR[(true_s2, init)] = MLR_matrice

			# write output?
			if write_txt:
				outpre = f'sim{seed}_init{init}_s{true_s2}_{numReplicates}reps'
				comments = f'## Simulation {seed}: true s2 = {true_s2}, initial condition {init}, deltaT = {deltaT}.\n'
				comments += f'## Ne_lon = {compParameters["lon13"][2]}; Ne_den = {compParameters["den13"][2]}; u01 = {llAlpha}; u10 = {llBeta}; default numStates & deltaT\n## london pre/post samp.times = {compParameters["lon13"][0]}, selection interval = {compParameters["lon13"][1]}\n## london pre/during samp.times = {compParameters["lon12"][0]}, selection interval = {compParameters["lon12"][1]}\n## denmark pre/post samp.times = {compParameters["den13"][0]}, selection interval = {compParameters["den13"][1]}\n'
				write_output(LL_matrice, LR_matrice, MLR_matrice, np.array(range(numReplicates)), sorted_s_pairs, outpre, comment = comments)
		# save stuff
		with open(LL_pklpath,"wb") as LLoutfile:
			pickle.dump( (allLL_matrice, allLR_matrice, all_MLR) , LLoutfile)

	# and then do stat stuff
	import matplotlib.pyplot as plt
	from Compute_CI_from_sims import interpolate_and_output_CIs
	# plot violins
	figname = f'sim{seed}_{strNe}_violins_linearIntpCI.png'
	# start plotting. Clear the panel
	plt.clf()

	xCoord = np.arange(len(trueS2_list))
	colors = plt.cm.coolwarm(np.linspace(0,1, len(trueS2_list)))
	#print(colors)
	theQuantiles = [[0.025,0.975]]*(len(trueS2_list))

	fig, axs = plt.subplots(4, 1, constrained_layout=True, figsize=(6, 6)) #

	## one panel for each init
	for i, init in enumerate(init_list):
		figtitle = f'Started with {init} bootstrapped init. freq. distn.'

		# initialize containers
		shats = []
		upperCIs = []
		lowerCIs = []
		medians = []
		means = []
		for true_s2 in trueS2_list:
			shat_pool = all_MLR[(true_s2, init)]["joint_all"][:,0]
			shats.append(shat_pool)
			upperCI, median, lowerCI = np.quantile(shat_pool, [0.975, 0.5, 0.025])
			mn = np.mean(shat_pool)
			# print quantile
			print(f'{init}, s2 = {true_s2}, [2.5% 50% 97.5%] = {(lowerCI, median, upperCI)}, mean = {mn}\n')
			upperCIs.append(upperCI)
			lowerCIs.append(lowerCI)
			medians.append(median)
			means.append(mn)

		# interpolate & write out table
		CIfilename = f'sim{seed}_{init}_freqDistn{strNe}_linear-interpolated_95CIs.txt'
		# header = "true_s2\tshat_means\tshat_median\tCI_upper\tCI_lower"
		CI_matrix = interpolate_and_output_CIs(means, upperCIs, medians, lowerCIs, trueS2_list_x, CIfilename)


		axs[i].violinplot(dataset = shats, positions = xCoord, showmedians=True, showextrema=False, quantiles=theQuantiles)
		#violinDict = plt.violinplot(shats, trueS2_list, showmedians=True, showextrema=False, quantiles=theQuantiles)
		#print (violinDict['cquantiles'])

		axs[i].set_xticks (xCoord, trueS2_list_2dg)
		#axs[i].set_yticks (fontsize = 8.5), fontsize = 8.5
		axs[i].tick_params(labelsize = 7.5)
		axs[i].set_xlabel ('True Selection Coefficient', fontsize=8.5)
		#axs[i].set_ylabel ('Max-Likelihood Estimates')
		axs[i].set_ylabel ('MLE', fontsize = 8.5)

		# some jittering for points
		for (idx, thisData) in enumerate (shats):
			jitterX = np.random.normal (loc=xCoord[idx], scale=0.02, size=len(thisData))
			jitterY = np.random.normal (loc=thisData, scale=0.01)
			axs[i].plot (jitterX, jitterY, color = colors[idx], marker = "o", linestyle = "", markersize=2, alpha = 0.3)

		# plot interpolated CI borders
		# true_s2\tshat_mean\tshat_median\tCI_upper\tCI_lower
		scaled_newx = np.linspace(min(xCoord), max(xCoord), CI_matrix.shape[0])
		axs[i].plot( scaled_newx, CI_matrix[:,0], color = "black", linestyle = ":", linewidth = 1, label = "true")#
		axs[i].plot( scaled_newx, CI_matrix[:,1], color = "green", linestyle = "--", linewidth = 0.8, label = "mean")#
		axs[i].plot( scaled_newx, CI_matrix[:,2], color = "m", linestyle = "-.", linewidth = 0.8, label = "median")#
		axs[i].plot( scaled_newx, CI_matrix[:,3], color = "blue", linestyle = "--", linewidth = 0.8, label = "95% CI-lower")#
		axs[i].plot( scaled_newx, CI_matrix[:,4], color = "blue", linestyle = "-.", linewidth = 0.8, label = "95% CI-upper")#

		axs[i].set_ylim (-1.05, 1.05)
		axs[i].set_xlim (min(xCoord) - 0.5, max(xCoord) + 0.5)
		# just some mins and maxes
		axs[i].hlines(y = 0, xmin = min(xCoord) - 0.5, xmax = max(xCoord) + 0.5, colors="red" , linestyle = "-", linewidth=0.5) #[0,1]
		axs[i].set_title(figtitle, loc = "left", fontsize = 9)
		#axs[i].title.set_text(figtitle)
		#axs[i].title.set_size(9)
		if i < 1:
			axs[i].legend(loc = "lower right", ncol = 3, facecolor = (1,1,1,0), edgecolor = (1,1,1,0), fontsize = "xx-small") #, bbox_to_anchor = (0.5,0)

	#fig.tight_layout()
	plt.savefig(figname, dpi=400)



if __name__ == '__main__':
	main()







