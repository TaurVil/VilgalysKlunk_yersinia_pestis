'''
This script bootstrap samples of the target SNPs in empirical data and simulate reps accordingly

usage: python Bootstrap_CIreps_simulation.py <numRepliates> [--seed <sd> --Ne <N1,N2> --write_txt --trueS2 <list of s2s to simulate> ]

'''
import sys, os, time, re, pickle
import numpy as np
import pandas as pd

import diplo_locus.simulate as simulate
import diplo_locus.likelihood as likelihood

maindir = os.getcwd()

from GTprobs_to_counts_wFilter import readGenolik_asMatrice


def filter_missing_and_extract_SNP_GTs(pop, panels, chrs, positions, flip_pos):
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
    for i in range(4):
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
        if pos in flip_pos:
            Filtered_samples[f'snp{i + 1}'] = (snpN - snpX, snpN)
        else:    
            Filtered_samples[f'snp{i + 1}'] = (snpX, snpN)

    return (Filtered_samples)


import scipy.special


# parameter for the inverse frequency change
# need to make sure it works with multiple reps
# logitLambda = 1
def inverse_delta_freq(startFreq, endFreq, invSampSize, logitLambda=1):
    preLogit = scipy.special.logit(startFreq)
    postLogit = scipy.special.logit(endFreq)
    # sanity checks
    if type(preLogit) is np.array:
        assert len(preLogit) == len(postLogit)

    deltaLogit = postLogit - preLogit
    invLogit = preLogit - logitLambda * deltaLogit
    postInvFreq = scipy.special.expit(invLogit)
    # print (postInvFreq)

    if type(preLogit) is np.ndarray:
        assert type(postInvFreq) is np.ndarray
        assert len(postInvFreq) == len(startFreq)
        ## make sure it makes at least some sense
        assert all((0 <= postInvFreq) & (postInvFreq <= 1)), postInvFreq
        ## change in opposing directions
        assert all(0 >= (endFreq - startFreq) * (postInvFreq - startFreq)), (endFreq - startFreq) * (
                postInvFreq - startFreq)
        invSample = np.random.binomial(invSampSize, postInvFreq.reshape((len(postInvFreq), 1)))
    else:
        # assert len(postInvFreq) == 1 #TypeError: object of type 'numpy.float64' has no len()
        assert (0 >= (endFreq - startFreq) * (postInvFreq - startFreq)), (endFreq - startFreq) * (
                postInvFreq - startFreq)
        invSample = np.random.binomial(invSampSize, postInvFreq, size=1)

    return (invSample)


# init_samples = {'lon': np.array(list-of-(P_Aa + 2P_AA)]) , 'den': np.array(list-of-samples)}"SFS",
def simulate_reps_from_samples(repNum, init_Samples, scenarioParameters, popGenParameters, s2s, sampSizeLib, h=0.5,
                               deltaT=1,
                               initList=('snp1', 'snp2', 'snp3', 'snp4')):
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
            sampTimes, selInterval, Ne = scenarioParameters[group]
            # retrieve init
            pop = group[:3]
            # get init samples
            for init in initList:
                initcond = f'{pop}_{init}'
                (initXs, initNs) = init_Samples[initcond]
                sampSizes = sampSizeLib[f'{group}_{init}']
                # sanity check
                assert len(sampTimes) == len(sampSizes) + 1
                # sanity check
                assert len(initXs) == len(initNs)
                numIndv = len(initNs)

                # initialize LL container
                ## dimension: repNum x #(sampling time points)
                samples[(s2, init)][group] = np.empty((repNum, len(sampSizes)))
                # go through 13 sims first
                if group[3:] == "13":
                    # container for trajectory
                    trajs[(s2, init)][group] = np.empty((repNum, int((sampTimes[-1] - sampTimes[0]) / deltaT) + 2))
                    # simulate (repNum - r) reps per while cycle
                    r = 0
                    while r < repNum:
                        # get bootstrapped init freqs
                        resample_idxs = np.random.choice(range(numIndv), (repNum - r, numIndv))
                        resampled_Xs = initXs[resample_idxs]
                        resampled_Ns = initNs[resample_idxs]
                        init_freqs = np.sum(resampled_Xs, axis=1) / np.sum(resampled_Ns, axis=1)
                        # sanity check
                        assert np.all(init_freqs >= 0)
                        thisSample, thisTraj = simulate.simulateSamples(Ne, [0, h * float(s2), 0], [0, float(s2), 0],
                                                                        mAlpha, mBeta, sampTimes, sampSizes,
                                                                        initCond='initFreq', numReplicates=(repNum - r),
                                                                        initFreq=init_freqs,
                                                                        selectionChangeTimes=selInterval, deltaT=deltaT)
                        # check this here
                        toKeep = np.any((thisSample != np.zeros((thisSample.shape[1]))) & (thisSample != sampSizes),
                                        axis=1)
                        numToKeep = sum(toKeep)
                        if numToKeep < (repNum - r):
                            if group not in thrownOut[(s2, init)]: thrownOut[(s2, init)][group] = 0
                            thrownOut[(s2, init)][group] += (repNum - r - numToKeep)
                            thisSample = thisSample[toKeep, :]
                            thisTraj = thisTraj[toKeep, :]
                        # store it
                        samples[(s2, init)][group][r:(r + numToKeep), :] = thisSample
                        trajs[(s2, init)][group][r:(r + numToKeep), :] = thisTraj
                        # update r: number of passed reps
                        r += numToKeep
                # end of while-loop for lon13 & den13

                # when it's lon12's turn, lon13 would've been generated already. Only re-simulate if needed
                elif group == "lon12":
                    # sanity check
                    assert (s2, init) in trajs
                    assert (s2, init) in samples
                    # get traj index
                    selStartIdx = int((selInterval[0] - sampTimes[0]) / deltaT)
                    selEndIdx = int((selInterval[1] - sampTimes[0]) / deltaT)
                    # get delta freq
                    thatTraj = trajs[(s2, init)]['lon13']
                    preSelFreq = thatTraj[:, selStartIdx]
                    postSelFreq = thatTraj[:, selEndIdx]
                    duringSampSize = sampSizes[1]
                    assert sampSizes[1] > 10
                    # get samples & redo if necessary
                    reps_to_redo = set(range(repNum))
                    good_num = 0
                    while len(reps_to_redo) > 0:
                        # sanity check
                        assert good_num + len(reps_to_redo) == repNum
                        # get samples
                        duringSample = inverse_delta_freq(preSelFreq[list(reps_to_redo)],
                                                          postSelFreq[list(reps_to_redo)], duringSampSize)
                        # pair with t1 for checking
                        pairedSample = np.hstack((samples[(s2, init)]['lon13'][list(reps_to_redo), 0].reshape(
                            (len(reps_to_redo), 1)), duringSample))
                        # check if good
                        toKeep = np.any(
                            (pairedSample != np.zeros((pairedSample.shape[1]))) & (pairedSample != sampSizes), axis=1)
                        good_reps = np.array(list(reps_to_redo))[toKeep]
                        # update samples
                        samples[(s2, init)][group][good_reps, :] = pairedSample[toKeep, :]
                        # update indicator
                        reps_to_redo = reps_to_redo - set(good_reps)
                        good_num += len(good_reps)
                        # update thrownOut
                        if good_num < repNum:
                            if group not in thrownOut[(s2, init)]: thrownOut[(s2, init)][group] = 0
                            thrownOut[(s2, init)][group] += len(reps_to_redo)

                    # end of while loop. sanity check:
                    assert good_num == repNum and len(reps_to_redo) == 0
    # end of for-loop for init
    # end of for loop for group/path
    # end of for loop for s2
    return (samples, trajs, thrownOut)


from scipy import interpolate


def interpolate_and_output_CIs(means, upperCIs, medians, lowerCIs, trueS2_list_x, CIfilename):
    x = np.array(trueS2_list_x)
    f_up = interpolate.interp1d(x, upperCIs, kind="linear")  # cubic
    f_low = interpolate.interp1d(x, lowerCIs, kind="linear")  # cubic
    f_med = interpolate.interp1d(x, medians, kind="linear")  # cubic
    f_mean = interpolate.interp1d(x, means, kind="linear")  # cubic

    newx = np.linspace(min(x), max(x), (int((max(x) - min(x)) * 1000) + 1))
    # print(len(newx))

    fullmatrix = np.array([newx, f_mean(newx), f_med(newx), f_low(newx), f_up(newx)]).transpose()
    # print(fullmatrix.shape)

    np.savetxt(CIfilename, fullmatrix, fmt="%g", delimiter="\t",
               header="true_s2\tshat_mean\tshat_median\tCI_lower\tCI_upper")

    return (fullmatrix)


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("NUM", nargs="+", help="Number of replicates.")
    parser.add_argument("--Ne", dest="Nes", default=None,
                        help="Effective pop sizes of london and denmark, separated by comma.")
    parser.add_argument("-H", dest="H", type=float, default=0.5, help="Dominance Coefficient.")
    parser.add_argument("--seed", dest="seed", type=int, default=np.random.randint(1e4))
    parser.add_argument("--init", dest="initial_samples", default=None,
                        help="Path to .pkl file of initial samples.")
    parser.add_argument("--trueS2", dest="trueS2_str", help="List of alternative s2 values to simulate.")
    parser.add_argument("--write_txt", dest="write_txt", action="store_true", default=False, 
                        help="Indicator for whether to write out text outputs of likelihoods for each scenario simulated. Default is faulse.")

    if len(sys.argv) == 1:
        parser.print_usage()
        sys.exit()
        
    args = parser.parse_args(sys.argv)
    print(args)

    numReplicates = int(float(args.NUM[1]))

    # simulation paramters
    # mutation rates for diffusion simulation
    # simAlpha = 1.25e-8
    simAlpha = 0
    simBeta = simAlpha
    # mutation rates for likelihood computation
    llAlpha = 0
    llBeta = llAlpha

    # H = 0.5
    H = args.H

    # we need to make this somewhat explicit
    # this only applies to simulation
    deltaT = 0.1

    # the following four SNPs are shorthandedly denoted as snp1, 2, 3, and 4
    chrs = [5, 5, 18, 2]
    positions = [96244549, 114915460, 77287776, 204738938]
    panels = ['gwas', 'exons', 'exons', 'gwas']
    sMLEs = [0.38887195407942104, -0.27907186611219886, 0.442118998767634, 0.258994285919875]
    rsIDs = ['rs2549794', 'rs17473484', 'rs1052025', 'rs11571319']

    # need to have separate dict for sample sizes
    n_lon = [  # t1, t2, t3
        [72, 75, 112],  # WG
        [78, 68, 108],  # snp 1, 5-96244549
        [70, 78, 116],  # snp 2, 5-114915460
        [64, 64, 80],  # snp 3, 18-77287776
        [74, 76, 118]  # snp 4, 2-204738938
    ]

    n_den = [  # t1, t3
        [52, 64],  # WG
        [44, 58],  # snp 1, 5-96244549
        [60, 68],  # snp 2, 5-114915460
        [42, 62],  # snp 3, 18-77287776
        [56, 60]  # snp 4, 2-204738938
    ]

    sampleSizes_library = {}
    for i, init in enumerate(("SFS", "snp1", "snp2", "snp3", "snp4")):
        sampleSizes_library[f'lon13_{init}'] = n_lon[i][::2]
        sampleSizes_library[f'lon12_{init}'] = n_lon[i][:2]
        sampleSizes_library[f'den13_{init}'] = n_den[i]

    if args.Nes is not None:
        Nes = args.Nes.split(",")
        strNe = f"_Ne{'-'.join(Nes)}"
        Ne_lon, Ne_den = map(lambda x: int(float(x)), Nes)
    else:
        Ne_lon, Ne_den = 5000, 3000
        strNe = ""

    # set up other paramters for simulation
    scenarioParameters = {  # group: (samplingTimes, selectionInterval, Ne)
        "lon13": ([0, 0, 15], [5, 8], Ne_lon),  # 5000, [72, 112]
        "den13": ([0, 0, 17], [6, 9], Ne_den),  # 3000, [52, 64]
        "lon12": ([0, 0, 8], [5, 8], Ne_lon)  # 5000, [72, 75]
    }

    # selection coefficients
    if args.trueS2_str is not None:
        trueS2_str = args.trueS2_str
        print("Check out simulations under s = ", trueS2_str)
        trueS2_list_x = list(map(float, trueS2_str.split(",")))
        # print((0 not in trueS2_list), trueS2_list)
        if 0 not in trueS2_list_x:
            trueS2_list_x = [0] + trueS2_list_x
            print(trueS2_list_x)
    else:
        trueS2_list_x = [-0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
                         0.45, 0.5]
        trueS2_str = "-".join(str(s) for s in trueS2_list_x)
    trueS2_list = [f"{x:.4f}" for x in trueS2_list_x]
    trueS2_list_2dg = [f"{x:g}" for x in trueS2_list_x]
    print(trueS2_list)

    import pickle
    # read init pikl if provided
    if args.initial_samples is not None:
        initial_samples = sys.argv[(sys.argv.index("--init") + 1)]
        print(f'Loading initial samples from {initial_samples}')
        with open(initial_samples, "rb") as ini:
            initSamples = pickle.load(ini)
        print(initSamples.keys())
    else:
        initial_samples = 'Init_samples_for_targetSNPs.pkl'
        if os.path.isfile(initial_samples):
            with open(initial_samples, "rb") as ini:
                initSamples = pickle.load(ini)
            print(initSamples.keys())
        else:
            # get initial frequency
            initSamples = {}
            for population in ("london", "denmark"):
                # read input and retrieve sites
                popSamples = filter_missing_and_extract_SNP_GTs(population, panels, chrs, positions, flip_pos = {114915460, 204738938})
                # reformat key names
                pop = population[:3]
                for snp in ("snp1", "snp2", "snp3", "snp4"):
                    initSamples[f"{pop}_{snp}"] = popSamples[snp]
                del (popSamples)

            with open(initial_samples, "wb") as ini:
                pickle.dump(initSamples, ini)
            print(initSamples.keys())

    # create simulations folder if not exist
    if os.path.isdir("simulations") is False:
        os.mkdir("simulations")
    # load if already simulated
    sim_path = os.path.join("simulations",
                            f"sim_u{simAlpha:g}_{numReplicates}reps{strNe}_seed{args.seed}_bootstrapSNPsamp_s{trueS2_str}_h{H:g}.pkl")
    print(sim_path)
    if os.path.isfile(sim_path):
        print('Loading simulations with seed', args.seed)
        with open(sim_path, "rb") as f:
            samples, trajs, thrownOut = pickle.load(f)
    else:
        print('Start simulating with seed', args.seed)
        np.random.seed(args.seed)
        samples, trajs, thrownOut = simulate_reps_from_samples(numReplicates, initSamples, scenarioParameters,
                                                               (simAlpha, simBeta), trueS2_list,
                                                               sampSizeLib=sampleSizes_library, h=H, deltaT=deltaT,
                                                               initList=('snp1', 'snp2', 'snp3', 'snp4'))
        print(thrownOut)
        # save stuff
        with  open(sim_path, "wb") as f:
            pickle.dump((samples, trajs, thrownOut), f)

    # now compute likelihood
    init_list = ['snp1', 'snp2', 'snp3', 'snp4']  # "SFS",
    ## load pkl if already exist
    LL_pklpath = f"LL_LR_MLR_uSim{simAlpha:g}_uLL{llAlpha:g}{strNe}_seed{args.seed}_{numReplicates}reps_bootstrapSNPsamp_s{trueS2_str}_h{H:g}.pkl"

    if os.path.isfile(LL_pklpath):
        print("Loading LLs from", LL_pklpath)
        with open(LL_pklpath, "rb") as LL:
            (allLL_matrice, allLR_matrice, all_MLR) = pickle.load(LL)
    else:
        from Compute_LLs_est_shat import get_lon_den_LLs, write_output

        # parameters
        compParameters = {
            # group: (samplingTimes, selectionInterval, Ne)
            "lon13": ([0, 0, 15], [5, 8], Ne_lon),
            "den13": ([0, 0, 17], [6, 9], Ne_den),
            "lon12": ([0, 0, 8], [5, 8], Ne_lon),
        }

        ## prep grid
        from Compute_LLs_est_shat import get_geom_grid

        s2_grid = get_geom_grid(-0.9, 0.9, 50)
        # exclude grid points that are too small
        s2_grid = s2_grid[np.where((np.abs(s2_grid) >= 8.9e-4) | np.isclose(s2_grid, 0))]
        ## for later convenience:
        s_pairs = [(s1, s2) for (s1, s2) in zip(H * s2_grid, s2_grid) if (abs(s1) < 1)]
        sorted_s_pairs = sorted(list(set(s_pairs)))

        import itertools

        allLL_matrice = {}
        allLR_matrice = {}
        all_MLR = {}
        multiplier = np.ones((numReplicates, 1))
        for true_s2, init in itertools.product(trueS2_list, init_list):

            sampleSizes = {
                "lon13": multiplier * sampleSizes_library[f"lon13_{init}"],
                "lon12": multiplier * sampleSizes_library[f"lon12_{init}"],
                "den13": multiplier * sampleSizes_library[f"den13_{init}"]
            }

            LL_matrice, LR_matrice, MLR_matrice = get_lon_den_LLs(samples[(true_s2, init)], sampleSizes, compParameters,
                                                                  (llAlpha, llBeta), s2_grid, H)

            allLL_matrice[(true_s2, init)] = LL_matrice

            allLR_matrice[(true_s2, init)] = LR_matrice

            all_MLR[(true_s2, init)] = MLR_matrice

            # write output?
            if args.write_txt:
                outpre = f'sim{args.seed}_init{init}-bootstrapSNPf_s{true_s2}_{numReplicates}reps'
                comments = f'## Simulation {args.seed}: true s2 = {true_s2}, initial condition {init}, deltaT = {deltaT}.\n'
                comments += f'## Ne_lon = {compParameters["lon13"][2]}; Ne_den = {compParameters["den13"][2]}; u01 = {llAlpha}; u10 = {llBeta}; default numStates & deltaT\n## london pre/post samp.times = {compParameters["lon13"][0]}, selection interval = {compParameters["lon13"][1]}\n## london pre/during samp.times = {compParameters["lon12"][0]}, selection interval = {compParameters["lon12"][1]}\n## denmark pre/post samp.times = {compParameters["den13"][0]}, selection interval = {compParameters["den13"][1]}\n'
                write_output(LL_matrice, LR_matrice, MLR_matrice, np.array(range(numReplicates)), sorted_s_pairs,
                             outpre, comment=comments)
        # save stuff
        with open(LL_pklpath, "wb") as LLoutfile:
            pickle.dump((allLL_matrice, allLR_matrice, all_MLR), LLoutfile)

    # and then do stat stuff
    import matplotlib.pyplot as plt

    # plot violins
    figname = f'sim{args.seed}_h{H:g}-bootStrapFreq{strNe}_violins_linearIntpCI.png'
    # start plotting. Clear the panel
    plt.clf()

    xCoord = np.arange(len(trueS2_list))
    colors = plt.cm.coolwarm(np.linspace(0, 1, len(trueS2_list)))
    # print(colors)
    theQuantiles = [[0.025, 0.975]] * (len(trueS2_list))

    fig, axs = plt.subplots(4, 1, constrained_layout=True, figsize=(6, 6))  #

    ## one panel for each init
    for i, init in enumerate(init_list):
        figtitle = f'Bootstrap for {rsIDs[i]}'

        # initialize containers
        shats = []
        upperCIs = []
        lowerCIs = []
        medians = []
        means = []
        for true_s2 in trueS2_list:
            shat_pool = all_MLR[(true_s2, init)]["joint_all"][:, 0]
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
        CIfilename = f'sim{args.seed}_{init}-h{H:g}-bootStrapFreq{strNe}_linear-interpolated_95CIs.txt'
        # header = "true_s2\tshat_means\tshat_median\tCI_upper\tCI_lower"
        CI_matrix = interpolate_and_output_CIs(means, upperCIs, medians, lowerCIs, trueS2_list_x, CIfilename)

        axs[i].violinplot(dataset=shats, positions=xCoord, showmedians=True, showextrema=False, quantiles=theQuantiles)

        axs[i].set_xticks(xCoord, trueS2_list_2dg)
        # axs[i].set_yticks (fontsize = 8.5), fontsize = 8.5
        axs[i].tick_params(labelsize=7.5)
        axs[i].set_xlabel('True Selection Coefficient', fontsize=8.5)
        axs[i].set_ylabel('MLE', fontsize=8.5)

        # some jittering for points
        for (idx, thisData) in enumerate(shats):
            jitterX = np.random.normal(loc=xCoord[idx], scale=0.02, size=len(thisData))
            jitterY = np.random.normal(loc=thisData, scale=0.01)
            axs[i].plot(jitterX, jitterY, color=colors[idx], marker="o", linestyle="", markersize=2, alpha=0.3)

        # plot interpolated CI borders
        # true_s2\tshat_mean\tshat_median\tCI_upper\tCI_lower
        scaled_newx = np.linspace(min(xCoord), max(xCoord), CI_matrix.shape[0])
        axs[i].plot(scaled_newx, CI_matrix[:, 0], color="black", linestyle=":", linewidth=1, label="true")  #
        axs[i].plot(scaled_newx, CI_matrix[:, 1], color="green", linestyle="--", linewidth=0.8, label="mean")  #
        axs[i].plot(scaled_newx, CI_matrix[:, 2], color="m", linestyle="-.", linewidth=0.8, label="median")  #
        axs[i].plot(scaled_newx, CI_matrix[:, 3], color="blue", linestyle="--", linewidth=0.8, label="95% CI-lower")  #
        axs[i].plot(scaled_newx, CI_matrix[:, 4], color="blue", linestyle="-.", linewidth=0.8, label="95% CI-upper")  #

        axs[i].set_ylim(-1.05, 1.05)
        axs[i].set_xlim(min(xCoord) - 0.5, max(xCoord) + 0.5)
        # just some mins and maxes
        axs[i].hlines(y=0, xmin=min(xCoord) - 0.5, xmax=max(xCoord) + 0.5, colors="red", linestyle="-",
                      linewidth=0.5)  # [0,1]
        axs[i].set_title(figtitle, loc="left", fontsize=9)
        if i < 1:
            axs[i].legend(loc="lower right", ncol=3, facecolor=(1, 1, 1, 0), edgecolor=(1, 1, 1, 0),
                          fontsize="xx-small")  # , bbox_to_anchor = (0.5,0)

    plt.savefig(figname, dpi=400)


if __name__ == '__main__':
    main()
