'''
This script bootstraps init freq distns from GTprobs, reads (interpolated) CIs from file, look up true s2, and simulate four SNPs with their own initial frequency and s2.
The 4 target snps are: positions = [96244549, 114915460, 77287776, 204738938] (chrs = [5,5,18,2], panels = [gwas, exons, exons, gwas])

usage: python Simulate_target_SNPs.py <numReps> <CIfile>  --init <init_distns.pkl> [--seed <sd> --Ne <N1,N2> ]
'''
import sys, os, time, re, pickle
import numpy as np
import pandas as pd

import diplo_locus.simulate as simulate
import diplo_locus.likelihood as likelihood

maindir = os.getcwd()


def get_trueS2_from_stat(CIfile, list_of_shat, stat="median"):
    '''
    Replace "SFS" in the CIfile to snp1/2/3/4, with the list_of_shat being of length 4
    # CIfile header: 0true_s2   1shat_mean  2shat_median    CI_lower    CI_upper
    if stat == "median", stat_idx = 2; elif "mean", then stat_idx = 1
    '''
    if stat == "median":
        stat_idx = 2
    elif stat == "mean":
        stat_idx = 1
    else:
        print("Please choose from stat=\"median\" or stat=\"mean\".")
        return (False)
    true_s2s = []
    for i, shat in enumerate(list_of_shat):
        snpCIfile = re.sub("SFS", f"snp{i + 1}", CIfile)
        CIs = np.loadtxt(snpCIfile, delimiter="\t", comments="#")
        # find the closest median
        idx = np.argmin(np.abs(CIs[:, stat_idx] - shat))
        true_s2s.append(CIs[idx, 0])
    return (true_s2s)


def get_tildeS2s_from_hatS2s(CIfile, shat_list, stat="median"):
    '''
    Find the matching "true s2"s in the CIfile for each shat in the shats list
    if method == "median", shats / hat_s2 will be the interpolated median (idx=2) of the true s2 / tilde_s2;
    if "mean", then idx=1
    file header/columns: # 0true_s2 1shat_mean  2shat_median    CI_lower    CI_upper
    '''
    # read file, make sure stat is monotone
    CIs = np.loadtxt(CIfile)
    if stat == "median":
        stat_idx = 2
    elif stat == "mean":
        stat_idx = 1
    else:
        print("Please choose from stat=\"median\" or stat=\"mean\".")
        return (False)
    stilde_list = CIs[:, 0]
    stat_list = CIs[:, stat_idx]
    # check for monotone
    assert np.all(np.diff(stat_list) > 0)
    # order shats, and get the indice from/to the sorted lists
    shat_list = np.sort(shat_list)
    # stilde_idx = np.apply_along_axis(lambda x: np.argmin(abs(stat_list - x)) , 0, shat_list)
    matching_tildeS2s = []
    stat_idx = 0
    correction_above = stilde_list[-1] - stat_list[-1]
    correction_bellow = stilde_list[0] - stat_list[0]
    for i, s in enumerate(shat_list):
        # if lower than lowest reference, use correction_bellow
        if stat_list[0] > s:
            try:
                assert stat_idx == 0
            except:
                print(stat_idx, len(stat_list), stat_list[stat_idx], s, i)
                sys.exit()
            stilde = s + correction_bellow
        # if higher than highest reference, use correction_above
        elif stat_list[-1] <= s:
            try:
                assert stat_idx >= len(stat_list)-1
            except:
                print(stat_idx, len(stat_list), stat_list[stat_idx], s, i)
                sys.exit()
            stilde = s + correction_above
        # if in between, find match
        else:
            while stat_list[stat_idx] <= s:
                stat_idx += 1
                if stat_idx >= len(stat_list):
                    stat_idx -= 1
                    break
            stilde = stilde_list[stat_idx]
        # save value
        matching_tildeS2s.append(stilde)
    # in the end, all slots should be filled
    matching_tildeS2s = np.array(matching_tildeS2s)
    try:
        assert all((matching_tildeS2s > -1) & (matching_tildeS2s < 1) )
    except:
        # get the unfilled ones
        unmatched_up = np.any(matching_tildeS2s > 1)
        unmatched_down = np.any(matching_tildeS2s < -1)
        print(shat_list[unmatched_up], shat_list[unmatched_down])
        print(np.quantile(shat_list,[0,0.025,0.5,0.975,1]))
        print(np.quantile(matching_tildeS2s,[0,0.025,0.5,0.975,1]))
        sys.exit()

    return (matching_tildeS2s)


from Bootstrap_CIreps_simulation import inverse_delta_freq


def simulate_snps(repNum, SNPs, init_Samples, scenarioParameters, popGenParameters, sampSizeLib, h=0.5, deltaT=1):
    '''simulate the four snps.
    SNPs = ((s2_snp1, 'snp1'), (s2_snp2, 'snp2'),...)
    '''
    # extract values
    mAlpha, mBeta = popGenParameters
    # initialize
    samples = {}
    trajs = {}
    thrownOut = {}

    # make containers for these keys
    # snp = (s2, init)
    for (s2, init) in SNPs:
        thrownOut[(s2, init)] = {}
        samples[(s2, init)] = {}
        trajs[(s2, init)] = {}

        # for group, parameters in scenarioParameters:
        for group in ('lon13', 'den13', 'lon12'):
            # unpack parameters
            sampTimes, selInterval, Ne = scenarioParameters[group]
            # retrieve init
            pop = group[:3]

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
                    toKeep = np.any((thisSample != np.zeros((thisSample.shape[1]))) & (thisSample != sampSizes), axis=1)
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
                    duringSample = inverse_delta_freq(preSelFreq[list(reps_to_redo)], postSelFreq[list(reps_to_redo)],
                                                      duringSampSize)
                    # pair with t1 for checking
                    pairedSample = np.hstack((samples[(s2, init)]['lon13'][list(reps_to_redo), 0].reshape(
                        (len(reps_to_redo), 1)), duringSample))
                    # check if good
                    toKeep = np.any((pairedSample != np.zeros((pairedSample.shape[1]))) & (pairedSample != sampSizes),
                                    axis=1)
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
        # end of for loop for group/path
    # end of for loop for s2, init
    return (samples, trajs, thrownOut)


def get_ROCs(neut_data, sel_data):
    # make sure things are sorted, default is ascending
    neut_data = np.sort(neut_data)
    sel_data = np.sort(sel_data)
    TP, FP, TN, FN = [], [], [], []
    for thr in range(len(neut_data)):
        TP.append(np.sum(sel_data >= thr))
        FP.append(np.sum(neut_data >= thr))
        TN.append(np.sum(neut_data < thr))
        FN.append(np.sum(sel_data < thr))
    return (TP, FP, TN, FN)


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("NUM", nargs="+", help="Number of replicates.")
    parser.add_argument("CIfile", nargs="+", help="Name of the interpolated CI file (with SFS in filename).")
    parser.add_argument("--Ne", dest="Nes", default=None,
                        help="Effective pop sizes of london and denmark, separated by comma.")
    parser.add_argument("-H", dest="H", type=float, default=0.5, help="Dominance Coefficient.")
    parser.add_argument("--seed", dest="seed", type=int, default=np.random.randint(1e4))
    parser.add_argument("--init", dest="init_file", required=True, default=None,
                        help="Path to .pkl file of initial samples.")

    if len(sys.argv) == 1:
        parser.print_usage()
        sys.exit()
        
    args = parser.parse_args(sys.argv)
    print(args)

    numReplicates = int(float(args.NUM[1]))
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

    if args.Nes is not None:
        Nes = args.Nes.split(",")
        strNe = f"_Ne{'-'.join(Nes)}"
        Ne_lon, Ne_den = map(lambda x: int(float(x)), Nes)
    else:
        Ne_lon, Ne_den = 5000, 3000
        strNe = ""

    print(f'Loading initial distributions from {args.init_file}')
    with open(args.init_file, "rb") as ini:
        init_Samples = pickle.load(ini)
    print(init_Samples.keys())

    scenarioParameters = {  # group: (samplingTimes, selectionInterval, Ne)
        "lon13": ([0, 0, 15], [5, 8], Ne_lon),  # 5000, [72, 112]
        "den13": ([0, 0, 17], [6, 9], Ne_den),  # 3000, [52, 64]
        "lon12": ([0, 0, 8], [5, 8], Ne_lon),  # 5000, [72, 75]
    }

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

    # loci info
    chrs = [5, 5, 18, 2]
    positions = [96244549, 114915460, 77287776, 204738938]
    panels = ['gwas', 'exons', 'exons', 'gwas']
    if H == 0.5:
        sMLEs = [0.38887195407942104, -0.27907186611219886, 0.442118998767634, 0.258994285919875]
    elif H == 0:
        sMLEs = [0.42793722197544276, -0.8999933687735545, 0.25248616107249156, 0.7060668542697011]
    elif H == 1:
        sMLEs = [0.34194682540089594, -0.15416217336815108, 0.8999933687735545, 0.15627053193909918]
    rsIDs = ['rs2549794', 'rs17473484', 'rs1052025', 'rs11571319']

    # get s2
    trueS2s = get_trueS2_from_stat(args.CIfile[0], sMLEs, stat="median")
    trueS2s_str = "-".join([f'{s:.4f}' for s in trueS2s])
    print(trueS2s_str)
    # pair with nulls
    trueS2s = trueS2s + [0, 0, 0, 0]
    trueS2s_str_list = [f'{s:.4f}' for s in trueS2s]
    SNPs = tuple(zip(trueS2s_str_list, ('snp1', 'snp2', 'snp3', 'snp4') * 2))
    print(SNPs)

    # load if already simulated
    pickle_path = os.path.join("simulations",
                               f"sim_u{simAlpha:g}_{numReplicates}reps{strNe}_seed{args.seed}_s{trueS2s_str}_h{H:g}.pkl")
    if os.path.isfile(pickle_path):
        print('Loading simulations with seed', args.seed)
        with open(pickle_path, "rb") as f:
            samples, trajs, thrownOut = pickle.load(f)
    else:
        # now simulate
        print('Start simulating with seed', args.seed)
        np.random.seed(args.seed)
        samples, trajs, thrownOut = simulate_snps(numReplicates, SNPs, init_Samples, scenarioParameters,
                                                  (simAlpha, simBeta), sampSizeLib=sampleSizes_library, h=H,
                                                  deltaT=deltaT)
        print(time.ctime(), thrownOut)
        # save stuff
        with open(pickle_path, "wb") as f:
            pickle.dump((samples, trajs, thrownOut), f)

    # get likelihood
    LLfilename = f"LL_LR_MLR_uSim{simAlpha:g}_uLL{llAlpha:g}{strNe}_seed{args.seed}_{numReplicates}reps_s{trueS2s_str}_h{H:g}.pkl"
    ## if file exist, load it
    if os.path.isfile(LLfilename):
        with open(LLfilename, 'rb') as f:
            allLL_matrice, allLR_matrice, all_MLRs = pickle.load(f)
    else:
        from Compute_LLs_est_shat import get_lon_den_LLs

        ## scan parameters
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

        allLL_matrice = {}
        allLR_matrice = {}
        all_MLRs = {}
        multiplier = np.ones((numReplicates, 1))
        for (true_s2, init) in SNPs:
            sampleSizes = {
                "lon13": multiplier * sampleSizes_library[f"lon13_{init}"],
                "lon12": multiplier * sampleSizes_library[f"lon12_{init}"],
                "den13": multiplier * sampleSizes_library[f"den13_{init}"]
            }
            LL_matrice, LR_matrice, MLR_matrice = get_lon_den_LLs(samples[(true_s2, init)], sampleSizes, compParameters,
                                                                  (llAlpha, llBeta), s2_grid, H)

            allLL_matrice[(true_s2, init)] = LL_matrice

            allLR_matrice[(true_s2, init)] = LR_matrice

            all_MLRs[(true_s2, init)] = MLR_matrice

        # save stuff
        with open(LLfilename, "wb") as LLoutfile:
            pickle.dump((allLL_matrice, allLR_matrice, all_MLRs), LLoutfile)

    # now do stat stuff
    import matplotlib.pyplot as plt

    # initialize containers
    true_s2s = []
    # shats = [];
    tildeS2s = []
    upperCIs = []
    lowerCIs = []
    medians = []
    means = []
    print('')
    for i,(true_s2, init) in enumerate(SNPs):
        # if it's neutral, skip
        if true_s2 == "0.0000": continue
        shat_pool = all_MLRs[(true_s2, init)]["joint_all"][:, 0]
        # stilde_pool = get_tildeS2s_from_hatS2s(re.sub("SFS", init, args.CIfile[0]), shat_pool)
        stilde_pool = shat_pool + (float(true_s2) - sMLEs[i])
        # shats.append(shat_pool)
        tildeS2s.append(stilde_pool)
        true_s2s.append(float(true_s2))
        upperCI, median, lowerCI = np.quantile(stilde_pool, [0.975, 0.5, 0.025])
        mn = np.mean(stilde_pool)
        # print quantile
        print(f'{init}, s2 = {true_s2}, [2.5% 50% 97.5%] = {(lowerCI, median, upperCI)}, mean = {mn}')
        upperCIs.append(upperCI)
        lowerCIs.append(lowerCI)
        medians.append(median)
        means.append(mn)

    # write out CIs
    CIfilename = f'sim{args.seed}_4SNPs-h{H:g}-bootStrapFreq{strNe}_tildeS2s_95CIs.txt'
    CImatrix = np.array([true_s2s, means, medians, lowerCIs, upperCIs]).transpose()
    np.savetxt(CIfilename, CImatrix, fmt="%g", delimiter="\t",
               header="true_s2\tshat_mean\tshat_median\tCI_lower\tCI_upper")

    # now plot violins
    # figtitle = 'Simulated s_MLE distributions'
    figname = f'sim{args.seed}_4SNPs-h{H:g}-bootStrapFreq{strNe}_tildeS2s_violins-ROC-presRecall.png'

    fig, axd = plt.subplot_mosaic([['top', 'top'], ['bottomL', 'bottomR']], figsize=(7, 6),
                                  gridspec_kw={'height_ratios': [1.5, 1]})  # , constrained_layout=True

    ## top: violins
    print(true_s2s)
    xCoord = np.arange(len(true_s2s))
    colors = plt.cm.Dark2(range(len(true_s2s)))[[2, 0, 1, 3], :]
    # plt.grid(visible = True)
    # violin_labels = [f'{rsIDs[i]}\nchr{chrs[i]}-{positions[i]}\n{true_s2s[i]: .3f}' for i in range(4)]
    violin_labels = [f's = {true_s2s[i]:.2f}\n{rsIDs[i]}' for i in range(4)]
    axd['top'].set_xticks(xCoord, violin_labels)
    # axd['top'].grid(visible = True)
    axd['top'].tick_params(labelsize=8)
    # axd['top'].set_xlabel ('True Selection Coefficient')
    axd['top'].set_ylabel('Unbiased MLEs')

    ### add scatters
    for (idx, thisData) in enumerate(tildeS2s):
        jitterX = np.random.normal(loc=xCoord[idx], scale=0.02, size=len(thisData))
        jitterY = np.random.normal(loc=thisData, scale=0.01)
        axd['top'].plot(jitterX, jitterY, color=colors[idx], marker="o", linestyle="", markersize=2, alpha=0.3)

    theQuantiles = [[0.025, 0.975]] * (len(true_s2s))
    axd['top'].violinplot(dataset=tildeS2s, positions=xCoord, showmedians=True, showextrema=False, quantiles=theQuantiles)

    axd['top'].set_ylim(-0.95, 0.95)
    axd['top'].set_xlim(min(xCoord) - 0.5, max(xCoord) + 0.5)
    # minSel, maxSel = s2_grid.min(), s2_grid.max()
    # axd['top'].hlines(y = [float(s) for s in trueS2_list], xmin = min(xCoord)-0.2, xmax = max(xCoord)+0.2, colors="0.6", linestyle = "--", linewidth=0.5) #
    # just some mins and maxes
    axd['top'].hlines(y=0, xmin=min(xCoord) - 0.5, xmax=max(xCoord) + 0.5, colors="black", linestyle="-",
                      linewidth=0.5)  # [0,1]
    axd['top'].hlines(y=true_s2s, xmin=min(xCoord) - 0.5, xmax=max(xCoord) + 0.5, colors=colors, linestyle="--",
                      linewidth=0.8)  # [0,1]
    # axd['top'].title.set_text(figtitle)
    axd['top'].legend(loc='best', ncol=3, facecolor=None, fontsize="small")

    ## down left: ROC
    # now plot ROC & prec
    ## page set up
    axd['bottomL'].title.set_text("ROC")
    axd['bottomL'].set_xlabel('False Positive Rate')
    axd['bottomL'].set_ylabel('True Positive Rate')
    axd['bottomR'].title.set_text("Precision-Recall")
    axd['bottomR'].set_xlabel('Recall')
    axd['bottomR'].set_ylabel('Precision')

    for thisAxe in (axd['bottomL'], axd['bottomR']):
        thisAxe.set_xlim([-0.04, 1.04])
        thisAxe.set_ylim([-0.04, 1.04])
        thisAxe.tick_params(labelsize=9)
        thisAxe.hlines([0, 1], -0.05, 1.05, "0.5", "-.", linewidth=0.5)
        thisAxe.vlines([0, 1], -0.05, 1.05, "0.5", "-.", linewidth=0.5)

    axd['bottomL'].axline((0, 0), (1, 1), linestyle=":", linewidth=0.5, color="0.5")
    # axd['bottomR'].axline((0,1), (1,0), linestyle = ":" , linewidth=0.5, color = "0.5")
    axd['bottomR'].hlines(0.5, -0.06, 1.06, linestyle=":", linewidth=0.8, color="red", alpha=0.8)

    s_idx = 0
    for true_s2, init in SNPs:

        if true_s2 == '0.0000': continue

        TP, FP, TN, FN = get_ROCs(all_MLRs[('0.0000', init)]["joint_all"][:, 1],
                                  all_MLRs[(true_s2, init)]["joint_all"][:, 1])

        TPR = np.array(TP) / (np.array(TP) + np.array(FN))
        FPR = np.array(FP) / (np.array(FP) + np.array(TN))

        # recalls = np.array(TP)/(np.array(TP) + np.array(FN))
        with np.errstate(divide='ignore', invalid='ignore'):
            precision = np.array(TP) / (np.array(TP) + np.array(FP))

        axd['bottomL'].plot(FPR, TPR, color=colors[s_idx], label=rsIDs[s_idx])
        axd['bottomR'].plot(TPR, precision, color=colors[s_idx], label=rsIDs[s_idx])
        s_idx += 1

    axd['bottomL'].legend(loc='center left', bbox_to_anchor=(1, 0.5), title="SNP", ncol=1, edgecolor=(1,1,1,0),
                          fontsize="small")

    fig.tight_layout()
    plt.savefig(figname, dpi=400)


if __name__ == '__main__':
    main()
