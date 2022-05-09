"""
Python module for HMM class to compute likelihood of temporal data under the Wright-Fisher diffusion.
"""
import numpy
import scipy.stats
import scipy.special


# local imports
# from . import diffusion_core
# from . import utility
# the python package management system is a pure delight
import pathlib
module_path = pathlib.Path(__file__).parent.resolve().__str__()
import sys
sys.path.insert(0, module_path)
import diffusion_core
import utility




EPSILON = 1e-12




# make a general class for power of square matrix
class MatrixPower:


    # @staticmethod
    # def highestPowerOfTwo (anInteger):
    #     assert (anInteger >= 1)
    #     aPower = -1
    #     while (anInteger >= 1):
    #         anInteger /= 2
    #         aPower += 1
    #     return (numpy.power (2, aPower), aPower)


    def __init__ (self, daMatrix):
        # should be square matrix, and perhps we want to remember the dimension
        assert (len(daMatrix.shape) == 2)
        assert (daMatrix.shape[0] == daMatrix.shape[1])
        self.matrixDimension = daMatrix.shape[0]
        # make a dict that remembers the powers we already computed
        self.powers = {}
        # power 0 is identiy
        self.powers[0] = numpy.eye (self.matrixDimension)
        # and power 1 is matrix itself
        self.powers[1] = daMatrix


    def get (self, p):
        # for now only positive integer powers
        assert (isinstance(p, (int, numpy.integer))), type(p)
        assert (p >= 0)

        # and then get the power
        if p in self.powers:
            matrix = self.powers[p]
        else:
            # see if we have half exponent
            halfexp = p//2
            if 2*halfexp not in self.powers:
                halfMatrix = self.get (halfexp)
                almost_matrix = numpy.dot (halfMatrix, halfMatrix)
                self.powers[halfexp*2] = almost_matrix
            else:
                almost_matrix = self.powers[2*halfexp]
            # if uneven power, we need to add one 
            if p%2 == 1:
                almost_matrix = numpy.dot (almost_matrix, self.get(1))
            # store it before returning
            matrix = almost_matrix
            self.powers[p] = matrix

        # return it
        return(matrix)


    # def get (self, p):
    #     # for now only positive integer powers
    #     assert (isinstance (p, (numpy.integer, int))), (type(p), p)
    #     assert (p >= 0)

    #     # and then get the power
    #     if p in self.powers:
    #         returnMatrix = self.powers[p]
    #     else:
    #         # now see what the highest power of two is
    #         (highestPower, exponent) = self.highestPowerOfTwo (p)
    #         if (highestPower == p):
    #             # it actually is a power of two not yet computed, so compute it
    #             halfexp = p//2
    #             halfMatrix = self.get (halfexp)
    #             # multiply the half exponents
    #             returnMatrix = numpy.matmul (halfMatrix, halfMatrix)
    #         elif (highestPower < p):
    #             # there is something left after the power of two, so get two things and compute
    #             powerMatrix = self.get (highestPower)
    #             restMatrix = self.get (p - highestPower)
    #             # and multiply stuff
    #             returnMatrix = numpy.matmul (powerMatrix, restMatrix)
    #         else:
    #             assert (False)
    #         # store it before returning
    #         self.powers[p] = returnMatrix

    #     # return it
    #     return returnMatrix


class SelectionGrid:


    # initializes the object to efficiently store the transition matrices for the given selection coefficient list
    def __init__ (self, Ne, s1s, s2s, mAlpha, mBeta, selGridSize=11, numStates=1001):

        self.Ne = Ne
        self.mAlpha = mAlpha
        self.mBeta = mBeta

        # the states
        self.numStates = numStates
        yGrid = utility.getLinearGrid (numStates)
        yBoundaries = utility.getLinearBounds (yGrid)

        # check a bunch of stuff
        assert (isinstance (s1s, (list, numpy.ndarray)))
        s1s = numpy.array (s1s)
        assert (isinstance (s2s, (list, numpy.ndarray)))
        s2s = numpy.array (s2s)
        assert (len(s1s.shape) == 1)
        assert (s1s.shape == s2s.shape)

        # now build the scaffold
        self.sOneBounds = (min(s1s), max(s1s))
        self.sTwoBounds = (min(s2s), max(s2s))

        # just linear for now (cause a bit easier with + and -)
        # TODO maybe something smarter in the future
        self.sOneGrid = numpy.linspace (self.sOneBounds[0], self.sOneBounds[1], selGridSize)
        self.sTwoGrid = numpy.linspace (self.sTwoBounds[0], self.sTwoBounds[1], selGridSize)

        # what combinations do we really need?
        allPairs = []
        for sIdx in numpy.arange (len(s1s)):
            # where on grid do they fall
            sOneIdx = (numpy.abs(self.sOneGrid - s1s[sIdx])).argmin()
            sTwoIdx = (numpy.abs(self.sTwoGrid - s2s[sIdx])).argmin()
            allPairs.append ((sOneIdx, sTwoIdx))
        self.idxToPair = list(set(allPairs))

        # make a map
        self.pairToIdx = {}
        for (theIdx, thePair) in enumerate(self.idxToPair):
            self.pairToIdx[thePair] = theIdx

        self.transitionMatrices = []
        # build all the matrix objects we need
        for (sOneIdx, sTwoIdx) in self.idxToPair:
            # get a single step matrix
            singleStep = SelHmm.TransitionProabilities.normalMatrix (1.0, yGrid, yBoundaries,
                                    self.Ne, self.sOneGrid[sOneIdx], self.sTwoGrid[sTwoIdx], self.mAlpha, self.mBeta)
            # and then remember this matrix power object
            self.transitionMatrices.append (MatrixPower (singleStep))
        assert (len(self.transitionMatrices) == len(self.idxToPair))

        # that should be everything in terms of setting up


    # this one puts together transition matrices to span the given selection coeficients
    def getTransitionMatrix (self, s1s, s2s):

        assert (isinstance (s1s, (list, numpy.ndarray)))
        s1s = numpy.array (s1s)
        assert (isinstance (s2s, (list, numpy.ndarray)))
        s2s = numpy.array (s2s)
        assert (len(s1s.shape) == 1)
        assert (s1s.shape == s2s.shape)

        # make a list where we indicate the index of the right selection object in each geeneration
        selectionIndex = numpy.zeros(len(s1s), dtype=int)
        for sIdx in numpy.arange (len(s1s)):
            sOneIdx = (numpy.abs(self.sOneGrid - s1s[sIdx])).argmin()
            sTwoIdx = (numpy.abs(self.sTwoGrid - s2s[sIdx])).argmin()
            thisPair = (sOneIdx, sTwoIdx)
            selectionIndex[sIdx] = self.pairToIdx[thisPair]

        # put together a matrix that makes this transition
        currentT = 0
        tEnd = len(s1s)
        returnMatrix = numpy.identity (self.transitionMatrices[0].matrixDimension)
        # this should all work properly, even if the list of selection coefficients is empty
        while (currentT < tEnd):

            # find how long we can go with the same selection object
            restIdxs = numpy.arange(tEnd - currentT)
            restSelection = selectionIndex[currentT:tEnd]
            differentRest = restIdxs[restSelection != restSelection[0]]
            
            if (len(differentRest) > 0):
                nextT = currentT + differentRest[0]
            else:
                nextT = tEnd

            stepMatrix = self.transitionMatrices[selectionIndex[currentT]].get (nextT - currentT)
            returnMatrix = numpy.matmul (returnMatrix, stepMatrix)

            # step as far as you can
            currentT = nextT

        return returnMatrix


# class for emission & transition matrices, initial distribution, and HMM computation
class SelHmm:
    """
    Objects of this class can be used to compute the likelihood of observing temporal allele frequency data under given parameters.
    """

    # make a dirac delta
    @staticmethod
    def discreteDiracDeltaAt (freq, yGrid):
        # this is to return
        p = numpy.zeros (len(yGrid))
        # where to put the mass
        diffs = numpy.abs (freq - yGrid)
        minIdx = numpy.argmin (diffs)
        if (numpy.isclose(diffs[minIdx], 0)):
            # only one
            p[minIdx] = 1
        else:
            # shoud be inbetween two, find other one
            minDiff = diffs[minIdx]
            # all diffs should be between 0 and 1 anyways
            diffs[minIdx] = 1e10
            nextIdx = numpy.argmin (diffs)
            nextDiff = diffs[nextIdx]
            # give some to both, has to be crossed like this so lower diff gives higher weight
            p[minIdx] = nextDiff / (minDiff + nextDiff)
            p[nextIdx] = minDiff / (minDiff + nextDiff)
        # that should be good
        return p
        

    # initial distribution with given initial frequency
    @staticmethod
    def givenInitFrequency (initFreq, yGrid):
        p = numpy.zeros (len(yGrid))
        # where to put the one?
        diffs = numpy.abs (initFreq - yGrid)
        # too how many hidden states do we belong?
        # THIS IS MOSTLY JUST TO SPLIT TIES WHEN WE EXACTLY ON THE BOUNDARY
        numStates = numpy.sum(diffs == diffs.min())
        p[diffs == diffs.min()] = 1/numStates
        # for now, no starting in the boundary, cause we have to fix the transition probabilities for that
        # assert (not numpy.isclose (p[0], 1) and not numpy.isclose (p[-1], 1)), "no starting in the boundary, maybe initFreq to small for given # of states"
        return p


    # inner interface for different emission probabilities
    class EmissionProbabilities:


        def checkSamples (self, samplesSizeMatrix, sampleMatrix):
            assert (False), "Derived class does not implement full interface."


        def getEmissionProbs (self, sampleSize, numAlleles):
            assert (False), "Derived class does not implement full interface."


    # inner class for the integer emission probabilities
    class IntegerEmissionProbabilities (EmissionProbabilities):


        # initialize object and precompute emission matrices for all sample sizes given
        def __init__ (self, yGrid, sampleSizesSet):
            assert (isinstance (sampleSizesSet, set))
            # remember the grid
            self.yGrid = yGrid
            self.numStates = len(self.yGrid)
            # and the sample sizes that we will encounter
            self.sampleSizesSet = sampleSizesSet
            # set up the object
            self.emissionMatrices = {}
            for n in self.sampleSizesSet:
                assert (isinstance (n, (int, numpy.integer))), "Only integer sample sizes allowed."
                assert (n >= 0), "Only non-negative sample sizes allowed."
                # get a matrix for this sample size
                daMatrix = numpy.zeros ((len(yGrid),n+1))
                for k in range(n+1):
                    daMatrix[:,k] = scipy.stats.binom.pmf (k, n, yGrid)
                # save this one
                self.emissionMatrices[n] = daMatrix


        def checkSamples (self, samplesSizeMatrix, sampleMatrix):
            assert (sampleMatrix.dtype == int), "Only integer samples allowed."
            assert (samplesSizeMatrix.dtype == int), "only integer sample sizes allowed."
            keySet = set(self.emissionMatrices.keys())
            sampleSizeSet = set(samplesSizeMatrix.flatten())
            assert (sampleSizeSet.issubset(keySet)), f"Did not pre-compute emissions for sample sizes: {sampleSizeSet.difference(keySet)}."
            return True


        # return the row of the EmissionProb matrix for the given sample size corresponding to numAlleles emitted
        def getEmissionProbs (self, sampleSize, numAlleles):
            return self.emissionMatrices[sampleSize][:,numAlleles]


    # inner class for the fractional emission probabilities
    class FractionalEmissionProbabilities (EmissionProbabilities):

        def __init__ (self, globalYGrid):
            # this is a slight hack for freqs zero, but probably works and gives a 30% speed increase
            logZero = -1e12
            largeNumber = 2
            assert (globalYGrid.max() <= 1)
            assert (globalYGrid.min() >= 0)

            # remember the grid
            # self.yGrid = yGrid
            yGrid = globalYGrid.copy()
            oneMYGrid = 1 - yGrid

            # need to work around zeros, becuase sometimes log(0) is not wanted
            zeroIdx = numpy.isclose (yGrid, 0)
            yGrid[zeroIdx] = largeNumber
            logYGrid = (numpy.log(yGrid)).clip(min=logZero)
            logYGrid[zeroIdx] = logZero
            self.logYGrid = logYGrid

            zeroIdx = numpy.isclose (oneMYGrid, 0) 
            oneMYGrid[zeroIdx] = largeNumber
            oneMLogYGrid = (numpy.log(oneMYGrid)).clip(min=logZero)
            oneMLogYGrid[zeroIdx] = logZero
            self.oneMLogYGrid = oneMLogYGrid


        def checkSamples (self, samplesSizeMatrix, sampleMatrix):
            assert (sampleMatrix.dtype in [int, float]), "Only fractional samples allowed."
            assert (samplesSizeMatrix.dtype in [int, float]), "Only fractional samples sizes allowed."
            return True


        # get emission probabilities for all hidden states for observing numAlleles focal alleles in a sample of size sampleSize
        # directly computed using the binomial
        def getEmissionProbs (self, sampleSize, numAlleles):
            # we might want to do a normal approximation for large sample sizes
            # and this could lead into all kinds of numerical trouble
            # we might actually just use this
            # it's probably stable enough, and I don' think we can really get it faster with some fancy interpolation
            # return scipy.special.binom (sampleSize, numAlleles) * numpy.power (self.yGrid, numAlleles) * numpy.power (1-self.yGrid, sampleSize - numAlleles)
            assert (numAlleles >= 0)
            assert (sampleSize >= numAlleles)
            return numpy.exp (numpy.log (scipy.special.binom (sampleSize, numAlleles)) + numAlleles * self.logYGrid + (sampleSize - numAlleles) * self.oneMLogYGrid)



    # inner base class for transition probabilities
    # if we have a general interface, we can replace what's under the hood
    class TransitionProabilities:


        # take a normal step for deltaT generations (can be fractional)
        @staticmethod
        def normalMatrix (deltaT, yGrid, boundaries, Ne, s1, s2, mAlpha, mBeta):
            # get the parameters for the gaussian

            # put infities into boundaries
            theseBoundaries = boundaries.copy()
            assert (numpy.isclose(theseBoundaries[0], 0))
            theseBoundaries[0] = float("-inf")
            assert (numpy.isclose(theseBoundaries[-1], 1))
            theseBoundaries[-1] = float("inf")
            # using full grid should be fine
            assert (len(yGrid)+1 == len(boundaries))
    
            daMu = diffusion_core.mu (yGrid, s1, s2, mAlpha, mBeta)
            # some safety check for now that hopefully catches stuff we can't do yet
            # again,should be ok, but check anyways
            assert (numpy.logical_not(numpy.isin (daMu, [float("-inf"),float("inf")])).all())

            daSigmaSq = diffusion_core.sigmaSq (yGrid, Ne)
            assert (len(daMu) == len(daSigmaSq))

            (daLoc, daScale) = diffusion_core.diffusionLocScale (yGrid, daMu, daSigmaSq, deltaT)

            # now we have da means and da variances, so get da gaussians
            numStates = len(yGrid)
            oneStepMatrix = numpy.zeros ((numStates,numStates))
            assert (len(daMu) == numStates)
            # fill the interior
            for idx in numpy.arange(len(daMu)):
                # preRow = numpy.concatenate (([0],scipy.stats.norm.cdf (theseBoundaries, daLoc[idx], daScale[idx]),[1]))
                if (numpy.isclose (daScale[idx], 0)):
                    # 0 variance, so deterministic step
                    # but where to?
                    oneStepMatrix[idx,:] = SelHmm.discreteDiracDeltaAt (daLoc[idx], yGrid)
                else:
                    # positive variance, so gaussian step
                    preRow = scipy.stats.norm.cdf (theseBoundaries, daLoc[idx], daScale[idx])
                    oneStepMatrix[idx,:] = numpy.diff (preRow) 

            # give it away now
            return oneStepMatrix


        def transitionMatrix (self, tStart, tEnd):
            assert (False), "Intance of ase class not allowed."


    # implementation of transition class for constant selection coefficient
    class ConstantTransitionProabilities(TransitionProabilities):


        def __init__(self, yGrid, boundaries, Ne, s1, s2, mAlpha, mBeta, deltaT=1):
            self.deltaT = deltaT

            # get a single step matrix for self.deltaT
            singleStep = self.normalMatrix (self.deltaT, yGrid, boundaries, Ne, s1, s2, mAlpha, mBeta)

            # and then set up the matrix power object
            self.oneStepPower = MatrixPower (singleStep)


        def transitionMatrix (self, tStart, tEnd):
            # get the probability matrix for transitioning from generation tStart to generation tEnd

            # first get the right power when self.deltaT is the stepsize for the one step matrix
            daPower = round((tEnd-tStart)/self.deltaT)

            # now we should have the right interger power, so raise the matrix
            return self.oneStepPower.get (int(daPower))
    

    # implementation of transition class for piece-wise constant selection coefficient
    class PiecewiseTransitionProabilities(TransitionProabilities):


        # initialize class
        # selectionChangeTimes is a list/np.array of >=1 numbers indicating the generations when (s1,s2) changes
        # s1s & s2s are list/np.arrays of length (len(selectionChangeTimes) + 1) specifying (s1,s2) for the respective generations
        def __init__(self, yGrid, boundaries, Ne, s1s, s2s, mAlpha, mBeta, selectionChangeTimes, deltaT=1):

            # some checks
            assert (type(selectionChangeTimes) in [list, numpy.ndarray])
            assert (len(s1s)-1 == len(selectionChangeTimes))
            assert (len(s2s)-1 == len(selectionChangeTimes))
            assert (all ([x <= y for (x,y) in zip (selectionChangeTimes[:-1], selectionChangeTimes[1:])]))

            # remember stuff
            self.deltaT = deltaT

            # get a matrix power for every unqiue selection coefficient
            # add infinities for convenience
            self.realSelectionChangeTimes = numpy.concatenate (([-float("inf")], selectionChangeTimes, [float("inf")]))
            self.transitionMatrices = []
            for sIdx in numpy.arange(len(s1s)):
                # get a single step matrix for self.deltaT
                singleStep = self.normalMatrix (self.deltaT, yGrid, boundaries, Ne, s1s[sIdx], s2s[sIdx], mAlpha, mBeta)

                # and then remember this matrix power object
                self.transitionMatrices.append (MatrixPower (singleStep))
            assert (len(self.transitionMatrices)+1 == len(self.realSelectionChangeTimes))



        def transitionMatrix (self, tStart, tEnd):
            # get the probability matrix for transitioning from generation tStart to generation tEnd
            # this could be several pieces

            # first find where we are
            pieces = numpy.searchsorted (self.realSelectionChangeTimes, [tStart, tEnd])
            assert (len(pieces) == 2)
            startPieceIdx = pieces[0] - 1
            endPieceIdx = pieces[1] - 1

            # some loop
            thisPieceIdx = startPieceIdx
            returnMatrix = None
            assert (startPieceIdx <= endPieceIdx)
            # plus one to have base case included (and does at lest one loop)
            while (thisPieceIdx < endPieceIdx + 1):
                # get the time in this interval
                realStartTime = max (self.realSelectionChangeTimes[thisPieceIdx], tStart)
                realEndTime = min (self.realSelectionChangeTimes[thisPieceIdx+1], tEnd)

                # get the right power in this interval when self.deltaT is the stepsize for the one step matrix
                daPower = round((realEndTime-realStartTime)/self.deltaT)

                # now we should have the right interger power, so raise the matrix
                thisMatrix = self.transitionMatrices[thisPieceIdx].get (int(daPower))
                if (returnMatrix is None):
                    returnMatrix = thisMatrix
                else:
                    returnMatrix = numpy.matmul (returnMatrix, thisMatrix)

                # update index
                thisPieceIdx += 1
            
            # and then all should be good
            assert (returnMatrix is not None)
            return returnMatrix
    

    # implementation of transition class for continuously varying selection coefficients
    # this class computes all the transition matrices for each generation
    class ContinuousTransitionProabilitiesAll(TransitionProabilities):


        # initialize object, that is, compute all the single step matrices for every generation
        def __init__(self, yGrid, boundaries, Ne, s1, s2, mAlpha, mBeta):

            assert (type(s1) in [list, numpy.ndarray])
            assert (type(s2) in [list, numpy.ndarray])
            assert (len(numpy.array(s1).shape) == 1)
            assert (len(s1) == len(s2))

            # prepare for future use
            self.previousTransitions = {}

            # print ("== prepare single steps")
            # get all single step matrices
            self.singleStepMatrices = []
            for sIdx in numpy.arange(len(s1)):
                self.singleStepMatrices.append (self.normalMatrix (1.0, yGrid, boundaries, Ne, s1[sIdx], s2[sIdx], mAlpha, mBeta))
                if (sIdx % 16 == 0):
                    pass
                    # print (sIdx)

            # should be all good
            assert (len(self.singleStepMatrices) == len(s1))


        # check if we have single step matrices for each generation, and that first time is 0
        def timesCompatible (self, times):
            assert (numpy.isclose (times[0], 0))
            # we should have one for every generation
            assert (numpy.isclose (times[-1], len(self.singleStepMatrices))), (times[-1], len(self.singleStepMatrices))
            # if we here, all good
            return True


        def transitionMatrix (self, tStart, tEnd):

            # discretize them so we can use them as indices
            tStart = int(tStart)
            tEnd = int(tEnd)

            # where we here before?
            if ((tStart, tEnd) in self.previousTransitions):
                return self.previousTransitions[(tStart, tEnd)]

            # if not, do normal stuff
            assert (tStart <= tEnd)

            # put together a matrix that makes this transition
            currentT = tStart
            # print (f"++ ({currentT}, {tEnd})")
            returnMatrix = numpy.identity (self.singleStepMatrices[0].shape[0])
            while (currentT < tEnd):
                # print (currentT)
                returnMatrix = numpy.matmul (returnMatrix, self.singleStepMatrices[currentT])
                currentT += 1
            
            # remember that we were here before
            self.previousTransitions[(tStart, tEnd)] = returnMatrix

            return returnMatrix

        
    # implementation of transition class for continuously varying selection coefficients
    # this class uses the selectionGridObject for the transitions, so it does not have to compute them at each step
    class ContinuousTransitionProabilitiesZero(TransitionProabilities):


        def __init__(self, s1, s2, selectionGridObject):

            assert (type(s1) in [list, numpy.ndarray])
            assert (type(s2) in [list, numpy.ndarray])
            assert (len(numpy.array(s1).shape) == 1)
            assert (len(s1) == len(s2))

            # I think that for now we just store the stuff
            self.selectionGridObject = selectionGridObject
            self.s1 = s1
            self.s2 = s2

            # its ok to remember some stuff
            self.previousTransitions = {}


        # check if we have selection for each generation, and that first time is 0
        def timesCompatible (self, times):
            assert (numpy.isclose (times[0], 0))
            # we should have one for every generation
            assert (numpy.isclose (times[-1], len(self.s1))), (times[-1], len(self.s1))
            # if we here, all good
            return True


        def transitionMatrix (self, tStart, tEnd):

            # discretize them so we can use them as indices
            tStart = int(tStart)
            tEnd = int(tEnd)

            # total ts are checked somewhere else
            assert (tStart >= 0)
            assert (tEnd <= len(self.s1))

            # where we here before?
            if ((tStart, tEnd) in self.previousTransitions):
                return self.previousTransitions[(tStart, tEnd)]

            # if not, do normal stuff
            assert (tStart <= tEnd)

            # what is the exact selection sequence that we have here?
            theseS1 = self.s1[tStart:tEnd]
            theseS2 = self.s2[tStart:tEnd]

            # get a matrix out of the grid object just for these coefficients
            returnMatrix = self.selectionGridObject.getTransitionMatrix (theseS1, theseS2)

            # remember that we were here before
            self.previousTransitions[(tStart, tEnd)] = returnMatrix

            return returnMatrix


    # parse initial condition and return corresponding initial distribution
    @staticmethod
    def getInitDistribution (yGrid, gridWeights, Ne, s1, s2, mAlpha, mBeta, initCond,
                                initFreq=None, initMAlpha=None, initMBeta=None, initS1=None, initS2=None):

        # if special parameters for initial distribution given, use those
        if (initMAlpha is None):
            initMAlpha = mAlpha
        if (initMBeta is None):
            initMBeta = mBeta
        if (initS1 is None):
            initS1 = s1
        if (initS2 is None):
            initS2 = s2

        # these need non-zero mutation rates
        if (initCond in ["statBeta", "statBetaSel"]):
            assert (initMAlpha > 0), "Stationary initial condition needs positive mutation rate. Or supply one just for initial with 'initMAlpha'"
            assert (initMAlpha > 0), "Stationary initial condition needs positive mutation rate. Or supply one just for initial with 'initMBeta'"

        # handle initial distributions
        if (initCond == "initFreq"):

            assert (initFreq is not None)

            # do we have one frequency or many frequencies
            if ((type(initFreq) == float) or (len(initFreq) == 1)):
                # one 
                if (type(initFreq) != float):
                    initFreq = initFreq[0]
                assert (initFreq >= 0)
                assert (initFreq <= 1)

                # have given init frequency
                initial = SelHmm.discreteDiracDeltaAt (initFreq, yGrid)
            else:
                # many
                assert (len(initFreq) > 1)
                initFreq = numpy.array(initFreq)
                assert ((initFreq >= 0).all())
                assert ((initFreq <= 1).all())

                # already good
                initial = numpy.zeros ((len(initFreq),len(yGrid)))
                for (idx, f) in enumerate(initFreq):
                    initial[idx,:] = SelHmm.discreteDiracDeltaAt (f, yGrid)

        elif (initCond == "statBeta"):
            initial = diffusion_core.Stationary.stationaryNeutralBeta (Ne, initMAlpha, initMBeta, yGrid, gridWeights)
            #print(self.initDist)
        elif (initCond == "statBetaSel"):
            initial = diffusion_core.Stationary.stationarySelectionBeta (Ne, initMAlpha, initMBeta, s1, s2, yGrid, gridWeights)
        # elif (initCond == "statNeutralPRF"):
        # elif (initCond == "statSelPRF"):
        # easy enough to have it here
        elif (initCond == "uniform"):
            initial = gridWeights
        else:
            # unknown initial conditions
            assert (False), f'Unknown initial condition {initCond}. Please choose between "statBeta", "statBetaSel", "uniform", and "initFreq"'

        return initial


    # initialize emission & transition matrices, initial distribution, and set up the HMM
    # TODO kinda weird to have sample sizes in __init__ and computeLogLikelihood, but we deal with it for now
    def __init__ (self, Ne, s1, s2, mAlpha, mBeta, initCond,
                    initFreq=None, initMAlpha=None, initMBeta=None, initS1=None, initS2=None, sampleSizesSet=None,
                    numStates=1001, deltaT=1, emissionType="integer", transitionType="constant", selectionChangeTimes=None,
                    selectionGridObject=None):
        r"""Construct object of class `SelHmm`.

        Initializes object to compute likelihood of observing temporal allele frequency data under the given parameters.

        Parameters
        ----------
        Ne : int or float
            Effective diploid population size (*i.e.*, for each locus, 2Ne alleles exist in total).

        s1 : float or array_like
            Selection coefficients of the heterozygote.

        s2 : float or array_like
            Selection coefficients of the homozygote.

        mAlpha : float
            Per-site per-generation forward mutation rate.

        mBeta : float
            Per-site per-generation backward mutation rate.

        initCond : {"uniform", "initFreq", "statBeta", "statBetaSel"}
            Specify the initial condition for the HMM at generation zero.
            Depends on the type of initial condition, other parameters need to be provided:
            - "uniform" :
                HMM starts with a uniform distribution on [0,1].
            - "initFreq" :
                Starts with a given allele frequency. Must also specify the frequency with `initFreq`.
            - "statBeta" :
                Starts with the stationary Beta distribution under mutation rates `mAlpha` and `mBeta`.
            - "statBetaSel" :
                Starts with the stationary Beta distribution under mutation rate `initAlpha`, `initBeta`,
                and selection coefficient `initS1`, `initS2`.

        emissionType : {"integer", "fractional"}, default="integer"
            Type of emission data, *i.e.* observed sample numbers, to consider. Choose "integer" for
            integer allele counts, "fractional" for when such count has been adjusted or a fractional estimate.

        transitionType : {"constant", "piecewise"}, default="constant"
            Type of transitions to consider in the HMM:
            - "constant" :
                Selection coefficients stay constant throughout the entire duration considered.
                `s1` and `s2` must be constants (not `array_like`) for this option.
            - "piecewise" :
                The entire duration can be considered as several pieces in tandem, where each piece has a
                different pair of selection coefficients. Must also specify `selectionChangeTimes` for this option.

        Other Parameters
        ----------------
        deltaT : int or float, optional, default=1
            Unit increment of time (in generations).

        numStates : int or float, optional, default=1001
            Number of discretized states with which to discretize the allele frequency space [0,1].

        initFreq : float, optional
            Required when ``initCond="initFreq"``. Must be between 0 and 1.

        initAlpha : float, optional, default=`mAlpha`
            Per-site per-generation forward mutation rate underlying the initial Beta distribution.

        initBeta : float, optional, default=`mBeta`
            Per-site per-generation backward mutation rate underlying the initial Beta distribution.

        sampleSizesSet : set of array_like
            Set of list, numpy array, or tuple objects of the same length as the number of sampling time points
            that summarizes all possible sample sizes

        selectionChangeTimes : int or array_like, optional
            Set the generation time when selection coefficients change.
            Must match the length of `allS1` and `allS2`.
        """
        # the states
        self.numStates = numStates
        yGrid = utility.getLinearGrid (numStates)
        boundaries = utility.getLinearBounds (yGrid)
        gridWeights = utility.getWeights (boundaries) 

        # initial distribution
        self.initDist = SelHmm.getInitDistribution (yGrid, gridWeights, Ne, s1, s2, mAlpha, mBeta, initCond,
                                                    initFreq=initFreq, initMAlpha=initMAlpha, initMBeta=initMBeta, initS1=initS1, initS2=initS2)

        # emissions for all possible sizes
        if (emissionType == "integer"):
            assert (sampleSizesSet is not None), "If integer emission, have to provide set of sample sizes to expect."
            assert (isinstance (sampleSizesSet, set))
            self.emissionObject = SelHmm.IntegerEmissionProbabilities (yGrid, sampleSizesSet)
        elif (emissionType == "fractional"):
            assert (sampleSizesSet is None)
            self.emissionObject = SelHmm.FractionalEmissionProbabilities (yGrid)
        else:
            assert (False), f"Unknown emission type: {emissionType}. Use 'integer' or 'fractional'."

        # transition
        if (transitionType == "constant"):
            assert (selectionChangeTimes is None)
            assert (selectionGridObject is None)
            assert (isinstance (s1, (int, numpy.integer, float, numpy.floating))), type(s1)
            assert (isinstance (s2, (int, numpy.integer, float, numpy.floating))), type(s2)
            self.transitionObject = SelHmm.ConstantTransitionProabilities (yGrid, boundaries, Ne, s1, s2, mAlpha, mBeta)
        elif (transitionType == "piecewise"):
            assert (selectionChangeTimes is not None)
            assert (selectionGridObject is None)
            assert (type(s1) in [list, numpy.ndarray])
            assert (type(s2) in [list, numpy.ndarray])
            self.transitionObject = SelHmm.PiecewiseTransitionProabilities (yGrid, boundaries, Ne, s1, s2, mAlpha, mBeta, selectionChangeTimes)
        elif (transitionType == "continuous"):
            # we don't want this
            assert (selectionChangeTimes is None)
            # but this
            assert (selectionGridObject is not None)
            # we can at least do this
            assert (selectionGridObject.numStates == self.numStates)
            assert (selectionGridObject.Ne == Ne)
            assert (selectionGridObject.mAlpha == mAlpha)
            assert (selectionGridObject.mBeta == mBeta)
            # but we want this
            assert (type(s1) in [list, numpy.ndarray])
            assert (type(s2) in [list, numpy.ndarray])
            # and this is ugly, but cannot help it for now
            assert (numpy.isclose (deltaT, 1))
            # self.transitionObject = SelHmm.ContinuousTransitionProabilitiesAll (yGrid, boundaries, Ne, s1, s2, mAlpha, mBeta)
            self.transitionObject = SelHmm.ContinuousTransitionProabilitiesZero (s1, s2, selectionGridObject)
        else:
            assert (False), f"Unknown transition type: {transitionType}. Use one of ['constant', 'piecewise', 'continuous']."

    # run the HMM & compute likelihoods for given data. Return (list of) log-likelihood
    # @numba.jit(nopython=True)
    def computeLogLikelihood (self, times, samplesSizes, samples):
        r"""Compute the likelihood of the given data under the HMM model.

        After specifying the underlying parameters of a `SelHmm` object,
        this function compute the log-likelihoods of observing the given samples
        at the sampling times `times`.
        Note that this model assumes all loci are independent and bi-allelic.

        Parameters
        ----------
        times : array_like
            The generation times when the samples were taken. Must start with zero and ascend from fast to present.

        samplesSizes : array_like
            An N by K matrix of the total numbers of alleles observed, *i.e.* sample sizes.
            N --> number of loci ;  K --> number of sampling times.
            ``sampleSizes[i,j]`` records the sample size of locus ``i`` at time point ``j``.

        samples : array_like
            An N by K matrix of the numbers of alleles.
            N --> number of loci ;  K --> number of sampling times.
            ``samples[i,j]`` records the number of a particular allele on locus ``i`` at time point ``j``.

        Returns
        -------
        numpy.array of log-likelihoods for the loci.
        """
        # get them in good shape
        samplesSizes = numpy.array(samplesSizes)
        samples = numpy.array(samples)

        # and check stuff
        assert (samplesSizes.shape[0] ==  samples.shape[0])
        assert (len(samplesSizes.shape) == len(samples.shape))
        if (len(samples.shape) == 1):
            # need to fake some stuff if only one replicate, to be compatible with the algorithm
            samplesSizes = numpy.array([samplesSizes])
            samples = numpy.array([samples])
        assert (len(samplesSizes.shape) == 2)
        assert (samplesSizes.shape[1] ==  samples.shape[1])
        assert (len(times)-1 == samplesSizes.shape[1])
        assert (self.emissionObject.checkSamples (samplesSizes, samples))

        numReplicates = samplesSizes.shape[0]

        # need some checking for possibly continuous selection
        # if (isinstance (self.transitionObject, SelHmm.ContinuousTransitionProabilitiesAll)):
        if (isinstance (self.transitionObject, SelHmm.ContinuousTransitionProabilitiesZero)):
            # check whether times are compatible with number of coefficients we have
            assert (self.transitionObject.timesCompatible (times))

        # set up container for scaling factors
        scalingFactors = numpy.zeros ((numReplicates, len(times)))

        # run the HMM (forward)
        for i in range(len(times)):
            # first iteration or later
            if (i == 0):
                # do the init
                if (len(self.initDist.shape) == 1):
                    # we just have one vector, so use it for all samples
                    currProbs = numpy.array([self.initDist,]*numReplicates)
                else:
                    # initDist already as some stuff for multiple samples
                    currProbs = self.initDist
                #print(f't=0, currProbs.shape = {currProbs.shape}')
                # this is what we ultimately need
                assert (currProbs.shape == (numReplicates, self.numStates)), f"Num init freqs might not match numReplicates: {currProbs.shape} == ({numReplicates}, {self.numStates})."
            else:
                # do the transition
                currProbs = numpy.matmul (currProbs, self.transitionObject.transitionMatrix (times[i-1], times[i]))

                # and then some emission
                sampleIdx = i-1
                thisEmissions = numpy.zeros ((numReplicates, self.numStates))
                for nIdx in range(numReplicates):
                    thisEmissions[nIdx,:] = self.emissionObject.getEmissionProbs (samplesSizes[nIdx,sampleIdx], samples[nIdx,sampleIdx])
                currProbs = thisEmissions * currProbs
                # print (elapsedTime)

            # normalize to sum to one and store the normalizing factor
            scalingFactors[:,i] = numpy.sum(currProbs,axis=1)
            currProbs = currProbs / scalingFactors[:,i].reshape((numReplicates,1))

        # these are the final likelihoods
        ll = numpy.sum (numpy.log(scalingFactors),axis=1)
        # however, any sample that had a scaling factor of 0 somewhere should get a - infty for ll, because prob = 0
        zeroProbsIndicator = (numpy.sum(scalingFactors <= 0, axis=1) >= 1)
        ll[zeroProbsIndicator] = float("-inf")
        return ll
