# Python library for core functions to compute HMM for WF diploid selection
import numpy
import scipy.special
import scipy.stats




# genetic drift term in the diffusion
def sigmaSq (y, Ne):
    return y * (1-y) / (2*Ne)


# infinitesimal change in mean under diploid selection of with coefficients (s1, s2)
# mAlpha/mBeta are 0->1/1->0 mutation probability (per site per generation)
def mu (y, s1, s2, mAlpha, mBeta):
    # general diploid selection pressure
    return (y*(1-y) * (s1*(1-2*y) + s2*y)) + (1-y)*mAlpha - y*mBeta


# mean fitness of population if allele freq is y
def meanScaledFitness (y, s1, s2, Ne):
    return 2*y*(1-y)*2*Ne*s1 + y*y*2*Ne*s2


# integral over the scale function in the Poisson Random Field model (no recurrent mutation)
def intScalePRF (y, s1, s2, Ne):
    if (numpy.isclose (s1, 0) and numpy.isclose (s2, 0)):
        # neutral case
        return 1-y
    if numpy.isclose (2*s1, s2):
        # genic selection
        # return the simplified expression
        return -(numpy.exp(-2*Ne*s2) - numpy.exp(-2*Ne*s2*y))/(2*Ne*s2)
    else:
        # general diploid selection
        # factor in the front
        factor1 = numpy.sqrt(numpy.pi/2)/(2*numpy.sqrt(Ne*complex(2*s1 - s2)))
        factor2 = numpy.exp(-2*Ne*s1*s1/(2*s1 - s2))
        factor = factor1 * factor2
        # first erfi
        arg1 = numpy.sqrt(2*Ne)*(s1-s2)/numpy.sqrt(complex(2*s1 - s2))
        erfi1 = scipy.special.erfi (arg1)
        # second erfi
        arg2 = numpy.sqrt(2*Ne)*((2*y-1)*s1-y*s2)/numpy.sqrt(complex(2*s1 - s2))
        erfi2 = scipy.special.erfi (arg2)
        # give it away now
        result = factor * (erfi1 - erfi2)
        assert (numpy.isclose(0,numpy.imag(result)).all())
        return numpy.real(result)


# return the mean & std of the normal distribution for the next "step" in the diffusion
def diffusionLocScale (currPos, daMu, daSigmaSq, deltaT):
    # this is how it works
    return (currPos + deltaT * daMu, numpy.sqrt(deltaT * daSigmaSq))


# collect the stationary distributions
class Stationary:


    # this normalizes stuff to one and adds the last mass on the left side
    @staticmethod
    def normalizeAndAddLeft (p, weights):
        assert (len(p) == len(weights)-1)
        # TODO maybe some epsilon in the future
        assert (numpy.min (p) >= 0)
        # now, the rest can either sum to something more than one or less than one
        daSum = numpy.sum (p)

        if (daSum <= 1):
            # it's less then one, so just add 1 - the sum to the left boundary
            return numpy.concatenate (([1-daSum], p))
        else:
            # so mass in the left boundary is just what's right next to it
            leftBoundary = p[0] * weights[0] / weights[1]
            preP = numpy.concatenate (([leftBoundary], p))
            # and then renormalize it
            return preP / numpy.sum(preP)

    # neutral
    # for now just evaluate f, maybe in future we can do integrals
    @staticmethod
    def stationaryNeutralPRF (Ne, m, yGrid, weights):
        assert (numpy.isclose (yGrid[0], 0))
        assert (numpy.isclose (yGrid[-1], 1))
        # all the mass missing from higher states goes into zero
        # density
        theta = 4 * Ne * m
        f = theta * 1/yGrid[1:]
        # probabilities
        preP = f * weights[1:]
        # get them in order
        p = Stationary.normalizeAndAddLeft (preP, weights)
        # give it away now
        return p

    # stationary with selection
    # for now just evaluate f, maybe in future we can do integrals
    @staticmethod
    def stationarySelectionPRF (Ne, m, s1, s2, yGrid, weights):
        assert (numpy.isclose (yGrid[0], 0))
        assert (numpy.isclose (yGrid[-1], 1))
        # all the mass missing from higher states goes into zero
        # hmm, here we have to remove both sides
        polyGrid = yGrid[1:-1]
        polyWeights = weights[1:-1]
        # first the inverse of the speed measure
        theta = 4 * Ne * m
        factor = theta * numpy.exp(meanScaledFitness(polyGrid, s1, s2, Ne)) / (polyGrid*(1-polyGrid))
        # and then the ratio of the speed integrals
        ratio = intScalePRF (polyGrid, s1, s2, Ne)
        ratio /= intScalePRF (0, s1, s2, Ne)
        # put it together
        preF = factor * ratio
        preP = preF * polyWeights
        # so now, how to divvy it up at the boundaries
        # for now, f at right boundary is just same value as closest one (or maybe zero?)
        # and left boundary is sum of all the other
        # IMPORTANT not use polyWeights here
        rightP = preF[-1] * weights[-1]
        preP = numpy.concatenate ((preP, [rightP]))
        # that's maybe it for now
        # get them in order
        p = Stationary.normalizeAndAddLeft (preP, weights)
        # give it away now
        return p

    # beta distribution
    # for now take the density, but we could also try getting at integrals with incomplete beta function
    @staticmethod
    def stationaryNeutralBeta (Ne, mAlpha, mBeta, yGrid, weights):
        # make sure about boundaries
        assert (numpy.isclose (yGrid[0], 0))
        assert (numpy.isclose (yGrid[-1], 1))
        assert (numpy.isclose (weights[0], weights[-1]))
        # popgen params
        alpha = 4 * Ne * mAlpha
        beta = 4 * Ne * mBeta
        # need to get boundary weights
        mean = alpha/(alpha+beta)
        # get beta for non-boundary
        preF = scipy.stats.beta.pdf (yGrid[1:-1], alpha, beta)
        # probabilities
        preP = preF * weights[1:-1]
        # now get the boundaries in order
        interiorMean = numpy.sum (yGrid[1:-1] * preP)
        # difference has to come from boundaries
        diffMean = mean - interiorMean
        # this comes from the linear system
        pOne = diffMean
        pZero = 1 - numpy.sum(preP) - pOne
        # get them in order
        p = numpy.concatenate (([pZero], preP, [pOne]))
        # give it away now
        return p    

    # beta with selection
    # bit more tricky, cause we don't know the mean and integrals are a bit complicated
    # maybe we cannot do much better than extrapolating from close to boundary
    @staticmethod
    def stationarySelectionBeta (Ne, mAlpha, mBeta, s1, s2, yGrid, weights):
        # make sure about boundaries
        assert (numpy.isclose (yGrid[0], 0))
        assert (numpy.isclose (yGrid[-1], 1))
        assert (numpy.isclose (weights[0], weights[-1]))
        # popgen params
        alpha = 4 * Ne * mAlpha
        beta = 4 * Ne * mBeta
        # get modified beta for non-boundary 
        # should be e^{meanFitness}
        daFunc = lambda y : numpy.exp(diffusion_core.meanScaledFitness(y, s1, s2, Ne))
        preF = numpy.power(yGrid[1:-1], alpha-1) * numpy.power(1-yGrid[1:-1], beta-1) * daFunc (yGrid[1:-1])
        # get the normalizing constant
        normConst = scipy.integrate.quad (daFunc, 0, 1, weight='alg', wvar=(alpha-1,beta-1))
        preF = preF / normConst[0]
        # get the mean using quadrature
        preMean = scipy.integrate.quad (lambda x : daFunc(x) * x, 0, 1, weight='alg', wvar=(alpha-1,beta-1))
        daMean = preMean[0] / normConst[0]
        # probabilities
        preP = preF * weights[1:-1]
        # get the mean so far
        interiorMean = numpy.sum (yGrid[1:-1] * preP)
        # difference has to come from boundaries
        diffMean = daMean - interiorMean
        # this comes from the linear system
        pOne = diffMean
        pZero = 1 - numpy.sum(preP) - pOne
        # get them in order
        assert (pZero >= -EPSILON)
        assert (pOne >= -EPSILON)
        pZero = numpy.clip (pZero, 0 ,1)
        pOne = numpy.clip (pOne, 0 ,1)
        p = numpy.concatenate (([pZero], preP, [pOne]))
        # give it away now
        return p
