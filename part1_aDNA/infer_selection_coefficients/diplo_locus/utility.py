import numpy




# get a linear grid
# include 0 and 1
def getLinearGrid (numStates):
    return numpy.linspace (0,1,numStates)


# get linear boundaries
# start from 0, got to 1
def getLinearBounds (yGrid):
    bounds = numpy.concatenate (([0], yGrid[:-1] + numpy.diff(yGrid)/2, [1]))
    assert (len(bounds) == len(yGrid)+1)
    return bounds


# get weights
def getWeights (bounds):
    weights = numpy.diff(bounds)
    assert (len(bounds) == len(weights)+1)
    return weights

