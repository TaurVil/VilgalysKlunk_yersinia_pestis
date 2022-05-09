# `diplo_locus`: Python libraries for simulating and computing log-likelihoods for time-stratified genetic samples based on diploid Wright-Fisher diffusion
> still in beta

# Table of Contents

* [`diffusion_core` module](#diffusion_core)
* [`likelihood` module](#likelihood)
  * [class `SelHmm`](#likelihood.SelHmm)
    * [Usage](#SelHmm.init) 
    * [Function](#SelHmm.compute)
* [`simulate` module](#simulate)


<a id="diffusion_core"></a>

## `diffusion_core` module

Python library for core functions to compute HMM for WF diploid selection

<a id="likelihood"></a>

## `likelihood` module

Python library for HMM functions to compute per-site likelihood based on `diffusion_core` module.

<a id="likelihood.SelHmm"></a>

### *class* `SelHmm`

A `SelHmm` object construct a diploid Wright-Fisher diffusion HMM, configures its population genetic parameters, and
compute log-likelihoods of samples observed under the presumed parameters.

<a id="SelHmm.init"></a>
### Usage
```python
>>> SelHmmObject = likelihood.SelHmm(Ne, s1, s2, mAlpha, mBeta, initCond, initFreq=None,
                                      initMAlpha=None, initMBeta=None, initS1=None, initS2=None,
                                      sampleSizesSet=None, numStates=1001, deltaT=1, 
                                      emissionType="integer", transitionType="constant", 
                                      selectionChangeTimes=None)
```

<details>
<summary>Click to see parameter details</summary>
 
#### Parameters

- **`Ne`** `int or float
`

  Effective diploid population size (*i.e.*, for each locus, 2Ne alleles exist in total).

- **`s1`** `float or array_like
`

   Selection coefficients of the heterozygote.

- **`s2`** `float or array_like
`

   Selection coefficients of the homozygote.

- **`mAlpha`** `float
`

   Per-site per-generation forward mutation rate.

- **`mBeta`** `float
`

   Per-site per-generation backward mutation rate.

- **`initCond`** `{"uniform"` `"initFreq"` `"statBeta"` `"statBetaSel"}
`

   Specify the initial condition for the HMM at generation zero. Depends on the type of initial condition, other parameters need to be provided:
    * `"uniform"` :
           HMM starts with a uniform distribution on [0,1].
    * `"initFreq"` :
           Starts with a given allele frequency. Must also specify the frequency with `initFreq`.
    * `"statBeta"` :
           Starts with the stationary Beta distribution under mutation rates `mAlpha` and `mBeta`.
    * `"statBetaSel"` :
           Starts with the stationary Beta distribution under mutation rate `initAlpha`, `initBeta`, and selection coefficient `initS1`, `initS2`.

- **`emissionType`** `{"integer"` `"fractional"}` `default="integer"
`

   Type of emission data, *i.e.* observed sample numbers, to consider. Choose `"integer"` for integer allele counts, `"fractional"` for when such count has been adjusted or a fractional estimate.

- **`transitionType`** `{"constant"` `"piecewise"` `"continuous"}` `default="constant"
`

   Type of transitions to consider in the HMM:
    - `"constant"` :
           Selection coefficients stay constant throughout the entire duration considered.
        `s1` and `s2` must be constants (not `array_like`) for this option.
    - `"piecewise"` :
           The entire duration can be considered as several pieces in tandem, where each piece has a different pair of selection coefficients. Must also specify `selectionChangeTimes` for this option.

#### Other Parameters

- **`deltaT`** `int or float` `optional` `default=1
`
  Unit increment of time (in generations).

- **`numStates`** `int or float` `optional` `default=1001
`
   Number of discretized states with which to discretize the allele frequency space [0,1].

- **`initFreq`** `float` `optional
`
   Required when ``initCond="initFreq"``. Must be between 0 and 1.

- **`initAlpha`** `float` `optional` `default=mAlpha
`
   Per-site per-generation forward mutation rate underlying the initial Beta distribution.

- **`initBeta`** `float` `optional` `default=mBeta
`
   Per-site per-generation backward mutation rate underlying the initial Beta distribution.

- **`sampleSizesSet`** `set of array_like
`
   Set of list, numpy array, or tuple objects of the same length as the number of sampling time points that summarizes all possible sample sizes

- **`selectionChangeTimes`** `int or array_like` `optional
`
   Set the generation time when selection coefficients change. Must match the length of `allS1` and `allS2`.

 </details>

<a id="SelHmm.compute"></a>
### Function

#### `likelihood.SelHmm.computeLogLikelihood()`

```python
>>> LLs = SelHmmObject.computeLogLikelihood(times, samples, sampleSizes)
```

After specifying the underlying parameters of a `SelHmm` object, this function compute the log-likelihoods of observing the given samples at the sampling times `times`.
Note that this model assumes all loci are independent and bi-allelic.

<details>
<summary>Click to see parameter details</summary>

#### Parameters

- **`times`** `array_like
`

The generation times when the samples were taken. Must start with zero and ascend from fast to present.

- **`samplesSizes`** `array_like
`

   An N by K matrix of the total numbers of alleles observed, *i.e.* sample sizes.
    N --> number of loci ;  K --> number of sampling times.
    ``sampleSizes[i,j]`` records the sample size of locus ``i`` at time point ``j``.

- **`samples`** `array_like
`

   An N by K matrix of the numbers of alleles. N --> number of loci ;  K --> number of sampling times.
    ``samples[i,j]`` records the number of a particular allele on locus ``i`` at time point ``j``.

#### Returns

   numpy.array of log-likelihood for each locus.

</details>


<a id="simulate"></a>

## `simulate` module

Python library for functions to perform simulations under a Wright-Fisher diffusion and taking samples at given times.

#### `simulate.simulateSamples()`

```python
>>> samples, trajectories = simulate.simulateSamples(Ne, allS1, allS2, mAlpha, mBeta, times, sampleSizes,
                                                     initCond=None, initFreq=None, numReplicates=1,deltaT=0.05,
                                                     condInitSeg=True, initGridResolution=None,
                                                     initProbs=None, initValues=None, initMAlpha=None, initMBeta=None,
                                                     selectionChangeTimes=None)
```

The `simulate.computeSamples()` function takes relevant population genetic parameters and simulates time-stratified samples of independent loci and reports the allele frequency trajectories of simulated loci.

<details>
<summary>Click to see parameter details</summary>
 
#### Parameters

- **`Ne`** `int` `float
`

  Effective diploid population size (*i.e.*, for each locus, 2Ne alleles exist in total).

- **`allS1`** `int` `float` `numpy numbers` `or array_like
`

   Selection coefficient of the heterozygote. When simulating piecewise time-varying selection, the length of `list` or `numpy.ndarray` should match that of `selectionChangeTimes` so that ``len(allS1) = len(selectionChangeTimes) + 1``.

- **`allS2`** `int` `float` `numpy numbers` `list` `or array_like
`

   Selection coefficient of the homozygote. Requirements are the same as `allS1`.

- **`mAlpha`** `float
`

   Per-site per-generation forward mutation rate.

- **`mBeta`** `float
`

   Per-site per-generation backward mutation rate.

- **`times`** `array_like` `int
`

   Numbers in forward generation times when samples were taken. Must start with zero and satisfy ``len(times) = len(sampleSizes) + 1``.

- **`samplesSizes`** `array_like
`

   A list or array of length K recording the total numbers of alleles observed, *i.e.* sample sizes, at each samping time point.
    ``sampleSizes[i]`` records the sample size at time point ``i``.
    
- **`initCond`** `{"initFreq"` `"contBeta"` `"discBeta"` `"discBetaSel"` `"choice"}
`

   Indicate how the initial condition will be decided for each replicate. Based on the selection and input from their related arguments, a probability density will be determined, by which an initial frequency will be sampled for each simulated replicate at t = 0. Below are details for these options:
- `"initFreq"` :
       Simulation starts with a fixed given frequency `initFreq`.
- `"contBeta"` :
       Initial frequency will be drawn from a stationary Beta distribution. Use `initAlpha` and `initBeta` to specify its parameters. By default, they are set to equate `mAlpha` and `mBeta`, respectively.
- `"discBeta"` :
       Initial frequency will be drawn from a discretized Beta distribution. Require `initGridResolution` to specify the number of bins to discretize (0,1). Use `initAlpha` and `initBeta` to specify its parameters.
- `"discBetaSel"` :
       Initial frequency will be drawn from the stationary Beta distribution under `initAlpha`, `initBeta`, and selection `initS1`, `initS2`
- `"choice"` :
       Draw initial frequency from the custom defined `initProbs`.

- **`numReplicates`** `int` `optional` `default=1
`

   Number of independent replicate loci to simulate.


#### Other Parameters

- **`initFreq`** `float` `optional
`

  Required when `initCond="initFreq"`. Must be between 0 and 1.

- **`initAlpha`** `float` `optional` `default=`mAlpha`
`

   Per-site per-generation forward mutation rate underlying the initial Beta distribution.

- **`initBeta`** `float` `optional` `default=`mBeta`
`

   Per-site per-generation backward mutation rate underlying the initial Beta distribution.

- **`deltaT`** `float` `default=0.05
`

   Unit increment of time (in generations).

- **`condInitSeq`** `bool` `default=True
`

   Set whether to condition simulations on the initial samples being segregating.

- **`initGridResolution`** `int or float
`

   Required when initial condition is set to be "discBeta" or "discBetaSel". Number of bins to discretize (0,1).

- **`selectionChangeTimes`** `int or array_like` `optional
`

   Set the generation time when selection coefficients change. Must match the length of `allS1` and `allS2`.

- **`initProbs`** `numpy.ndarray` `optional
`

   Required when initial condition is set to be `"choice"`. This is an array of allele frequencies for which `initValue` specifies their corresponding probability densities.

- **`initValues`** `numpy.ndarray` `optional
`

   Required when initial condition is set to be `"choice"`. Provides probability densities for the allele frequencies in `initProb`.

#### Returns


- **`samples`** `numpy.ndarray
`

   A `numReplicates` x K matrix of simulated samples. K --> number of sampling times.
   ``samples[i,j]`` records the number of derived alleles observed on replicate ``i`` at sampling time point ``j``.


- **`wfDiffReplicates`** `numpy.ndarray
`

    A `numReplicates` x steps matrix of simulated allele frequency trajectories, with steps = #(generations spanned) / `deltaT`.

 </details>


