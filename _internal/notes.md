# Notes

## 2025-05-19

Finalized the calibration tests

## 2025-05-20

Today I am going to test the stochastic ALP field where $\tau_a \ll T_2$.

With the data, $\sqrt{\tau_a T_2}$ and $\sqrt{\tau_a \Delta T}$ are proven to be true for estimating tipping angle.

## 2025-05-26

#### Quantitive analysis of the M_transvers over time

$M_t = M_0 \, \sin\theta \approx  M_0 \,\theta$ assuming $\theta\ll1$. We have $M_0=1$, and

$\theta=\Omega \Delta t = \gamma B_a \Delta t$.

The keypoint here is the strength of the magnetic field $B_a$.


## 2025-06-02

### Study the effect of the frequency detuning on the tipping angle

Using $\tau_a \ll T_2$.

I used 40 $\nu_a$'s and wanted to see the tipping angle evolution over time. It is like scanning over the frequencies.

Findings:

1. When on resonance ($\nu_a\approx\nu_\mathrm{Larmor}$), the tipping angle growes up first ($t\ll T_2 \text{ or } \tau_a$), then saturates after $t\gg T_2 \text{ or } \tau_a$.

2. When off resonance ($|\nu_a-\nu_\mathrm{Larmor}|\gg \Delta\nu_a$), the tipping angle growes up first ($t\ll T_2 \text{ or } \tau_a$), then declines by a small amount, and becomes steady after $t\gg T_2 \text{ or } \tau_a$.

The data and figures are in the folder /20250602-tau_a_《_T2.

One of the conclusions: the difference between $\nu_a$ and $\nu_\mathrm{Larmor}$ matters in the magnitude of the tipping angle.

## 2025-06-03

what to do today...

I can use the lineshape of the tipping angle magnitude to quantitvely evaluate the effect of axion field. It should be something like this:

$\theta^2 \propto \lambda(\nu)$

or more precisely:

$\theta^2 \approx (\gamma B_a)^2\, \lambda(\nu)$ nono this is not right. check its unit!

## 2025-08-08

### magnitude test

How to quantitvely evaluate the effect of axion field? Here I do not consider the T2*, but only detuning $\Delta\nu$, $\tau_a$, and $T_2$. Assuming measurement time much longer than any one of the times above.

We expected

$\theta^2 \approx (\gamma B_a T_2)^2\, \lambda(\nu)$ when $T_2 \ll \tau_a$

$\theta^2 \approx^? (\gamma B_a )^2 \tau_a T_2\, \lambda(\nu)$ when $\tau_a\ll T_2\ll T_\mathrm{meas}$.

Let us verify these.

## 2025-09-19

### Brms magnitude test

in [src\tests\20250919-Bamp-Brms-fft-ifft-magnitude-test\check_B_rms.py] we test the relationship between Bamp (input axion Ba field amplitude) and the output Brms (rms of the simulated axion B field).

In principle Brms is equal to Bamp. However, due to numpy.fft and numpy.ifft, the amplitude of simulated Ba field should be Bamp $\times$ some N, where N is simulation rate or array length.

Test 0:
------------------------
(input) Bamp = 1e-10

(output) np.mean(B_rms_from_simu_arr) : float64(3.0370544560166723e-12)

np.std(B_rms_from_simu_arr) : float64(3.161938491415789e-13)

simuRate * duration = 500 * 20 = 10000

--------------------
np.mean(B_rms_from_simu_arr) : float64(2.7800346877931807e-12)

np.std(B_rms_from_simu_arr) : float64(8.35364629679566e-13)

simuRate * duration = 500 * 2 = 1000

---------------------------------------
np.mean(B_rms_from_simu_arr) : float64(2.733154499609955e-12)

np.std(B_rms_from_simu_arr) : float64(1.2322121825485608e-12)

simuRate * duration = 500 * 2 = 1000

---------------------

np.mean(B_rms_from_simu_arr) : float64(1.0460985972196068e-11)

np.std(B_rms_from_simu_arr) : float64(3.5039550658350867e-12)

simuRate * duration = 50 * 2 = 100

--------------

line: ax_FFT = Bamp * ax_lineshape * rvs_phase * **simuRate * np.sqrt(duration)**

fixed the problem.

Note that Brms of Bx or By is = Bamp / np.sqrt(2)

## 2025-10-21

### Plan

I am going to creatre a few classes for better clarifying the simulation conditions, including:

1. SQUID
2. Pickup
3. Sample
4. AxionWind
5. MagField
6. Simulation
7. Experiment: LF or HF -> give fluxpower, exclusions, etc.
8.


## 2025-11-13

### wiggles in the Larmor precession 

I found that, in the RF-pulse.py, where a pulsed-NMR decay signal was simulated, I can find wiggles in the signal after a few Tdelta (T2 and T1 are much longer). I tried to elminate such non-physical wiggles. (in principle we should expect a flat curve of Mx and My after a few Tdelta. )

In the end I found that 
    1. nFWHM
    2. numPt 
    3. the distribution of Larmor frequencies
are important for eliminating the wiggles. 
Before, we used:
    1. nFWHM=5
    2. numPt < int(
        11
        + duration.value_in("s")
        * 2
        * magnet_det.nFWHM
        * magnet_det.FWHM_T
        * sample.gamma.value_in("Hz/T") * 1
    )
    3. the distribution of Larmor frequencies: uniform across [-nFWHM *FWHM, nFWHM*FWHM]

Later we used:
    1. nFWHM=20
    2. numPt = int(
        11
        + duration.value_in("s")
        * 2
        * magnet_det.nFWHM
        * magnet_det.FWHM_T
        * sample.gamma.value_in("Hz/T") * 1
    )
    3. the distribution of Larmor frequencies: square of uniform across [-nFWHM * FWHM, nFWHM * FWHM]: self.B_vals_T = (
                self.nFWHM * np.sign(u) * np.abs(u) ** 2
            ) * self.FWHM_T + self.B0.value_in(
                "T"
            )

Then the problem of wiggles is gone. This indicates that:
    1. The frequency band needs to be large enough. Those high-frequency signals are also important. 
    2. Enough data points! Though I believe 7/12 of the current number of data points may be sufficient. 
    3. Give more data points to low-frequency band. 


## 2025-11-23

### simulation optimization 

---

`@nb.jit` is less efficient than using an explicit signature, such as:

```python
@nb.jit(
    [
        "void(float64[:,:],...)"
    ],
    nopython=True,
)
```

for `generateTrajectoryLoop()`.

---

### Comparison of simulation approaches

The comparison script is located at
`src/tests/optimization-test/RF-CW-compare-speed.py`.

When the bias field is inhomogeneous (requiring many spin packets to be simulated), the performance differs significantly between the implementations.

* `generateTrajectory` and `generateTrajectoryLoop`
  Outer loop over spin packets in Python, inner loop accelerated with
  `@nb.jit(nopython=True)`

  Time consumption: `0.214855 s`

* `generateTrajectory_doubleLoop` and `generateTrajectoryLoop_doubleLoop`
  Outer loop over spin packets and inner loop both in `generateTrajectoryLoop_doubleLoop`, accelerated with
  `@nb.jit(nopython=True)`

  Time consumption: `0.071009 s` 

* `generateTrajectory_vectorized` and `generateTrajectoryLoop_vectorized`
  vectorize over the outer loop (different spin paets) and keep the inner loop, accelerated with
  `@nb.jit(nopython=True)`

  Time consumption: `0.06 s` 


another test:

magnet_det.numPt = 2811

generateTrajectory time consumption = 0.736619 s

generateTrajectory_doubleLoop time consumption = 0.210296 s

generateTrajectory_vectorized time consumption = 0.172490 s

vectorization method wins. 
---

### test of time consumption on dell XPS

* 1 GHz
    individual step time consumption = **1.731e-08 s**

    individual step time consumption = **1.473e-08 s**

    individual step time consumption = **1.692e-08 s**

* Axion Compton frequency = **1.0 megahertz**

    simulation duration = 5.269622e+00 (s).

    numPt = 673.2005774553354

    magnet_det.numPt = 673

    simuRate = 207.079727562906 hertz

    Number of simulation steps = 1091.

    individual step time consumption = **1.833e-08 s**

    individual step time consumption = **1.364e-08 s**

    individual step time consumption = **1.787e-08 s**

* Axion Compton frequency = **1.0 kilohertz**

    simulation duration = 5.269622e+03 (s).

    numPt = 673.2005774553355

    magnet_det.numPt = 673

    simuRate = 50.1570797275629 hertz

    Number of simulation steps = 264308.

    individual step time consumption = **1.282e-08 s**

    individual step time consumption = **1.280e-08 s**

    individual step time consumption = **1.300e-08 s**


### Compare ifft methods in generating axion time series in terms of efficiency

Array size: 200000000

Non-zero band width: 2000

NumPy ifft               : 12.357230 s

Manual pad + NumPy ifft  : 13.801974 s

SciPy ifft               : 12.042085 s

FFTW ifft                : 18.020005 s

conclusion: numpy and scipy are sufficient

### Cutoff / upper limit for nu in axion_lineshape() to increase robustness when nu is large

`frequencies = nu_a + np.linspace(0, nu_a * 2e4 / 1e6, 10_000, endpoint=True)`

On dell XPS: 

before: time consumption = 0.001209 s. 

after cutting off at 100 axion linewidth (1e-4 / 0.02 = 0.5 % of previous array length): time consumption = 0.000341 s. 

On Ryzen PC : time consumption = 0.000199 s. 

### data analysis

data generated by `src\tests\20251123-signal-at-frequencies\methanol.py`. 

analyzed by `src\tests\20251123-signal-at-frequencies\data-analysis-and-plot.py`. 

For a weak drive (B_rms = 1e-15 T), reasonable measurement time (>= 10 * shortest coherence time), and thermally polarized methanol sample, the signal in time domain grows up at first, and goes to a platform later. This is true for 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9 Hz. 

### Compare estimation and simulation

`src\tests\20251123-signal-at-frequencies\compare-estimation-simulation.py`

tau = 1 / (1/tau_a +1/T2star) is a good estimation for simulation results > 1 kHz, but it is not very good for 1 kHz. Could it be simulation problem? 

## 2025-11-24

### test of differeny ways of operating on 3D arrays in C++

`src\tests\20251124-C++\test_speed.py`

add_3d_arrays_3loops  runtime = 0.0019998550415039062 s: relative run time = 1.8218940052128585

add_3d_arrays_parallel  runtime = 0.003582000732421875 s: relative run time = 3.263249348392702

add_3d_arrays_flattern_parallel  runtime = 0.0009975433349609375 s: **relative run time = 0.9087749782797567**

add_3d_arrays_flattern_SIMD  runtime = 0.0020017623901367188 s: relative run time = 1.8236316246741964

time_Cpy: 0.39560532569885254

Conclusion: use try the best to flattern loops and use parallel

```C++
void _add_3d_arrays_flattern_parallel(const double *A, const double *B, double *C,
                                      std::size_t X, std::size_t Y, std::size_t Z)
{
    std::size_t N = X * Y * Z;

#pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N; ++idx)
    {
        C[idx] = A[idx] + B[idx];
    }
}
```

## 2025-11-27

### Optimize setALPField efficiency

`src\tests\20251127-bloch-C++\test-setField-speed.py`

#### Optimize rng

before optimization:

axion_lineshape time consumption = 3.402e-04 s

rng time consumption = **7.521e-04 s**

ifft and allocation time consumption = 2.582e-03 s

setALP_Field total time consumption = 0.004287 s = 0.000071 (min)

---
axion_lineshape time consumption = 3.060e-04 s

rng time consumption = **7.059e-04 s**

ifft and allocation time consumption = 2.379e-03 s

setALP_Field total time consumption = 0.003944 s = 0.000066 (min)

---

after optimization:

axion_lineshape time consumption = 3.254e-04 s

rng time consumption = **6.227e-04 s**

ifft and allocation time consumption = 2.414e-03 s

setALP_Field total time consumption = 0.004003 s = 0.000067 (min)

---

#### optimize ifft

change from numpy to scipy.fft.ifft and saved ifft run time by ca. 0.2/1.6 = 12.5%. 

----

optimize when we need a lot of fields in one simulation

by adding one more dimension in arrays

before optimization

1 fields

ifft total time consumption = **1.471e-03 s**

---

after optimization

100 fields: 

axion_lineshape time consumption = 3.484e-04 s

rng time consumption = **3.869e-02 s**

ifft total time consumption = **5.255e-02 s**

individual runtimes = [0.002660500002093613, 0.005015600007027388, 0.02148019999731332, 0.002693000016734004, 0.02068309998139739]

array-asignment time consumption = 7.059e-02 s

10_000 fields:

axion_lineshape time consumption = 3.263e-04 s

rng time consumption = **4.270e+00 s**

ifft total time consumption = **6.430e+00 s**

individual runtimes = [0.33017550001386553, 0.7834994000149891, 2.649957799934782, 0.4694043999770656, 2.1972203999757767]

array-asignment time consumption = 1.604e+01 s

[New] setALP_Fields time consumption = **30.307611 s** = 0.505127 (min)

saved about 50 % time

note: when it comes to many dimensions and long arrays, scipy wins over numpy. 

#### array-asignment optimization

saved about 88 % by 

1. viewing but not copying (.real instead of np.real())
1. no np.outer but np.zeros((,,))
1. more dimensions

#### summary of optimizations

1 field:

setALP_Field total time consumption = 0.003444 s = 0.000057 (min)

1000 fields:

[New] setALP_Fields time consumption = 12.520021 s = 0.208667 (min)

12.520021/0.003444/10000 = 0.363531387921022

1-12.520021/0.003444/10000 = 0.636468612078978

**64 %** time saved. 

## 2025-11-30

### optimization of `Simulation.generateTrajectories`

`src\tests\20251130-C-and-trjry\test-trjry-speed.py`

with new method `Simulation.generateTrajectories`: 

1 MHz axion

simulation duration = 1.000000e+01 (s).

magnet_det.numPt = 1267

simuRate = 300.00018121674657 hertz

[setALP_Fields] time consumption estimated = 0.0886 min

[generateTrajectories] time consumption estimated = 0.75 min

Total time consumption estimated = **0.8409638006666668 min**

[setALP_Fields] time consumption = 4.389461 s = **0.073 (min)**

[setALP_Fields] individual step time consumption = 4.389e-04 s

[generateTrajectories] time consumption = 44.771881 s = **0.746 (min)**

[generateTrajectories] individual step time consumption = 1.178e-09 s

-----

with the best old methods `Simulation.generateTrajectory_vectorized`:

Axion Compton frequency = 1000000.0 hertz

simulation duration = 1.000000e+01 (s).

numPt = 1267.637820503248

magnet_det.numPt = 1267

simuRate = 300.00018121674657 hertz

Number of simulation steps = 3000.

Single run time consumption estimated = 0.0008866043666666667 min

Total time consumption estimated = **8.87 min**

-----

Conclusion: even if we do not consider the run time of setting ALP fields in the old method, the new method can still win by saving >90% of time. 

### attempts of new simulations

not really successful because not all information are saved. Tomorrow I will try once more. 


## 2025-12-16

### a complete simulation with proper data saving

Let us check if we save everything rightly:

------

for h5 level-1 group axion_wind:

done. only physical quantities are saved. 

-----------

for sample: done. only physical quantities are saved. 

-----

for magnet_det: done. all quantities are saved. 

------------

for simulation: done. all quantities are saved. 

--------

### camel and snake case

use camel case when it dicates one word / phrase

use snake case if there are multiple words or if there is break in it

-------------

## 2025-12-18

### data saving and loading realized by Simulations

with pickle, now we can save and load instances of Simulations easily

## 2025-12-19

### test s(t) at all frequencies

use 0.01 ppm, nu_a = 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9 Hz

`src\tests\20251219-Simulations\save-Simulations.py`

`src\tests\20251219-Simulations\plot-Simulations.py`

it looks like my estimation in the `src\tests\20251219-Simulations\plot-Simulations.py` is very close to truth. The discrepancy happens when the magnet homogeneity goes to 0.01 ppm. Problems from detuning when lw is narrow? Further investigation is needed. 

### CPMG

added CPMG in the calibrations

I need to refine the pi half and pi pulses in the calibrations

## 2025-12-22

### test 0.01 ppm homogeneity

I have found that, for T2star << tau_a, the signal strictly follows

$s(t) = (1 - e^{-t/T_2^*}) \gamma B T_2^*$. 

When $T_2^*\gg\tau_a$, the final signal amplitude $\approx \gamma B (T_2^*\tau_a)^{1/2}$. And we also found that the signal grows up with $\sim (t\tau_a)^{1/2}$ at first. So we can use

$s(t) = \left(1 - \exp{\left(\dfrac{-t}{\frac{1}{2}T_2^*}\right)}\right)^{1/2} \gamma B (T_2^*\tau_a)^{1/2}$

to describe the signal. 

Besides, for the sack of simplicity in numerical calculations, we can use this formula to approximate the signal under any circumstances

$s(t) = \left(1 - \exp{\left(\dfrac{-t}{x T_2^*}\right)}\right)^{x} \gamma B T_2^* \left(\dfrac{\tau_a}{\tau_a + T_2^*}\right)^{1/2}\,.$

## 2025-12-??

### test detuning?


## 2026-01-23

todo: connect our tests to physical experiments, with real coupling stength for example. 

