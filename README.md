# CosmoCov (Cibola Edition) ... for more compplex covariances, please go to the CosmoCov_Saguaro_Edition repo 
Xiao Fang, Elisabeth Krause, Tim Eifler

Configuration Space Covariances for Projected Galaxy 2-Point Statistics, built on **CosmoLike**. We provide a flat sky covariance module, computed with the 2D-FFTLog algorithm, and a curved sky covariance module. We also provide a response covariance module written by Alex Barreira and Elisabeth Krause.

The non-Gaussian (NG) covariances (including connected NG and super-sample covariance) are by default evaluated using halo model.

[-> Papers to Cite](#papers-to-cite)

For more details, see [CosmoCov_Notes.pdf](CosmoCov_Notes.pdf).

## Quick Guide

CosmoCov require a recent <span>gcc</span> compiler,
as well as the <span>gsl</span> and <span>FFTW</span> libraries. 
To get started with computing covariances, follow the steps below:

1.  Clone this repository to the directory you like;

2.  Navigate to the <span>covs</span> directory and run command:  
    ```shell
    $ make covs
    ```

3.  Several executables are created: `cov_flat_fft` is
    for the flat sky covariances (using 2D-FFTLog algorithm);
    `cov` is for the curved sky covariances (using Limber
    approximation for galaxy-galaxy lensing and cosmic shear angular
    spectra, but non-Limber for galaxy clustering angular spectra);

4.  As an exmaple, run command:  
    ```shell
    $ ./cov 1 ini_files/cov_test_g.ini
    ```
    to compute the 1st block of the curved sky 3x2pt Gaussian
    covariance of the test example, specified by the ini file
    `ini_files/cov_test_g.ini`. The result of this
    block is output to the directory:
    `output/out_cov_test`.

5.  There are 66 blocks in total for the example run above. One can
    compute all of them by running command:  
    ```shell
    $ for i in {1..66}; do ./cov $i ini_files/cov_test_g.ini; done
    ```
    or
    ```shell
    $ echo {1..66} | xargs -n 1 -I{} ./cov {} ini_files/cov_test_g.ini
    ```

    **Warning**: non-Gaussian covariances can take hours to compute for
    each block, so one may consider computing different blocks in
    parallel.

6.  After all the blocks are computed, you can make a plot of the precision
    matrix by first combining all the blocks and then running the
    provided plotting script <span>plot.py</span>:  
    ```shell
    $ f="cov_test"; cat output/out_cov_test/t* > $f; python plot.py $f
    ```
    The combined covariance file is specified by variable
    `f`, and the plot (in <span>.pdf</span>) will be saved in
    the same directory.

### Input
The ini files contain all the settings, including

  - `Omega_m, Omega_v, sigma_8, n_spec, w0, wa, omb, h0`: the cosmological parameters,

  - `area`: the survey area (in square degrees),

  - `c_footprint_file`: (optional) a footprint file containing the mask power spectrum, which is read-in in `C_survey_window` in `cosmolike_core/theory/covariances_3D.c`; note that the normalization of the power spectrum is automatically adjusted,

  - `clustering_REDSHIFT_FILE, shear_REDSHIFT_FILE, lens_tomobins, source_tomobins, lens_n_gal, source_n_gal`: the lens and source galaxy samples (file paths, the numbers of tomographic bins, the number densities in each bin). The redshift file has (number of tomo bin + 1) columns, in which the 1st column is the z_min of each z bin.

  - `sigma_e`: the total shape noise of the weak lensing measurement,

  - `lens_tomogbias`: the linear galaxy bias parameter of each lens galaxy bin,
  
  - `lens_tomo_bmag`: the magnification bias parameter of each lens galaxy bin (with `b_mag` described in Section 5.1.3 of [Fang et al. (arXiv:1911.11947)](https://arxiv.org/abs/1911.11947)),

  - `IA`: 0 or 1, the switch of running the intrinsic alignment NLA model,

  - `A_ia`, `eta_ia`: the parameters of the NLA model (see Eq. 4.9 of [Fang et al. (arXiv:1911.11947)](https://arxiv.org/abs/1911.11947), but with `A_ia` represented by `a_IA` in the equation),

  - `tmin, tmax, ntheta`: min and max of the angles in arcmins, and the number of
    logarithmically spaced bins, specifying the binning of the angular correlation functions,

  - `ng`: 0 or 1, the switch of running the non-Gaussian covariances,

  - `cng`: 0 or 1, the switch of including the connected non-Gaussian contribution in
    the non-Gaussian computation,

  - `outdir, filename, ss, ls, ll` : the path and filename prefix of the output, the options of computing
    blocks of the covariance involving the shape-shape (ss),
    position-shape (ls), position-position (ll) angular correlation
    functions. Computing 3x2pt covariance means setting all of
    them as `true`,

  - `linear_binning`: 0 (default) or 1, the optional switch of computing covariances in linear angular binning, currently only supported in curved sky covariance routine,

  - `full_tomo`: 0 (default) or 1, the optional switch of including the full cross tomographic bin clustering correlations, rather than the default auto clustering correlations, currently only supported in curved sky covariance routine.

### Output
The covariances will be output as separate blocks in `output/out_cov_.../`, with each block representing the covariance matrix of two 2-point functions. The header of each file contains a list of papers to be cited based on the module used.

The ordering of the data vector this covariance corresponds to is also output in various `order_...` files. The columns are
  
  - column 0: index i;
    
  - column 1: bin-averaged angular scale (in radians);
    
  - column 2: the type of the 2-point function `w, gammat, xi+, xi-`;
    
  - column 3, 4: the {s: source, l: lens} galaxy tomographic bin index 1 and 2.
  
The columns of each covariance block (in `output/out_cov_.../`) are

  - column 0, 1: matrix indices of the element in the full covariance matrix;

  - column 2, 3: the corresponding bin-averaged angular separations (in
    radians) of the element;

  - column 4, 5, 6, 7: the tomographic bins involved;

  - column 8, 9: the Gaussian part and the non-Gaussian part of the
    element. The total value is the sum of the two.

## Papers to Cite

Please cite the following papers if you use covariances in your
research:

1.  [E. Krause, T. Eifler; *CosmoLike - Cosmological Likelihood Analyses
    for Photometric Galaxy Surveys*;
    arXiv:1601.05779](https://ui.adsabs.harvard.edu/abs/2017MNRAS.470.2100K/exportcitation)

2.  [X. Fang, T. Eifler, E. Krause; *2D-FFTLog: Efficient computation of
    real space covariance matrices for galaxy clustering and weak
    lensing*; arXiv:2004.04833](https://arxiv.org/abs/2004.04833)

In addition, <span class="underline">if you use the the non-Limber
galaxy clustering power spectra in the Gaussian covariance</span>
included in the curved sky covariance module, please also cite:

  - [X. Fang, E. Krause, T. Eifler, N. MacCrann; *Beyond Limber:
    Efficient computation of angular power spectra for galaxy clustering
    and weak lensing*;
    arXiv:1911.11947](https://ui.adsabs.harvard.edu/abs/2019arXiv191111947F/exportcitation).

<span class="underline">If you include an approximate shape/shot noise
correction for survey geometry</span>, which requires a power spectrum
of the survey mask with sufficiently high resolution (see [\[Quick Guide\]](#quick-guide)), please also cite:

  - [M. Troxel, E. Krause, et al. (DES Collaboration); *Survey geometry
    and the internal consistency of recent cosmic shear measurements*;
    arXiv:1804.10663](https://ui.adsabs.harvard.edu/abs/2018MNRAS.479.4998T/exportcitation);

<span class="underline">If you calculate weak lensing covariances using
the response model</span>, including the non-Limber super-sample
covariance (available for weak lensing only), please also cite:

  - [C. Wagner, F. Schmidt, C.-T. Ching, E. Komatsu; *The angle-averaged
    squeezed limit of nonlinear matter N-point functions*;
    arXiv:1503.03487](https://ui.adsabs.harvard.edu/abs/2015JCAP...08..042W/exportcitation)

  - [A. Barreira, F. Schmidt; *Responses in Large-Scale Structure*;
    arXiv:1703.09212](https://ui.adsabs.harvard.edu/abs/2017JCAP...06..053B/exportcitation);

  - [A. Barreira, F. Schmidt; *Response Approach to the Matter Power
    Spectrum Covariance*;
    arXiv:1705.01092](https://ui.adsabs.harvard.edu/abs/2017JCAP...11..051B/exportcitation);

  - [A. Barreira, E. Krause, F. Schmidt; *Complete super-sample lensing
    covariance in the response approach*;
    arXiv:1711.07467](https://ui.adsabs.harvard.edu/abs/2018JCAP...06..015B/exportcitation);
    
  - [A.S. Schmidt, S.D.M. White, F. Schmidt, J. Stücker; *Cosmological
    N-body simulations with a large-scale tidal field*;
    arXiv:1803.03274](https://ui.adsabs.harvard.edu/abs/2018MNRAS.479..162S/exportcitation)

_A list of papers to be cited will be autogenerated based on the modules used, and printed in the header of each covariance file._

If you use the given DES Y3 ini files, they will assume public DES Y1
galaxy redshift distributions from
<http://desdr-server.ncsa.illinois.edu/despublic/y1a1_files/redshift_bins/y1_redshift_distributions_v1.fits>
