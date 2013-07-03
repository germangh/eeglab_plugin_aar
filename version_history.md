Version history
===

## Version 1.3

### Release 05-05-2009
  - A minor bug in function fd.m has been corrected.


### Release 06-04-2009, version 1.3
  - Added the possibility of adding labels to chunks of EEG data via GUI.
    This option works only for continuous datasets (not for epoched datasets)
    In future releases epoched datasets will also be supported. This feature
    is not yet described in the current version of the toolbox documentation.
  - Few minor bugs corrected.

### Release 31-10-2008
  - IMPORTANT: The order of the input parameters of functions pop_autobssemg
    pop_autobsseog has changed.
  - Corrected few bugs related with running pop_autobssemg and pop_autobsseog
    script mode (i.e. from the command line).

### Release 09-12-2007, minor changes
  - SOBI is again the default algorithm for BSS-based EOG correction
  - New documentation with several correction examples.

### Release 04-12-2007, few major bugs corrected
  - A few minor and major bugs were still left in the BSS interface
    functions and have now been corrected.
  - The interface functions use now pinv instead of inv to avoid
    numerically unstable results when the data covariance matrix is close
    to singular.
  - The default window length in the automatic EMG correction method is
    now twice as long as it was before.

### Release 03-12-2007, major bugs corrected
  - Several major bugs have been corrected in the functions that act as
    interface to the BSS algorithms.
  - Parameters passed to the BSS algorithms were ignored in previous
    releases due to a bug in function autobss.m. This has been corrected.

### Release 29-11-2007
  - Minor bugs corrected.

### Release 28-11-2007, major bugs corrected
  * Thanks to Dr. Ranjit Thuraisingham for noting several of the important
    bugs that have been corrected in this new release.
  - A major bug in function emg_psd.m has been corrected. This bug was
    causing the emg_psd criterion to incorrectly detect the EMG-related
    components.
  - Minor update of function pop_autobssemg.m to adapt to the new version
    of function emg_psd.m
  - The criterion emg_psd.m now accepts a new parameter (opt.estimator)
    that allows the user to select the spectral estimator to be used to
    compute the EEG average power and the EMG average power of each esti-
    mated component.
  - A compatibility issue related to function emg_psd.m has been solved.
    Previouly to this change, the criterion emg_psd was producing
    different results when using MATLAB's Signal Processing Toolbox 6.6
    and previous when using older versions of the Signal Processing Toolbox.
    Now it produces "almost" the same results. We are not yet sure of the
    reason for the small remaining differences.
  - Several other minor bugs corrected.

### Release 07-11-2007, minor update
  - pop_autobssemg now automatically passes the sampling rate to emg_psd
  - first draft of the toolbox documentation

### Release 31-10-2007, minor update
  - Minor bugs corrected
  - EFICA v1.9 was updated with EFICA v2.0 (implemented by Z. Koldovsk˝)

### Release 29-10-2007, new features
  - Few minor bugs corrected
  - Slightly improved documentation
  - A sample EEG dataset has been included in the release.
  - Added several new algorithms for BSS: iWASOBI, EFICA, COMBI, FCOMBI.

## Version 1.2

### Release 01-07-2006, new features
  - Automatic EOG correction using Least Mean Squares (LMS) adaptive filtering.
    This is implemented in lms_regression.m.
  - Automatic EOG correction using the conventional Recursive Least
    Squares (RLS) adaptive algorithm. This is implemented in crls_regression.m.
  - Automatic EOG correction using an stable version of the RLS algorithm.
    This is implemented in scrls_regression.m.
  - Automatic EOG correction using two H-infinity regression methods. These
    are implemented in hinftv_regression.m and hinfew_regression.m.
  - Automatic EOG correction using Principal Component Analysis (PCA).
  - Some bugs corrected.


## Version 1.1

###  Release 10-04-2006, new features
  - Automatic correction of EMG artifacts using canonical correlation analysis.
    This is implemented in functions bsscca.m and emg_psd.m
  - Function BSS has been modified to separate the actual BSS algorithm
    from the components selection criteria.
  - Interfaces to BSSCCA, FastICA, JADE and RUNICA (implementation of Infomax
    included in the EEGLAB toolbox)  have been included. Please see the
    comments of fastica.m, jader.m and runica.m for copyright information.
  - Ill-conditioning of the data covariance matrix is now handled properly.
  - Some minor bugs corrected.

## Version 1.0

### Release 27-02-2006
  - Automatic correction of EOG artifacts using Blind Source Separation (BSS)
  - Only an interface to SOBI is provided

