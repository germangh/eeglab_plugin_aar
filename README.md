eeglab_plugin_aar
=================

This repository stores the code of the AAR plug-in for [EEGLAB][eeglab] that was
released on 31-10-2008. I am not developing this code anymore but if you report
any major bug I will try to fix it whenever I have time. If you are able to,
I would really appreciate if you fork the repo, fix the bug and make a pull
request.

If you are interested in removing artifacts from EEG time-series, you might find
interesting the [meegpipe][meegpipe] project. You can use [meegpipe][meegpipe]
to do everything you can do with the AAR plug-in and much much more.

[meegpipe]: http://germangh.com/meegpipe


## Overview

The AAR plug-in is a collection of [MATLAB][matlab] scripts that implement
several state-of-art (in 2008 anyways...) methods for automatic correction of
ocular and muscular artifacts in the [EEG]. The capabilities of the toolbox are
briefly described in a fairly outdated [tutorial document][tutorial] in .pdf
format.

[matlab]: http://www.mathworks.nl/products/matlab/

[tutorial]: http://kasku.org/pubs/aardoc07.pdf

The toolbox is implemented as an [EEGLAB][eeglab] plug-in, but can also be used
without EEGLAB if you don't need or want to use EEGLAB's GUI. The current
version implments several fully automatic methods to correct ocular ([EOG][eog])
artifacts, and one automatic method to correct muscle ([EMG][emg]) artifacts.

For suggestions, comments and bug reports, please contact [German
Gomez-Herrero][ggh].

[eog]: http://emedicine.medscape.com/article/1140247-overview#aw2aab6b3

[emg]: http://emedicine.medscape.com/article/1140247-overview#aw2aab6b3

[ggh]: http://germangh.com

## Installation instructions

1. Install [EEGLAB][eeglab] for MATLAB, if you haven't done so already.

[eeglab]: http://sccn.ucsd.edu/eeglab/


## Additional resources

The [BSS][bss]-based EOG correction procedure is based on the following
scientific publication:

* [Gomez-Herrero, G.] et al., _Automatic Removal of Ocular Artifacts in the EEG
  without an EOG Reference Channel_, In Proceedings of the 7th Nordic Signal
  Processing Symposium, 2006. DOI: [10.1109/NORSIG.2006.275210][eog-doi]. The
  article is freely available [here][aar-tut]. A slightly extended version of
  the manuscript is available at [my homepage][aar-extended].

[aar-extended]: http://germangh.com/papers/norsig06.pdf
[eog-doi]: http://dx.doi.org/10.1109/NORSIG.2006.275210
[aar-tut]: http://sp.cs.tut.fi/publications/archive/Gomez-Herrero2006-Automatic.pdf

Please __cite the publication above__ if you use the AAR plug-in in any of your
scientific publications.

The automatic EMG correction method is based on the following reference:

* De Clercq, W. et al., __Canonical Correlation Analysis Applied to Remove
  Muscle Artifacts from the Electroencephalogram__, IEEE Trans. Biomed. Eng 53
  (12), pp. 2583-2587. DOI: [10.1109/TBME.2006.879459][doi-emg].

[doi-emg]: http://dx.doi.org/10.1109/TBME.2006.879459

You can also  get some of the [datasets][datasets] that were used to evaluate
the performance of some of the methods included in the AAR toolbox.

[datasets]: http://germangh.com/datasets/epilepsy




