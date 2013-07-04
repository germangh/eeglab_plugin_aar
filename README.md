AAR plug-in for EEGLAB
=================

This repository stores the code of the AAR plug-in for [EEGLAB][eeglab] that was
released on 31-10-2008. I am not developing this code anymore but if you report
any major bug I will try to fix it whenever I have time. For suggestions,
comments and bug reports, please contact [German Gomez-Herrero][ggh].


[ggh]: http://germangh.com


## Overview

The AAR plug-in is a collection of [MATLAB][matlab] scripts that implement
several state-of-art (in 2008 anyways...) methods for automatic correction of
ocular and muscular artifacts in the [EEG][eeg]. The capabilities of the toolbox
are briefly described in a fairly outdated tutorial document, which you can get
in [.pdf][tut-pdf], or in [html][tut-html] format.

[matlab]: http://www.mathworks.nl/products/matlab/
[eeg]: http://en.wikipedia.org/wiki/Electroencephalography
[tut-pdf]: http://germangh.com/pubs/aardoc07.pdf
[tut-html]: http://germangh.com/aar/aardoc/

The toolbox is implemented as an [EEGLAB][eeglab] plug-in, but can also be used
without EEGLAB if you don't need or want to use EEGLAB's GUI. The current
version implements several fully automatic methods to correct ocular ([EOG][eog])
artifacts, and one automatic method to correct muscle ([EMG][emg]) artifacts.


[eog]: http://emedicine.medscape.com/article/1140247-overview#aw2aab6b3

[emg]: http://emedicine.medscape.com/article/1140247-overview#aw2aab6b3


## Installation instructions

1. Install [EEGLAB][eeglab] for MATLAB, if you haven't done so already. Ensure
   that EEGLAB is in your MATLAB's search path.

2. Copy and paste the following code in the MATLAB command window and press
   `Enter`:

    ````matlab
    eeglabRoot = fileparts(which('eeglab'));
    url = 'https://github.com/germangh/eeglab_plugin_aar/archive/master.zip';
    unzip(url, [eeglabRoot filesep 'plugins']);
    eeglab redraw;
    ````
[git]: http://git-scm.com/
[masterzip]: https://github.com/germangh/eeglab_plugin_aar/archive/master.zip
[eeglab]: http://sccn.ucsd.edu/eeglab/

3. Optionally, you may want to install additional [BSS][bss] algorithms. Both
   [JADE][jade] and [FastICA][fastica] will be automatically detected by the AAR
   plug-in as long as they are in the MATLAB's search path:

[jade]: http://www.tsi.enst.fr/~cardoso/Algo/Jade/jadeR.m
[fastica]: http://www.cis.hut.fi/projects/ica/fastica/


## Additional resources

The [BSS][bss]-based EOG correction procedure is based on the following
scientific publication:

[bss]: http://en.wikipedia.org/wiki/Blind_signal_separation

* [Gomez-Herrero, G.][ggh] et al., _Automatic Removal of Ocular Artifacts in the
  EEG without an EOG Reference Channel_, In Proceedings of the 7th Nordic Signal
  Processing Symposium, 2006. [10.1109/NORSIG.2006.275210][eog-doi]. The
  article is freely available from [TUT's archives][aar-tut], and from
  [my homepage][aar-home].

[aar-home]: http://germangh.com/papers/norsig06.pdf
[eog-doi]: http://dx.doi.org/10.1109/NORSIG.2006.275210
[aar-tut]: http://sp.cs.tut.fi/publications/archive/Gomez-Herrero2006-Automatic.pdf

The automatic EMG correction method is based on the following reference:

* De Clercq, W. et al., _Canonical Correlation Analysis Applied to Remove
  Muscle Artifacts from the Electroencephalogram_, IEEE Trans. Biomed. Eng 53
  (12), pp. 2583-2587. [10.1109/TBME.2006.879459][doi-emg].

[doi-emg]: http://dx.doi.org/10.1109/TBME.2006.879459

You can also  get some of the [datasets][datasets] that were used to evaluate
the performance of some of the methods included in the AAR toolbox.

[datasets]: http://germangh.com/datasets/epilepsy


If I find major bugs or have important announcements, I will post them to
[this Google group][group]. Please join the group if you are using the AAR
plug-in. Only I can post to the forum so you can be sure of receiving emails
very rarely, if you ever receive any.

[group]: https://groups.google.com/forum/#!forum/aartoolbox

## Known issues

- The algorithm for EMG correction which is based on the criterion
`emg_psd` requires MATLAB's Signal Processing Toolbox v.6.2 or newer.
We do not expect to solve this issue in the near future.

- There exist small differences between the correction results obtained
under MATLAB v7.4 and Signal Processing Toolbox v6.6 and the results
obtained under previous MATLAB releases. Nevertheless, the differences
found so far are very small (negligible with respect to typical EEG
noise levels). The probable cause are the changes that were introduced
in MATLAB's SPT toolbox v6.6.

- There exist very small differences between the correction results
obtain under MATLAB v7.4 for Windows and MATLAB v7.4 for Linux. Again,
the differences are well below typical EEG noise levels. The causes of
these differences are unknown and I have no plans of investigating this further.


## Version history

See [version history](./version_history.md).


## Credit to third parties

See [credits][credits].

[credits]: ./credits.md


## License

The AAR plug-in is released under the
[Creative Commons Attribution-NonCommercial-ShareAlike licence](http://creativecommons.org/licenses/by-nc-sa/3.0/).
Note that third-party dependencies shipped together with the AAR plug-in may
have their own licenses. If you use this software in any of your publications
you must cite the following article:


* [Gomez-Herrero, G.][ggh] et al., _Automatic Removal of Ocular Artifacts in the
  EEG without an EOG Reference Channel_, In Proceedings of the 7th Nordic Signal
  Processing Symposium, 2006. [10.1109/NORSIG.2006.275210][eog-doi]. The
  article is freely available from [TUT's archives][aar-tut], and from
  [my homepage][aar-home].


