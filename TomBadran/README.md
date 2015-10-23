# spire-hires-testing

/doc/ contains output figures and the presentation document given in the telecon.

/stats/ contains the data analysis of the entire set of observations in the archive
    -   {band}.csv contains the full set of statistics calculated
    -   Results.csv is set of information used as the threshold tests including pass/fail

/obs-lists/ contains dataset output from some scripts, and data classifying the various observations by proposal type

/obs-to-hires/ contain lists of the observation ids that pass the final agreed tests

Some scripts are designed to be run in hipe, denoted as (HIPS), others are pure python and use the standard scipy stack with astropy and aplpy where appropriate.

 beam-rotation.py (HIPE) - Run HiRes using a set of rotated beam images

 beam-size-profiles-SPIRE.py - Take a set of fits of files, named by hires beam size and create an interactive plot of a cut through the image - uses files output by beam-size-SPIRE.py and regird.py

 beam-size-SPIRE.py (HIPE) - Perform a hires mapping using beams cropped to various sizes, and output these to fits files

 regrid.py (HIPE) - Regrid nominal maps from beam-size-SPIRE.py to match the pixel grid of the HiRes maps

 beam-size-profiles.py - Similar to beam-size-profiles-spire.py but using truth sources rather than SPIRE observations

 common.py - Some helper functions

 display-thumbs.py - Sets of scripts to generate thumbnail figures for various outputs

 gen-thumbs.py (HIPE) - Script that generates images from lists of obsevation ids

 iter-compare.py - Slice into an observation set to see the difference between hires output for different iterations

 iter_truth.py (HIPE) - Generate hires output for varying iteration counts on spitzer truth sources

 iter_truth_compare.py - Plot the RMS pixel difference for output of iter_truth.py

 iteration-runtime.py (HIPE) - Measure the runtime for HiRes as a function of iteration Number

 iteration-tests.py (HIPE) - Output HiRes data for a an observation at selected iterations

 obs-stats-analysis.py - Perform analysis on the full statistical data set for the HSA

 obs-stats.py (HIPE) - Generate statistics on the full archive

 radial-beam-test.py (HIPE) - Run HiRes on observations using both the fine and radial beams

 regrid.py (HIPE) - Regrid certain nominal maps to be the same as HiRes output
