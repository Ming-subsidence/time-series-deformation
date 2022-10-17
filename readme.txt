#BWS-DIE MT InSAR - Time-series InSAR processing flow based on BWS-DIE point selection algorithm

## Description
BWS-DIE MT InSAR is a toolkit developed based on Matlab software, including BWS-DIE homogenous point selection algorithm, phase optimization algorithm based on covariance matrix,
Integrate PS points and DS points for deformation calculation and data export.

## Support environment
Linux system
MATLAB
StaMPS 4.1

## Import Data
The GAMMA software format was selected as the input format, including intensity image sequence (real number) and interference difference sequence (complex number).

## Homogeneous sample selection algorithm - BWS-DIE algorithm
[SHP]=SHP_BWSDIE(mlistack,CalWin,Alpha);
enter
mlistack: SAR intensity image sequence
CalWin: window size [azimuth, slant range]
Alpha: hypothesis test significance level
output
SHP: structure: includes homogeneous sample set, sample number, window size

## Phase optimization
[pcoh,optintf]=optphase(intfstack,SHP);
enter
intfstack: SAR image differential interferometry sequence
SHP: structure: includes homogeneous sample set, sample number, window size
output
pcoh: posterior coherence (goodness-of-fit value)
optintf: Phase-optimized SAR image differential interferometry sequence

## Filter DS points and merge DS points and PS points
[DS,DPS]=selection(SHP,pcoh,pcoh_th,BroNum_ph);
enter
SHP: structure: includes homogeneous sample set, sample number, window size
pcoh: posterior coherence (goodness-of-fit value)
pcoh_th: posterior coherence (goodness-of-fit value) (0.7)
BroNum_ph: Threshold for the number of homogenous pixels (20)
output
DS: DS point coordinates (in StaMPS format)
DPS: Fuse DS point and PS point coordinates (in StaMPS format)

## Data output
[disp,v]=sortout; (under StaMPS working path)
output
disp: los-direction settlement corresponding to the shooting time of the SAR image
v: annual average subsidence rate

## Monte Carlo simulation experiment
[h1,h2,h3,h4,meanh1,stdh1,meanh2,stdh2,meanh3,stdh3,meanh4,stdh4,]=Monte_Carlo(stacksize,CalWin,sigma);
enter
stacksize: number of simulated samples
CalWin: Simulate window size
sigma: analog noise level
output
hi, meanhi, stdhi: Algorithm simulation power value, power mean and power value standard deviation
