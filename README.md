<p align="center">
  <img width="350" src="https://github.com/PlanetHunters/watson/blob/main/images/watson.png?raw=true">
</p>

<b>WATSON</b> (<b>V</b>isual <b>V</b>etting and <b>A</b>nalysis of <b>T</b>ransits from <b>S</b>pace Observatio<b>N</b>s) is a lightweight software package
that enables a comfortable visual vetting of a transiting signal candidate from Kepler, K2 and TESS missions.


Any transiting candidate signal found in a space-based mission could have been potentially generated by 
different sources or even be instrumental artifacts induced into a target's light curve. 
To rule-out these scenarios, the Science Processing Operation Center (SPOC) of the NASA implemented 
the Data Validation (DV) Reports, which are one or two pages sheets showing different metrics to qualify 
or discard an analyzed candidate. These scenarios are mainly 

* Transit shape model fit
* Odd-even transits checks, 
* Centroids shifts
* Optical ghost effects
* Transit source offsets
* Rolling band contamination histograms

WATSON is also implementing similar but more simplified checks for all of those scenarios 
(SPOC fits transits models and we just compute the SNR of the possible candidate signal) except the 
rolling band contamination. In addition, we included a new check comparing the transit SNRs in the different 
available cadences and also all the single-transit plots computed with the official pipeline aperture and a 
smaller one. With all of these data, we compute metrics that might alert the scientist about problematic 
signals not complying with any of the thresholds.

# WATSON
Visual Vetting and Analysis of Transits from Space ObservatioNs
