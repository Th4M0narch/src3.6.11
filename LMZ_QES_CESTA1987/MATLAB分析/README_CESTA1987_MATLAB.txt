CESTA 1987 QES-Dense MATLAB analysis suite
===========================================

Target folder
-------------
I:\LMZ\QES\LMZ_QES_CESTA1987\MATLAB分析

The scripts resolve the project path automatically:
  MATLAB分析 -> parent = I:\LMZ\QES\LMZ_QES_CESTA1987
  Output     -> I:\LMZ\QES\LMZ_QES_CESTA1987\Output

Main command
------------
In MATLAB, cd to MATLAB分析 and run:

    results = run_all_CESTA1987;

Main files expected
-------------------
Output\observation_sampler.csv
Output\component_history.csv
Output\parcel_mass_diagnostic.csv
Output\lmz_windsOut.nc
Output\lmz_windsWk.nc
Output\lmz_turbOut.nc
Output\lmz_plumeOut.nc

Generated analysis
------------------
1. Input/source/meteorology consistency summary
2. Direct-entry CESTA wind profile + reference log/power-law fits
3. Component/reaction 20-particle history and conservation
4. A-stage qMass accepted-step residual diagnostics
5. Integrated sampler U/F model-vs-experiment metrics
6. Integrated crosswind profiles and downwind maxima
7. Sequential 6-minute U/F model-vs-experiment metrics
8. Sequential time histories for every distance/azimuth/height group
9. S0 42-min sampler vs duration-weighted S1-S7 comparison
10. Optional plume NetCDF inventory

Performance metrics
-------------------
FAC2, FAC5, MG, VG, FB, NMSE, NAD, CC are written to CSV files under Results.

Important interpretation note
-----------------------------
The current D2 implementation keeps UF6/HF/UO2F2 species in the sparse observation
sampler CSV. It intentionally does not expand the frozen plume NetCDF schema.
Therefore do not label a legacy NetCDF "conc" variable as a post-reaction species field.

CESTA baseline encoded in the supplied XML
------------------------------------------
Source: (550,100,3.15) m
UF6 mass: 146.2 kg
Release: model t=2 to 1807 s = 1805 s
Release rate: 80.9 g/s
Exit velocity: 3.66 m/s
Ambient: 286.15 K, RH=0.722, P=101300 Pa
Stability: Class C
Direct wind profile:
  2 m  -> 2.8 m/s
 10 m  -> 3.3 m/s
 18 m  -> 3.6 m/s
 28 m  -> 3.8 m/s
QES wind-from direction: 122.7 deg
CESTA downwind/plume azimuth: 302.7 deg
