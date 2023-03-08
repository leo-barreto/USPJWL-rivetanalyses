# USPJWL-rivetanalyses
Rivet analyses for JEWEL with realistic hydrodynamics

All analyses are written for JEWEL's custom version of Rivet 3 with the Constituent Subtraction methodology (see https://jewel.hepforge.org/subtraction.html). 

They were intended for JEWEL coupled with realistic hydro, but they will work for out-of-the-box JEWEL as well. For other MC generators, check and modify the uses of `SubtractedJewelEvent` and `SubtractedJewelFinalState` projections.


---

## Implemented Observables
 - Inclusive jet spectrum/nuclear modification factor $R_{AA}$: USPJWL_JETSPEC, USPJWL_EXTRASPEC.
 - Leading and subleading jet spectrum/nuclear modification factor $R_{AA}$: USPJWL_JETSPEC.
 - Dijet yield $x_J$: USPJWL_JETSPEC.
 - Jet azimuthal distribution (for $v_n$): USPJWL_PHIDIST
 - Leading/inclusive subjet fragmentation: USPJWL_SUBFRAG
 - Jet mass $M_{jet}$ : USPJWL_JET_MASS.
 - Semi-inclusive hadron+jet correlation spectrum: USPJWL_HJET.
