# Used for: MLD with the 0.125 kg/m3 threshold, per profile (Supplementary Figs. 1, 2)
# Original: I:\Prj\MLD_1M_NEMO\python code\findMLD_d125.py (data paths updated to I:\Prj; computation unchanged)
import numpy as np
import xarray as xr
import gsw

# --- Reviewer-1 sensitivity check: same algorithm as findMLD_d03, but with a
# 0.125 kg/m^3 sigma0 threshold instead of 0.03. Separate file/function so the
# original 0.03 pipeline (findMLD.py) is untouched. ---
def findMLD_d125( SA, CT, pres ):
    """
    Calculate MLD using 0.125 kg/m^3 sigma0 threshold from 10 dbar.
    Inputs: 1D arrays of SA, CT, pres (depth).
    Returns: mixed layer depth (MLD), or np.nan.
    """
    if np.all( np.isnan( SA ) ) or len( SA ) < 2:
        return np.nan

    starti = np.argmin( (pres - 10) ** 2 )
    SA = SA[starti:]
    CT = CT[starti:]
    pres = pres[starti:]

    if len(SA) < 2 or np.isnan( SA[0] ):
        return np.nan

    sigma0 = gsw.sigma0( SA, CT )
    ref_sigma0 = sigma0[0]

    for j in range( 1, len( sigma0 ) ):
        if abs( sigma0[j] - ref_sigma0 ) > 0.125:
            presseg = [pres[j - 1], pres[j]]
            denseseg = [ref_sigma0 - sigma0[j - 1], ref_sigma0 - sigma0[j]]
            if np.any( np.isnan( denseseg ) ):
                return np.nan
            P = np.polyfit( presseg, denseseg, 1 )
            presinterp = np.arange( presseg[0], presseg[1] + 0.5, 0.5 )
            interp_val = np.polyval( P, presinterp )
            idx = np.where( np.abs( interp_val ) < 0.125 )[0]
            return presinterp[idx[-1]] if len(idx) > 0 else pres[j]
    return np.nan
