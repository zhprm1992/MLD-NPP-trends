# Used for: depth-weighted mean N2 (0-200 m) from daily NEMO T/S at 1/12 deg, run on the PACE cluster (Supplementary Table 1)
# Original: I:\Prj\MLD_1M_NEMO\python code\compute_bvf_depth_weighted.py (data paths updated to I:\Prj; computation unchanged)
# compute_bvf_depth_weighted.py
#
# Recomputes 0-200m Brunt-Vaisala frequency squared (N^2 / BVF) with a
# proper DEPTH-WEIGHTED (trapezoidal) vertical average, replacing the
# unweighted level-mean used in mldD03.py (line 24: np.nanmean(N2_col[...])).
#
# Why this matters: CMEMS GLORYS12 (cmems_mod_glo_phy_my_0.083deg_P1D-m)
# has 26 depth levels within 0-200m, spaced from ~1m near the surface to
# ~30m near 200m (a ~30x range). An unweighted mean over these levels
# heavily overweights the shallow ocean. A quick empirical test (coarse
# spatial subsample, 1 day/season, 0-20N band) found depth-weighted means
# ~30% higher than unweighted, and the SIGN of the 20-year Spring trend
# flipped between the two methods. This script redoes the full daily/
# monthly computation properly so that result can be checked rigorously.
#
# Scope (deliberately narrower than mldD03.py to keep this fast):
#   - Only computes bvf_200 (NOT mld_d03 -- MLD itself doesn't have this
#     averaging issue, findMLD_d03 is a threshold crossing, not a level
#     mean, so it's left alone / not recomputed here).
#   - Only fetches depth 0-200m (mldD03.py fetches full depth for MLD;
#     skipping that here makes each day's API request much smaller/faster).
#   - Only latitude -60 to 60 (matches the 6 lat bands used in
#     correlation_bvf_npp_nemo.m / Table S5-S6 of the response letter --
#     no need for polar latitudes for this check).
#
# Usage (run once per year, e.g. from 10 separate terminal sessions):
#   nohup python compute_bvf_depth_weighted.py --year 2003 > log_2003.log 2>&1 &
#
# Output: one NetCDF per month, bvf_200_{year}_{month:02d}.nc, in
# OUT_BASE_DIR/{year}/, containing just:
#   - bvf_200_weighted : depth-weighted (trapezoidal) 0-200m average
# (The unweighted equivalent already exists from the original mldD03.py
# run, so it is NOT recomputed/stored here -- halves both compute and
# storage. Grid is 1441 x 4320 at 0.083deg -> ~23.7 MiB/month/field ->
# ~5.7 GiB total for all 20 years.)
# Re-running is safe -- a month is skipped if its output file already
# exists, so a killed/restarted session resumes where it left off.
#
# DIAGNOSTIC UPDATE 2026-07-27: the 2003-2019 output produced by an earlier
# run of this exact script on the cluster came back with valid file shapes
# but ZERO finite values anywhere (confirmed by direct inspection of the
# saved .nc files). Reproducing this exact script locally (same file,
# byte-for-byte) against real 2003 CMEMS data produced correct, physically
# sensible results -- so the script logic itself is not the bug. The
# leading suspect is the bare `except Exception: return (i, j, np.nan)` in
# fast_process_profile, which silently swallows any error (e.g. a gsw/
# numpy version mismatch between this cluster environment and the one
# tested locally) and converts it to NaN with zero trace. This version
# adds two things to find the real cause on the next run, without changing
# any computed values when things ARE working:
#   1. Package versions printed once at startup.
#   2. The first exception hit by each worker process is printed once
#      (not every time, to avoid flooding the log across millions of
#      profiles) via a per-process flag on the function object.
#   3. A per-day sanity check: if a day's fetch succeeds but the resulting
#      wtd_2d ends up entirely NaN, this is now printed explicitly with
#      SA/CT finite-value counts, so we can tell if the fetched DATA was
#      already empty (an API/fetch problem) vs. the fetch was fine but the
#      per-profile gsw computation failed for every point (a computation
#      problem).
# Also fixes a real (separate) bug: the early-return on an all-zero
# salinity column returned a 4-tuple `(i, j, np.nan, np.nan)` while every
# other path returns a 3-tuple; process_month unpacks with exactly 3
# variables (`i, j, n2_wtd = result`), so if that branch were ever hit it
# would crash with ValueError. Fixed to return a 3-tuple like every other
# path.
#
# Haipeng Zhao
# 2026-07-24, 2026-07-27

import argparse
import calendar
import os
import time

import copernicusmarine as cm
import gsw
import numpy as np
import xarray as xr
from concurrent.futures import ProcessPoolExecutor

# ---------------------------------------------------------------------------
OUT_BASE_DIR = r"<cluster_dir>/brunt_vaisala/bvf_200_depth_weighted"
MAX_WORKERS = 8           # CPU workers for the per-profile N^2 step -- matches the
                          # 8-core PACE OnDemand allocation this is meant to run under
MAX_RETRIES = 5           # API fetch retries before giving up on a day
RETRY_BASE_DELAY_SEC = 30  # exponential backoff base (30s, 60s, 120s, 240s, 480s)
LAT_MIN, LAT_MAX = -60.0, 60.0
DEPTH_MAX = 200.0
# ---------------------------------------------------------------------------

# ROOT CAUSE 2026-07-27, confirmed via cluster [DIAG] traceback: np.trapezoid
# is a NumPy 2.0+ name (replacing the older np.trapz, deprecated but still
# present as an alias in 2.x). The cluster's numpy==1.24.3 predates this
# rename entirely and has no np.trapezoid attribute, so every profile threw
# AttributeError, was silently caught by fast_process_profile's bare
# except, and became NaN -- this is the actual cause of the all-NaN
# 2003-2019 output. np.trapz and np.trapezoid compute the identical
# trapezoidal-rule integral (same algorithm, just renamed), so this shim
# is a pure compatibility fix with no effect on results. Local testing
# (numpy 2.3.5, which has trapezoid) is unaffected either way.
_trapz = getattr(np, 'trapezoid', None) or np.trapz


def fast_process_profile(args):
    sa_col, ct_col, p_col, lat_val, i, j = args

    if not np.any(np.isfinite(sa_col)) or np.all(sa_col == 0):
        return (i, j, np.nan)  # fixed 2026-07-27: was a 4-tuple, mismatched the 3-tuple unpack in process_month

    try:
        N2_col, p_mid_col = gsw.Nsquared(sa_col, ct_col, p_col, lat_val)
        z_mid = -gsw.z_from_p(p_mid_col, lat_val)
        in200 = (z_mid >= 0) & (z_mid <= DEPTH_MAX)
        n2_sub = N2_col[in200]
        z_sub = z_mid[in200]

        if not np.any(np.isfinite(n2_sub)):
            return (i, j, np.nan)

        # Depth-weighted (trapezoidal) mean over 0-200m
        order = np.argsort(z_sub)
        z_o, n2_o = z_sub[order], n2_sub[order]
        valid = np.isfinite(n2_o)
        if valid.sum() < 2:
            n2_weighted = np.nan
        else:
            z_v, n2_v = z_o[valid], n2_o[valid]
            n2_weighted = _trapz(n2_v, z_v) / (z_v[-1] - z_v[0])

        return (i, j, n2_weighted)
    except Exception as e:
        # DIAGNOSTIC 2026-07-27: print the first exception hit by each worker
        # process only (not every occurrence -- there can be millions of
        # profiles/day), so a re-run's log reveals the real cause instead of
        # silently going to NaN. Safe to leave in permanently: negligible
        # cost, and it stays silent once triggered.
        if not getattr(fast_process_profile, '_warned', False):
            import traceback
            print(f"    [DIAG pid={os.getpid()}] fast_process_profile failed at ({i},{j}): "
                  f"{type(e).__name__}: {e}", flush=True)
            traceback.print_exc()
            fast_process_profile._warned = True
        return (i, j, np.nan)


def fetch_with_retry(date_str, executor_workers):
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            ts = cm.open_dataset(
                dataset_id='cmems_mod_glo_phy_my_0.083deg_P1D-m',
                variables=['thetao', 'so'],
                start_datetime=date_str, end_datetime=date_str,
                minimum_latitude=LAT_MIN, maximum_latitude=LAT_MAX,
                minimum_depth=0, maximum_depth=DEPTH_MAX,
            )
            return ts
        except Exception as e:
            wait = RETRY_BASE_DELAY_SEC * (2 ** (attempt - 1))
            print(f"    [retry {attempt}/{MAX_RETRIES}] fetch failed for {date_str}: {e} "
                  f"-- waiting {wait}s", flush=True)
            if attempt < MAX_RETRIES:
                time.sleep(wait)
    print(f"    GIVING UP on {date_str} after {MAX_RETRIES} attempts", flush=True)
    return None


def process_month(year, month, executor):
    out_dir = os.path.join(OUT_BASE_DIR, str(year))
    os.makedirs(out_dir, exist_ok=True)
    outfile = os.path.join(out_dir, f"bvf_200_{year}_{month:02d}.nc")

    if os.path.exists(outfile):
        print(f"  {year}-{month:02d}: already done, skipping", flush=True)
        return

    num_days = calendar.monthrange(year, month)[1]
    sum_wtd = count = None
    lat = lon = None

    for day in range(1, num_days + 1):
        date_str = f"{year}-{month:02d}-{day:02d}"
        t0 = time.time()
        ts = fetch_with_retry(date_str, executor)
        if ts is None:
            continue  # skip this day, keep going -- monthly mean uses whatever days succeeded

        lat = ts['latitude'].values.astype(np.float32)
        lon = ts['longitude'].values.astype(np.float32)
        depth = ts['depth'].values.astype(np.float32)
        z = -depth

        Z, LAT = np.meshgrid(z, lat, indexing='ij')
        pres_2d = gsw.p_from_z(Z, LAT)
        PRES = np.repeat(pres_2d[:, :, np.newaxis], len(lon), axis=2)

        so_vals = ts['so'].values[0].astype(np.float32)
        thetao_vals = ts['thetao'].values[0].astype(np.float32)

        SA = gsw.SA_from_SP(so_vals, PRES,
                             lon[np.newaxis, np.newaxis, :],
                             lat[np.newaxis, :, np.newaxis]).astype(np.float32)
        CT = gsw.CT_from_pt(SA, thetao_vals).astype(np.float32)

        nz, nlat, nlon = SA.shape
        wtd_2d = np.full((nlat, nlon), np.nan)

        args = []
        for i in range(nlat):
            for j in range(nlon):
                sa_col = SA[:, i, j]
                if np.all(np.isnan(sa_col)):
                    continue
                args.append((sa_col, CT[:, i, j], PRES[:, i, j], lat[i], i, j))

        for result in executor.map(fast_process_profile, args, chunksize=5000):
            i, j, n2_wtd = result
            wtd_2d[i, j] = n2_wtd

        # DIAGNOSTIC 2026-07-27: if the fetch succeeded (we got here) but
        # every single profile came back NaN, that's exactly the failure
        # mode seen in the 2003-2019 output -- print fetch-side stats so we
        # can tell whether SA/CT themselves were already empty (fetch/API
        # problem) vs. populated but every gsw call failed (computation
        # problem, should also show up as [DIAG] lines above from the
        # per-worker exception print).
        n_finite = np.isfinite(wtd_2d).sum()
        if n_finite == 0 and len(args) > 0:
            print(f"    [DIAG] {date_str}: fetch OK but ALL {len(args)} attempted profiles "
                  f"came back NaN. SA finite={np.isfinite(SA).sum()}/{SA.size}, "
                  f"CT finite={np.isfinite(CT).sum()}/{CT.size}", flush=True)

        if sum_wtd is None:
            sum_wtd = np.zeros_like(wtd_2d)
            count = np.zeros_like(wtd_2d)

        sum_wtd += np.nan_to_num(wtd_2d)
        count += np.isfinite(wtd_2d)

        print(f"  {date_str}: done in {(time.time()-t0)/60:.2f} min", flush=True)

    if sum_wtd is None:
        print(f"  {year}-{month:02d}: ALL days failed, no output written", flush=True)
        return

    bvf_weighted = np.divide(sum_wtd, count, out=np.full_like(sum_wtd, np.nan), where=(count > 0))

    ds_out = xr.Dataset(
        {
            "bvf_200_weighted": (("lat", "lon"), bvf_weighted),
        },
        coords={"lat": lat, "lon": lon},
    )
    ds_out.to_netcdf(outfile)
    print(f"  Saved {outfile}", flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--year', type=int, required=True,
                         help='Single year to process, e.g. --year 2003')
    args = parser.parse_args()

    print(f"Starting year {args.year} (PID {os.getpid()})", flush=True)
    print(f"Package versions: gsw={gsw.__version__}, numpy={np.__version__}, "
          f"copernicusmarine={cm.__version__ if hasattr(cm, '__version__') else 'unknown'}, "
          f"xarray={xr.__version__}", flush=True)
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        for month in range(1, 13):
            process_month(args.year, month, executor)

    print(f"Finished year {args.year}", flush=True)
