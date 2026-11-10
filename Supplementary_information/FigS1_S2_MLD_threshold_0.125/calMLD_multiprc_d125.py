# Used for: monthly MLD-NEMO with the 0.125 threshold from daily 1/12 deg T/S (SI Figs. S1, S2)
# Original: I:\Prj\MLD_1M_NEMO\python code\calMLD_multiprc_d125.py (data paths updated to I:\Prj; computation unchanged)
import copernicusmarine as cm
import xarray as xr
import gsw
import numpy as np
import os
import calendar
from findMLD_d125 import findMLD_d125
from concurrent.futures import ProcessPoolExecutor
import time


def fast_process_profile(args):
    # bvf/N2 deliberately not recomputed here: it's threshold-independent (same
    # profiles, same gsw.Nsquared calc as the 0.03 pipeline), so the existing
    # one_d03/bvf output already covers it -- only MLD depends on the threshold.
    sa_col, ct_col, p_col, lat_val, i, j = args
    if not np.any(np.isfinite(sa_col)) or np.all(sa_col == 0):
        return (i, j, np.nan)
    try:
        mld_val = findMLD_d125(sa_col, ct_col, p_col)
        return (i, j, mld_val)
    except Exception:
        return (i, j, np.nan)


def mld_month_nemo_d125(year, month, out_base_dir=r"."):
    data_request = dict(
        dataset_id='cmems_mod_glo_phy_my_0.083deg_P1D-m',
        variables=['thetao', 'so']
    )

    num_days = calendar.monthrange(year, month)[1]
    sum_mld = None

    with ProcessPoolExecutor(max_workers=20) as executor:
        for day in range(1, num_days + 1):
            date_str = f"{year}-{month:02d}-{day:02d}"
            print(f"Processing {date_str}", flush=True)
            start_time = time.time()
            try:
                ts = cm.open_dataset(**data_request, start_datetime=date_str, end_datetime=date_str)
            except Exception as e:
                print(f"    Skipping {date_str} due to error: {e}", flush=True)
                continue

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
            mld = np.full((nlat, nlon), np.nan)

            args = [(SA[:, i, j], CT[:, i, j], PRES[:, i, j], lat[i], i, j)
                    for i in range(nlat) for j in range(nlon)
                    if not np.all(np.isnan(SA[:, i, j]))]

            for result in executor.map(fast_process_profile, args, chunksize=10000):
                i, j, mld_val = result
                mld[i, j] = mld_val

            if sum_mld is None:
                sum_mld = np.zeros_like(mld)
                count_mld = np.zeros_like(mld)

            sum_mld += np.nan_to_num(mld)
            count_mld += np.isfinite(mld)

            elapsed_min = (time.time() - start_time) / 60
            print(f"    Finished {date_str} in {elapsed_min:.2f} minutes", flush=True)

    mld_d125 = np.divide(sum_mld, count_mld, out=np.full_like(sum_mld, np.nan), where=(count_mld > 0))

    month_dir = os.path.join(out_base_dir, str(year))
    os.makedirs(month_dir, exist_ok=True)

    ds_out = xr.Dataset(
        {
            "mld_d125": (("lat", "lon"), mld_d125)
        },
        coords={
            "lat": lat,
            "lon": lon
        }
    )
    outfile = os.path.join(month_dir, f"mld_d125_{year}_{month:02d}.nc")
    ds_out.to_netcdf(outfile)
    print(f"Saved monthly average to {outfile}", flush=True)
