# Used for: daily 1/12 deg NEMO (CMEMS GLOBAL_MULTIYEAR_PHY_001_030, cmems_mod_glo_phy_my_0.083deg_P1D-m) T/S -> MLD per day -> monthly mean, mld_N2_YYYY_MM.nc (one12th_d03); the year and month loop was set per run (2003-2022)
# Original: I:\Prj\MLD_1M_NEMO\python code\mldD03.py (data paths updated to I:\Prj; computation unchanged)
# Calculate MLD based on density 0.03 threshold
# Haipeng Zhao
# 5/10/2025

import copernicusmarine as cm
import xarray as xr
import gsw  # TEOS-10 library for seawater calculations
import numpy as np
import os
import calendar
from findMLD import findMLD_d03
from concurrent.futures import ProcessPoolExecutor
import time

def fast_process_profile(args):
    sa_col, ct_col, p_col, lat_val, i, j = args

    if not np.any(np.isfinite(sa_col)) or np.all(sa_col == 0):
        return (i, j, np.nan, np.nan)

    try:
        N2_col, p_mid_col = gsw.Nsquared(sa_col, ct_col, p_col, lat_val)
        z_mid = -gsw.z_from_p(p_mid_col, lat_val)
        N2_avg = np.nanmean(N2_col[(z_mid >= 0) & (z_mid <= 200)])
        mld_val = findMLD_d03(sa_col, ct_col, p_col)
        return (i, j, N2_avg, mld_val)
    except Exception:
        return (i, j, np.nan, np.nan)

g = 9.81  # Gravity (m/s^2)
data_request = dict(
    dataset_id='cmems_mod_glo_phy_my_0.083deg_P1D-m',
    variables=['thetao', 'so']
)

out_base_dir = r"F:\\MLD_NEMO\\d03_1M_008dg"

if __name__ == '__main__':
    executor = ProcessPoolExecutor(max_workers = 10)
    for year in range(2003, 2004):  # Loop over years
        for month in range(11, 13):  # Loop over months
            num_days = calendar.monthrange(year, month)[1]
            sum_N2 = None
            sum_mld = None

            for day in range(1, num_days + 1):  # Loop over days
                start_time = time.time()
                date_str = f"{year}-{month:02d}-{day:02d}"
                print(f"Processing {date_str}")
                try:
                    ts = cm.open_dataset(**data_request, start_datetime=date_str, end_datetime=date_str)
                except Exception as e:
                    print(f"    Skipping {date_str} due to error: {e}")
                    continue

                # Convert variables to float32 to save memory
                lat = ts['latitude'].values.astype(np.float32)
                lon = ts['longitude'].values.astype(np.float32)
                depth = ts['depth'].values.astype(np.float32)
                z = -depth

                Z, LAT = np.meshgrid(z, lat, indexing='ij')
                pres_2d = gsw.p_from_z(Z, LAT)
                PRES = np.repeat(pres_2d[:, :, np.newaxis], len(lon), axis=2)  # shape: (depth, lat, lon)

                # Extract 3D salinity & temperature (time index 0)
                so_vals = ts['so'].values[0].astype(np.float32)
                thetao_vals = ts['thetao'].values[0].astype(np.float32)

                # Compute Absolute Salinity and Conservative Temperature
                SA = gsw.SA_from_SP(so_vals, PRES,
                                    lon[np.newaxis, np.newaxis, :],
                                    lat[np.newaxis, :, np.newaxis]).astype(np.float32)

                CT = gsw.CT_from_pt(SA, thetao_vals).astype(np.float32)
                nz, nlat, nlon = SA.shape
                N2_200 = np.full((nlat, nlon), np.nan)
                mld = np.full((nlat, nlon), np.nan)

                # Prepare arguments for parallel processing
                args = []
                for i in range(nlat):
                    for j in range(nlon):
                        sa_col = SA[:, i, j]
                        if np.all(np.isnan(sa_col)):
                            continue  # skip invalid column
                        ct_col = CT[:, i, j]
                        p_col = PRES[:, i, j]
                        lat_val = lat[i]
                        args.append((sa_col, ct_col, p_col, lat_val, i, j))

                size = round(len(args) / executor._max_workers)
                # Process using parallel processing
                for result in executor.map( fast_process_profile, args, chunksize = 10000 ):
                    i, j, n2_val, mld_val = result
                    N2_200[i, j] = n2_val
                    mld[i, j] = mld_val

                if sum_N2 is None:
                    sum_N2 = np.zeros_like(N2_200)
                    sum_mld = np.zeros_like(mld)
                    count_N2 = np.zeros_like(N2_200)
                    count_mld = np.zeros_like(mld)

                sum_N2 += np.nan_to_num(N2_200)
                sum_mld += np.nan_to_num(mld)
                count_N2 += np.isfinite(N2_200)
                count_mld += np.isfinite(mld)

                end_time = time.time()
                elapsed_min = (end_time - start_time) / 60
                print(f"    Finished {date_str} in {elapsed_min:.2f} minutes")

            # Calculate monthly average of buoyancy frequency (N2) and mixed layer depth (MLD)
            bvf_200 = np.divide(sum_N2, count_N2, out=np.full_like(sum_N2, np.nan), where=(count_N2 > 0))
            mld_d03 = np.divide(sum_mld, count_mld, out=np.full_like(sum_mld, np.nan), where=(count_mld > 0))

            # Save the output to NetCDF
            month_dir = os.path.join(out_base_dir, str(year))
            os.makedirs(month_dir, exist_ok=True)

            ds_out = xr.Dataset(
                {
                    "bvf_200": (("lat", "lon"), bvf_200),
                    "mld_d03": (("lat", "lon"), mld_d03)
                },
                coords={
                    "lat": lat,
                    "lon": lon
                }
            )
            outfile = os.path.join(month_dir, f"mld_N2_{year}_{month:02d}.nc")
            ds_out.to_netcdf(outfile)
            print(f"  Saved monthly average to {outfile}")

    executor.shutdown()