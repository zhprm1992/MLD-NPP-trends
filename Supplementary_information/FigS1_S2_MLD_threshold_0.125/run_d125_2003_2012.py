# Used for: driver for calMLD_multiprc_d125.py, 2003-2012 (2013-2022 run on a second computer with the same code)
# Original: I:\Prj\MLD_1M_NEMO\python code\run_d125_2003_2012.py (data paths updated to I:\Prj; computation unchanged)
import os
import time
import traceback
from calMLD_multiprc_d125 import mld_month_nemo_d125

OUT_BASE_DIR = r"I:\Prj\MLD_1M_NEMO\one12th_d125"

# Desktop covers 2003-2012; laptop covers 2013-2022 in parallel
# (see python code\laptop_2013_2022\). Capped here so the two runs don't
# duplicate each other's work.
if __name__ == '__main__':
    for year in range(2003, 2013):
        for month in range(1, 13):
            outfile = os.path.join(OUT_BASE_DIR, str(year), f"mld_d125_{year}_{month:02d}.nc")
            if os.path.exists(outfile):
                print(f"Skipping {year}-{month:02d}, already exists: {outfile}", flush=True)
                continue
            print(f"=== Starting {year}-{month:02d} ===", flush=True)
            t0 = time.time()
            try:
                mld_month_nemo_d125(year, month, out_base_dir=OUT_BASE_DIR)
            except Exception as e:
                print(f"!!! FAILED {year}-{month:02d}: {e}", flush=True)
                traceback.print_exc()
                continue
            print(f"=== Finished {year}-{month:02d} in {(time.time()-t0)/60:.1f} min ===", flush=True)
    print("ALL DONE (2003-2012)", flush=True)
