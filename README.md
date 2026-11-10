# Code for "Conflicting trends in mixed layer depth obscure changes in global ocean productivity"

Haipeng Zhao, Nicolas Cassar, M. Susan Lozier (manuscript in revision, Nature Climate Change)

This repository contains the scripts used to compute the mixed layer depth (MLD) products, run the productivity models, and produce every figure and table of the manuscript, its Supplementary Information, and the two response letters. Each script is the version that produced the published result. Apart from data paths, each file is unchanged, except where its two-line header says otherwise.

Every script starts with two comment lines:

```
% Used for: <figure or table, and what the script produces>
% Original: <path of the script in the authors' working folders>
```

## Software

- MATLAB R2023b, with the Statistics and Machine Learning Toolbox (`fitlme`, `fitlm`) and the Mapping Toolbox (`shaperead`).
- Python 3.12 for the NEMO MLD computation and two figure scripts (packages in `1_MLD_computation/NEMO/requirements.txt`: `copernicusmarine`, `xarray`, `gsw`, `numpy`; the figure scripts also use `pandas`, `scipy`, `h5py`, `matplotlib`).
- Third-party MATLAB code, not included here (please obtain from the original sources):
  - `ktaub` (Mann-Kendall test and Sen's slope, J. Burkey, MATLAB File Exchange)
  - `m_map` (R. Pawlowicz), which also provides `m_topo.mat`
  - Gibbs SeaWater (GSW) Oceanographic Toolbox (TEOS-10)
  - `gcmfaces` (MITgcm/ECCO), for the ECCO LLC270 grid
  - `crameri` (Scientific colour maps), `slanCM`, `brewermap` (colour maps)
  - `dscatter` (density scatter plot), `heatmap` (A. Deoras), `radar_chart` (radar chart plot), all from MATLAB File Exchange

## Data

The scripts read gridded data from one data root, written in the scripts as `I:\Prj\...`. To run them elsewhere, replace `I:\Prj` with your own data folder. Input data sources:

| Data | Source |
|---|---|
| NEMO (GLORYS12V1) daily temperature and salinity, 1/12 deg | Copernicus Marine Service, GLOBAL_MULTIYEAR_PHY_001_030, https://doi.org/10.48670/moi-00021 |
| ECCO-Darwin LLC270 daily temperature and salinity | ECCO (D. Carroll, H. Zhang); https://ecco.jpl.nasa.gov/drive/files/ECCO2 |
| HYCOM GOFS 3.1 3-hourly temperature and salinity (GLBv0.08, GLBy0.08) | https://www.hycom.org/ (THREDDS server) |
| Argo MLD, monthly 1 deg (MILA GPV) | JAMSTEC, https://www.jamstec.go.jp/argo_research/dataset/milagpv/mila_en.html |
| Ocean color (Chl, Kd490, PAR, bbp, IOPs, SST) | NASA Ocean Color, https://oceancolor.gsfc.nasa.gov/ |
| Nitrate climatology | World Ocean Atlas 2023, https://www.ncei.noaa.gov/products/world-ocean-atlas |
| Open-ocean biomes | Fay and McKinley (2014), Earth System Science Data 6, 273-284 (shapefile `meanBiom_FAY.shp`) |
| Photosynthesis-irradiance parameters (MAPPS) | Bouman et al. (2018), Earth System Science Data 10, 251-266; the subset used is `Extended_data_figures/FigS11_PE_curve/mapps_masked.csv` |
| BATS and HOT bottle data | BATS (BCO-DMO) and HOT-DOGS; read by `Response_letter_round1/FigR11_BATS_HOT_distance` |

## Folder structure and run order

The folders follow the order of the analysis. Scripts that write intermediate files (`.mat`, `.csv`) save them next to themselves or under the data root, and later scripts read them from there.

To run a script, add `helper_functions/` and the third-party toolboxes to the MATLAB path, then run the script from its own folder. Do not add the whole repository with `genpath`: several folders contain different scripts with the same name (for example one `latSeason_ot.m` per MLD product).

1. `1_MLD_computation/`: MLD from temperature and salinity profiles with a 0.03 kg m-3 density threshold referenced to 10 m (TEOS-10), then regridded to 1 deg.
   - `NEMO/`: `findMLD.py` (one profile) and `mldD03.py` (daily 1/12 deg fields to monthly means; the year and month range was set for each run). Then `prj_grid008to1GIdx_ot.m` (regridding index) and `prj_grid008to1G_ot_v2.m` (1 deg monthly MLD).
   - `ECCO/`: daily MLD on the LLC270 grid computed on the NASA NAS cluster (`compute_MXLDEPTH_NAS_wrapper.m`, `compute_MXLDEPTH_NAS.m`), monthly means (`monthly_ot.m`), and regridding (`proj_llc270to1Gidx.m`, `proj_llc270to1G_ot.m`).
   - `HYCOM/`: one script per year downloads the 3-hourly fields and computes monthly MLD at 1/12 deg (`mld_monthly_2003_ot.m` on a PC, `mld_monthly_2010_ot.m` on the PACE cluster; the year is set in the file). Regridding for 2003 to 2019 (GLBv0.08) and 2020 to 2022 (GLBy0.08) uses `proj_GLBv008to1G*.m` and `proj_GLBy008to1G*.m`.
   - Argo MLD is used as distributed by JAMSTEC.
2. `2_NPP_models/`: CbPM, CAFE and VGPM net primary production per 1 deg cell and month, integrated over the euphotic zone. CbPM and CAFE were run once for each MLD product by changing the MLD input and output folders; the header of each script lists the settings. `biomes_*` scripts split the output into the open-ocean biomes.
3. `3_latitude_and_biome_series/`: MLD and NPP by biome (`biomes_1M_ot.m`, `biomes_1Y_ot.m`), latitudinal series by season and month (`latSeason_ot.m`, `seasonLat_ot.m`, `monthlyLat_*_ot.m`), and per-cell Sen's slopes of annual MLD (`sen_yearly_ot.m`).
4. `4_trends_and_sensitivity/`: trends per latitude row and month (`fs6_heatmap_monthlyNPPtrends_ot.m`, which also draws Fig. S6) and the sensitivity of NPP to MLD (`senMLD2_ot.m`, Eq. 13).
5. Figure folders, below. Shared helper functions are in `helper_functions/`.

## Figures and tables

Figure numbers follow the revised manuscript. Extended Data numbers may change at publication.

| Figure or table | Folder | Main script(s) |
|---|---|---|
| Fig. 1 | none (schematic, drawn) | |
| Fig. 2 | `Main_figures/Fig2_MLD_trends` | `map_senTrends.m` (maps, run once per product; settings in the header); latitudinal bars from `Main_figures/Fig5_Argo_comparison` |
| Fig. 3 | `Main_figures/Fig3_latitude_season_trends` | `f2_plot_combine_NPP_{nemo,ecco,hycom}.m` (one row each) |
| Fig. 4 | `Main_figures/Fig4_sensitivity` | `f4_plot_sensitivity_NPP_MLD.m` (after `4_trends_and_sensitivity`) |
| Fig. 5 | `Main_figures/Fig5_Argo_comparison` | `LAT_vs_MLD_mkTrends_sub_full_annual_v2.m` |
| Table 1 | none (text) | |
| Fig. S1 | `Extended_data_figures/FigS1_CAFE_NPP_vs_aph443` | `LAT_vs_npp_aph443_group.m` |
| Fig. S2 | `Extended_data_figures/FigS2_CAFE_latitude_season_trends` | `f2_plot_combine_NPP_{nemo,ecco,hycom}.m` |
| Fig. S3 | `Extended_data_figures/FigS3_CAFE_month_latitude_trends` | `heatmap_monthlyNPPtrends_ot.m` |
| Fig. S4 | `Extended_data_figures/FigS4_Argo_sub_vs_full_sampled` | `fs4_lat_mld_sub_fully.m` |
| Fig. S5, Tables S1-S4 | `Extended_data_figures/FigS5_TablesS1-S4_MLD_validation` | `scatterMLD_ARGOvsModels_point_season.m` (tables), `s_table1to4.m` (highlighting), `spiderplot_v3.m` (Fig. S5) |
| Fig. S6 | `4_trends_and_sensitivity` | `fs6_heatmap_monthlyNPPtrends_ot.m` (NPP row) |
| Fig. S7 | `Extended_data_figures/FigS7_grid_cell_scatter` | `compute_grid_seasonal_trends.m`, `plot_grid_mld_npp_nemo_dscatter_v5.m` |
| Fig. S8 | `Extended_data_figures/FigS8_biome_sensitivity` | `compute_grid_monthly_trends.m`, `compute_biome_sensitivity.m`, `plot_biome_sensitivity_iqr('sig005')` |
| Figs. S9, S10 | `Extended_data_figures/FigS9_S10_nitracline_distance` | `map_nitracline_mld_pacific_pc_tall.m` (S9), `plot_nitracline_par_mld_climatology_range_combined_v3.m` (S10) |
| Fig. S11 | `Extended_data_figures/FigS11_PE_curve` | `compute_par_mld_mldonly.m`, `fig_pe_final_markers_ownmld_v2.py` |
| SI Table S1 | `Supplementary_information/TableS1_N2_vs_NPP` | `correlation_bvf_npp_nemo_weighted_compare.m` |
| SI Figs. S1, S2 | `Supplementary_information/FigS1_S2_MLD_threshold_0.125` | `plot_mldLat_trend_comparison_d125_vs_003_seasonal.m`, `f4_plot_sensitivity_NPP_MLD_d125_vs_003.m` |
| SI Figs. S3-S5 | `Supplementary_information/FigS3_S4_S5_ensemble_NPP` | `fs3_plot_sensitivity_NPP_MLD_ensemble.m`, `plot_grid_mld_npp_nemo_ensemble_dscatter_v5.m`, `plot_biome_sensitivity_iqr_ensemble('sig005')` |
| First response, Figs. R1, R2 | `Response_letter_round1/FigR1_R2_phytoplankton_carbon` | `f4_plot_sensitivity_Carbon_MLD.m`, `LAT_vs_carbon_mkTrends_allSeasons_nemo.m` |
| First response, Figs. R3-R10 | same as SI Figs. S1-S3 and Figs. S7-S11 | |
| First response, Fig. R11 | `Response_letter_round1/FigR11_BATS_HOT_distance` | `fig_nitracline_mld_distance_bats_hot.py` |
| Second response, Fig. R1 | `Response_letter_round2/FigR1_Fig4_at_p0.05` | `run_fig4_p005.m` |
| Second response, Figs. R2, R4 | `Response_letter_round2/FigR2_R4_interannual_sensitivity` | `plot_trend_vs_interannual_v2.m` (R2), `plot_trend_vs_interannual_biome.m` (R4) |
| Second response, Fig. R3 | `Response_letter_round2/FigR3_nitracline_distance` | `plot_distance_r9style_onepanel.m` |
| Second response, Fig. R5 | `Response_letter_round2/FigR5_Argo_sampling` | `fig_s4_seasonal_nc.m` |

Some figures were finished by hand after export (panel letters, legend position, titles removed, rows of Fig. 3 and Fig. S2 stacked, biome map added to Fig. S8). The data shown are those produced by the scripts.

## Notes

- Several scripts were reused for each MLD product or variable by editing their input and output paths. The repository copy holds the settings used for one product, and its header lists the others.
- The CAFE script (`2_NPP_models/CAFE/grid_npp_1M_1dg_toED_ot_modelMLD.m`) is the version that produced the stored CAFE NPP. Line 266 is marked: a later change to that line (light below the mixed layer, made on 4 November 2026) was not used for any result.
- `error_metrics_season.mat`, read by `spiderplot_v3.m`, was saved from the MATLAB workspace after running `scatterMLD_ARGOvsModels_point_season.m`. It holds the four seasonal tables (`error_table_sp`, `error_table_su`, `error_table_au`, `error_table_wi`) and their medians over 2003 to 2022 used in Fig. S5 (`error_sp = varfun( @median, error_table_sp )`, likewise for `su`, `au`, `wi`; these lines are commented out near the end of that script).
- `mapps_masked.csv` is the MAPPS table with a column `inmask` marking samples inside the study mask; it is provided as the input of Fig. S11.
- Large intermediate files are not part of this repository: the daily ECCO LLC270 MLD (computed on the NASA NAS cluster), the HYCOM monthly MLD at 1/12 deg, and the daily NEMO downloads. The 1 deg monthly MLD and NPP fields can be shared on request and will be archived with the final version of this repository.
- The scripts in `Response_letter_round1` and `Response_letter_round2`, and several analysis scripts written during the revision (grid-cell, biome, nitracline and interannual analyses), were developed with the help of an AI coding assistant under the authors' direction and checked against the published figures.

## Contact

Haipeng Zhao
