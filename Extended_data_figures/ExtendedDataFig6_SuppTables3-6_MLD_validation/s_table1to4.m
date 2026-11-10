% Used for: Supplementary Tables 3-6: product closest to slope 1 in each row (highlighting); run after scatterMLD_ARGOvsModels_point_season.m
% Original: Code\s_table1to4.m (data paths updated to I:\Prj; computation unchanged)
% Used for S. table 1 - 4 in supplementary materials 
% Run this script: scatterMLD_ARGOvsModels_point_season.m

metrics = [error_table_sp.NEMO_Slope, error_table_sp.ECCO_Slope, error_table_sp.HYCOM_Slope];
[~, best_idx_sp] = min(abs(metrics - 1), [], 2);


metrics = [error_table_su.NEMO_Slope, error_table_su.ECCO_Slope, error_table_su.HYCOM_Slope];
[~, best_idx_su] = min(abs(metrics - 1), [], 2);


metrics = [error_table_au.NEMO_Slope, error_table_au.ECCO_Slope, error_table_au.HYCOM_Slope];
[~, best_idx_au] = min(abs(metrics - 1), [], 2);

metrics = [error_table_wi.NEMO_Slope, error_table_wi.ECCO_Slope, error_table_wi.HYCOM_Slope];
[~, best_idx_wi] = min(abs(metrics - 1), [], 2);
