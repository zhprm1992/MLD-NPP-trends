% Used for: driver for daily MLD from ECCO-Darwin LLC270 T/S, run on the NASA NAS cluster
% Original: Code\mld_ECCO\Compute_MLD\compute_MXLDEPTH_NAS_wrapper.m (data paths updated to I:\Prj; computation unchanged)
clear 
close all;

addpath(genpath('<cluster_dir>/MATLAB'));

dataDir = 'I:\Prj\MLD_ECCO\LLC_270\';
% codeDir = '<cluster_dir>/debugging/m_files/compute_MXLDEPTH/LLC_270/';
% saveDir = '<cluster_dir>/debugging/bin/MXLDEPTH/LLC_270/';
saveDir = 'I:\Prj\MLD_ECCO\LLC_270';

% cd(codeDir);

compute_MXLDEPTH_NAS( dataDir, saveDir );

