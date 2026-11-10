# Used for: Fig. S11 (panel letters added by hand)
# Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\session_2026-09-12_scripts\fig_pe_final_markers_ownmld_v2.py (data paths updated to I:\Prj; computation unchanged)
"""Subtropical P-E curves with the study domain marked. Version 2.

Identical to fig_pe_final_markers.py except for the depth selection of the
MAPPS samples behind the median curve.

WHY     Version 1 kept samples from the surface to a fixed 30 m, close to the
        summer mixed layer of the 20-40 band (median 24 m) and about half its
        annual-mean mixed layer (median 54 m). Agreed with Susan and Nicolas
        (2026-09-17) to use a climatological MLD instead.

CURVES  MAPPS (Bouman et al. 2018). Each record is one water sample, so
        one curve, fully specified by its alpha and PBmax:
        P = PBmax * tanh(alpha*I/PBmax). Selection: the paper's biome +
        topo mask, 20 <= |lat| < 40, and sample depth <= that sample's own
        climatological MLD, where
            own MLD = MLD_clim[sample month, nearest 1 deg cell]
        MLD_clim (par_mld_mldonly.mat, from compute_par_mld_mldonly.m) is the
        2003-2022 mean of NEMO one_d03 monthly MLD per calendar month and
        cell, the same MLD behind the par_mld markers. n = 163 (v1: 139),
        median Ek = 226 (v1: 248).
        Caveats: climatology, not the MLD on the sampling day (samples are
        1979-2012, 69 of 382 within 2003-2022); monthly, 1 deg resolution;
        hard cutoff at the MLD.

MARKERS Unchanged from v1. Climatological par_mld for the same 20-40 band,
        converted per cell and month to a daytime-mean instantaneous value,
        then seasonal mean (hemisphere aware), then PER-CELL MEDIAN over the
        band:
            Winter 68   Autumn 256   Spring 249   Summer 470
        umol quanta m-2 s-1.

Paths are relative to this file, so it runs on either machine without the
portable drive.
"""
import os
import numpy as np, pandas as pd, h5py, scipy.io as sio, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16, 'axes.labelsize': 16,
                     'xtick.labelsize': 14, 'ytick.labelsize': 14})

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                                   # code for revision
OUT = HERE

KMLD = 1.94
XMAX = 800
XMAX_B = 550
SEASONS = ['Winter', 'Autumn', 'Spring', 'Summer']
INST_VALS = {'Winter': 68.0, 'Autumn': 256.0, 'Spring': 249.0, 'Summer': 470.0}
COL = {'Spring': '#C87F00', 'Summer': '#009E73', 'Autumn': '#4D4D4D', 'Winter': '#0072B2'}
INST = INST_VALS

# ---------------- each sample's own climatological MLD ---------------------
M = sio.loadmat(os.path.join(ROOT, 'FigS7_grid_cell_scatter', 'topo_biome_mask_nemo.mat'))
glat = M['lat'].ravel().astype(float)
glon = M['lon'].ravel().astype(float)
with h5py.File(os.path.join(OUT, 'par_mld_mldonly.mat'), 'r') as f:
    mld_clim = np.array(f['MLD_clim'], dtype=float)            # (12 month, 360 lon, 180 lat)
    assert np.allclose(np.array(f['lat']).ravel(), glat), 'latitude vectors differ'

d = pd.read_csv(os.path.join(HERE, 'mapps_masked.csv'))
band = d[d.inmask & (d.Latitude.abs() >= 20) & (d.Latitude.abs() < 40)].copy()
month = pd.to_datetime(band['Date/Time'], errors='coerce').dt.month.values
assert np.all(np.isfinite(month)), 'sample without a date'
ilat = np.abs(band.Latitude.values[:, None] - glat[None, :]).argmin(axis=1)
ilon = np.abs((band.Longitude.values[:, None] - glon[None, :] + 180) % 360 - 180).argmin(axis=1)
band['mld_own'] = mld_clim[month.astype(int) - 1, ilon, ilat]
assert np.all(np.isfinite(band.mld_own)), 'sample without a climatological MLD'

s = band[band.depth <= band.mld_own].copy()
ek = s.Ek.values
EK = np.median(ek)

I = np.linspace(0.5, XMAX, 500)
lI = np.log(I)
nz = nt = 260
zf = (np.arange(nz) + .5) / nz
tf = (np.arange(nt) + .5) / nt
ker = np.outer(np.exp(-KMLD * zf), np.sin(np.pi * tf))
ker = ker / ker.mean()
_x = np.logspace(-4, 3, 4000)
_f = np.tanh(np.multiply.outer(_x, ker)).mean(axis=(1, 2))
_f = _f / _f[-1]

ratio = I[None, :] / ek[:, None]
inst = np.tanh(ratio)
day = np.interp(ratio, _x, _f)
e_i = np.gradient(np.log(inst), lI, axis=1)
e_d = np.gradient(np.log(day), lI, axis=1)
med_inst, med_day = np.median(inst, axis=0), np.median(day, axis=0)
med_ei, med_ed = np.median(e_i, axis=0), np.median(e_d, axis=0)

fig, ax = plt.subplots(1, 2, figsize=(15.0, 6.8))

ORDER = ['Winter', 'Spring', 'Autumn', 'Summer']


def marks(a, label=False):
    h = []
    for k in ORDER:
        ln, = a.plot([INST[k], INST[k]], [0, 1.0], color=COL[k], lw=1.6, alpha=.85,
                     zorder=1, label=('%s (%d)' % (k, round(INST[k])) if label else None))
        h.append(ln)
    return h


# ---------------- (a) the median curve ----------------------------------
ax[0].plot(I, med_inst, color='#B03030', lw=3.0, zorder=5)
hseason = marks(ax[0], label=True)
ax[0].set_xlim(0, XMAX)
ax[0].set_ylim(0, 1.0)
ax[0].set_yticks(np.arange(0, 1.01, .2))
ax[0].set_xlabel('Irradiance ($\\mu$mol quanta m$^{-2}$ s$^{-1}$)')
ax[0].set_ylabel('Production / $P^B_{max}$')
ax[0].set_title('Subtropical P-E curve (20° to 40° latitude)', fontsize=17)

# ---------------- (b) responsiveness ------------------------------------
ax[1].plot(I, med_ei, color='#B03030', lw=2.8, zorder=3)
marks(ax[1])
ax[1].set_xlim(0, XMAX_B)
ax[1].set_ylim(0, 1.0)
ax[1].set_yticks(np.arange(0, 1.01, .2))
ax[1].set_xlabel('Irradiance ($\\mu$mol quanta m$^{-2}$ s$^{-1}$)')
ax[1].set_ylabel('Change in production (%)\nper change in light (%)')
ax[1].set_title('Sensitivity of production to a change in light', fontsize=17)

fig.legend(hseason,
           ['par_mld %s' % k.lower() for k in ORDER],
           loc='lower center', ncol=4, frameon=False, fontsize=15,
           bbox_to_anchor=(.5, .012))
fig.tight_layout(rect=[0, .105, 1, .96])
out_png = os.path.join(OUT, 'fig_pe_subtropical_marked_v2.png')
fig.savefig(out_png, dpi=200)
print('saved %s' % out_png)
print('in-band samples %d, kept (depth <= own climatological MLD) %d, median Ek %.0f (IQR %.0f-%.0f)\n'
      % (len(band), len(s), EK, *np.percentile(ek, [25, 75])))

print('%-8s %11s | %-22s | %s' % ('', '', 'panel (a) P-E curve', 'panel (b)'))
print('%-8s %11s | %10s %11s | %s' %
      ('season', 'umol/m2/s', 'P/Pmax', 'integrated', '% gets through'))
for k in SEASONS:
    j = np.argmin(abs(I - INST[k]))
    print('%-8s %11.0f | %9.0f%% %10.0f%% | %13.0f%%'
          % (k, INST[k], 100 * med_inst[j], 100 * med_day[j], 100 * med_ed[j]))
