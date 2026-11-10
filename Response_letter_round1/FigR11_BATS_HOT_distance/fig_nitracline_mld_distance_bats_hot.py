# Used for: round-1 letter Fig. R11
# Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\session_2026-09-28_scripts\fig_nitracline_mld_distance_bats_hot.py (data paths updated to I:\Prj; computation unchanged)
"""Observed nitracline-to-mixed-layer distance at BATS and HOT, 2003 to 2022.

WHY     Reviewer 3's reply holds the nitracline at climatology. Nicolas asks
        (comments #269 and #273) how much the nitracline actually varies within
        and between seasons, and whether time-series or BioArgo data can test it.
        This is the observational counterpart of Figure R9: the same quantity
        (nitracline minus MLD, positive when the nitracline lies below the
        mixed layer) but measured, at the two subtropical time-series stations
        that fall inside the 20 to 40 degree band of the reply.

DATA    Both quantities come from the SAME bottle cast, so the distance is a
        single observed number with no model in it.
          BATS  3782_v8_bats_bottle.csv   (31.7N, 64.2W)
          HOT   3773_v3_niskin_hot001...  (22.75N, 158W, station ALOHA only)
        Nitracline: first crossing of NO3 = 3 umol/kg, linearly interpolated,
        the same criterion and method as get_nitracline.m in the study.
        MLD: first depth where sigma_theta exceeds its 10 m value by
        0.03 kg/m3, linearly interpolated, the same threshold as MLD one_d03.

QC      A cast is used only if it resolves both quantities: a density sample at
        or above 15 m, a sub-threshold nitrate sample in the upper 40 m, at
        least 6 valid levels, and the crossing bracketed by real samples no
        more than 120 m apart. Casts that start below the nitracline (common in
        the raw files) would otherwise report the depth of their first bottle.

CAVEAT  HOT bottle spacing in the upper 250 m is 17 m (BATS is 4 m), so the HOT
        summer MLD of roughly 30 to 40 m carries an uncertainty of order 10 m.
"""
import os
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                                   # code for revision
OUT = os.path.join(ROOT, 'mld_npp_relationship_grid')
INSITU = (r"E:\OneDrive - Indiana State University\Research\Projects"
          r"\3D bio model\Code\npp_insitu")

NO3_THR = 3.0          # umol/kg, as in the study
SIG_THR = 0.03         # kg/m3, as in MLD one_d03
COL = {'Spring': '#C87F00', 'Summer': '#009E73',
       'Autumn': '#4D4D4D', 'Winter': '#0072B2'}
ORDER = ['Winter', 'Spring', 'Summer', 'Autumn']
SEAS = {12: 'Winter', 1: 'Winter', 2: 'Winter', 3: 'Spring', 4: 'Spring',
        5: 'Spring', 6: 'Summer', 7: 'Summer', 8: 'Summer', 9: 'Autumn',
        10: 'Autumn', 11: 'Autumn'}


def _clean(depth, val):
    o = np.argsort(depth)
    d, x = np.asarray(depth, float)[o], np.asarray(val, float)[o]
    k = np.isfinite(d) & np.isfinite(x)
    return d[k], x[k]


def nitracline(depth, no3):
    """First crossing of NO3 = 3 umol/kg, as get_nitracline.m but QC'd."""
    d, x = _clean(depth, no3)
    if len(d) < 6:
        return np.nan
    shallow = x[d <= 40]
    if shallow.size == 0 or shallow.min() >= NO3_THR or x.max() < NO3_THR:
        return np.nan
    i = int(np.argmax(x >= NO3_THR))
    if i == 0 or d[i] - d[i - 1] > 120:
        return np.nan
    return d[i] if x[i] == x[i - 1] else np.interp(NO3_THR, [x[i - 1], x[i]],
                                                   [d[i - 1], d[i]])


def mld(depth, sigma):
    """Depth where sigma_theta exceeds its 10 m value by 0.03 kg/m3."""
    d, s = _clean(depth, sigma)
    if len(d) < 6 or d.min() > 15:
        return np.nan
    s_ref = np.interp(10.0, d, s)
    over = s > s_ref + SIG_THR
    if not over.any():
        return np.nan
    i = int(np.argmax(over))
    if i == 0:
        return d[0]
    return np.interp(s_ref + SIG_THR, [s[i - 1], s[i]], [d[i - 1], d[i]])


def load_bats():
    f = pd.read_csv(os.path.join(INSITU, '3782_v8_bats_bottle.csv'),
                    low_memory=False)
    f['t'] = pd.to_datetime(f.ISO_DateTime_UTC, format='%Y-%m-%dT%H:%MZ',
                            errors='coerce')
    f = f[(f.t.dt.year >= 2003) & (f.t.dt.year <= 2022)]
    for c in ['Depth', 'Sigma_theta', 'NO3_plus_NO2']:
        f[c] = pd.to_numeric(f[c], errors='coerce').replace(-999, np.nan)
    rec = []
    for _, g in f.groupby(['Cruise_num', 'Cast_num']):
        rec.append((g.t.iloc[0],
                    nitracline(g.Depth.values, g.NO3_plus_NO2.values),
                    mld(g.Depth.values, g.Sigma_theta.values)))
    return pd.DataFrame(rec, columns=['t', 'nit', 'mld']).dropna()


def load_hot():
    f = pd.read_csv(os.path.join(INSITU,
                    '3773_v3_niskin_hot001_yr01_to_hot348_yr35.csv'),
                    low_memory=False)
    f['t'] = pd.to_datetime(f.Sampling_ISO_DateTime_UTC,
                            format='%Y-%m-%dT%H:%MZ', errors='coerce')
    f = f[(f.t.dt.year >= 2003) & (f.t.dt.year <= 2022)]
    for c in ['CTDPRS', 'SIGMA', 'NO2_NO3']:
        f[c] = pd.to_numeric(f[c], errors='coerce')
    f.loc[f.NO2_NO3 > 900, 'NO2_NO3'] = np.nan       # 9999 fill values
    f.loc[f.SIGMA > 100, 'SIGMA'] = np.nan
    lat = pd.to_numeric(f.Latitude, errors='coerce')
    f = f[(lat > 22.4) & (lat < 23.1)]               # station ALOHA only
    rec = []
    for _, g in f.groupby(['Cruise_number', 'STNNBR', 'CASTNO']):
        rec.append((g.t.iloc[0],
                    nitracline(g.CTDPRS.values, g.NO2_NO3.values),
                    mld(g.CTDPRS.values, g.SIGMA.values)))
    return pd.DataFrame(rec, columns=['t', 'nit', 'mld']).dropna()


def decyear(t):
    return t.dt.year + (t.dt.dayofyear - 1) / 365.25


def panel(ax, D, title, label):
    """One line per season: the seasonal-mean distance in each year."""
    D = D.copy()
    D['sep'] = D.nit - D.mld
    D['season'] = D.t.dt.month.map(SEAS)
    D['yr'] = D.t.dt.year
    for k, s in enumerate(ORDER):
        g = D[D.season == s].groupby('yr').sep.mean()
        ax.plot(g.index, g.values, '-o', color=COL[s], lw=2.0, ms=5.0,
                mfc=COL[s], mec='none', label=s, zorder=3)
        # climatology +- interannual SD, and the Theil-Sen trend with its MK p
        x = g.index.values.astype(float)
        ts = stats.theilslopes(g.values, x, .95)
        _, p = stats.kendalltau(x, g.values)
        ax.text(.013, .975 - .075 * k,
                '%-7s %3.0f $\\pm$ %2.0f m   %+.2f m yr$^{-1}$ (p = %.2f)'
                % (s, g.mean(), g.std(), ts[0], p),
                transform=ax.transAxes, color=COL[s], fontsize=15.5,
                family='DejaVu Sans Mono', va='top', ha='left', zorder=6)
    ax.set_xlim(2002.4, 2022.6)
    ax.set_ylim(0, 430)
    ax.set_ylabel('Nitracline $-$ MLD (m)', fontsize=19)
    ax.set_title(title, fontsize=20)
    # panel letter at the upper-left, outside the axes, as in Figure R10
    ax.text(-.085, 1.06, label, transform=ax.transAxes, fontsize=21,
            va='bottom', ha='left')
    ax.tick_params(labelsize=17)
    ax.grid(alpha=.25)
    return D


def report(name, D):
    print('\n=== %s ===  n = %d casts with both quantities resolved' % (name, len(D)))
    print('  distance: mean %.0f m, median %.0f m, min %.0f m, 5th pct %.0f m'
          % (D.sep.mean(), D.sep.median(), D.sep.min(), np.percentile(D.sep, 5)))
    print('  casts with distance <= 0 : %d' % (D.sep <= 0).sum())
    print('  %-8s %5s %8s %9s %9s %8s' % ('season', 'n', 'mean sep', 'TheilSen', 'MK p', '20 yr'))
    for s in ORDER:
        g = D[D.season == s].groupby(D.t.dt.year).sep.mean()
        if len(g) < 8:
            continue
        x = g.index.values.astype(float)
        ts = stats.theilslopes(g.values, x, .95)
        tau, p = stats.kendalltau(x, g.values)
        print('  %-8s %5d %8.0f %+9.2f %9.3f %+8.0f' % (s, (D.season == s).sum(),
              D[D.season == s].sep.mean(), ts[0], p, ts[0] * 19))
    ann = D.groupby(D.t.dt.year).sep.mean()
    x = ann.index.values.astype(float)
    ts = stats.theilslopes(ann.values, x, .95)
    tau, p = stats.kendalltau(x, ann.values)
    print('  ANNUAL   %5d %8.0f %+9.2f %9.3f %+8.0f   [95%% CI %+.2f to %+.2f]'
          % (len(D), D.sep.mean(), ts[0], p, ts[0] * 19, ts[2], ts[3]))


if __name__ == '__main__':
    B, H = load_bats(), load_hot()
    plt.rcParams.update({'font.size': 17})
    fig, ax = plt.subplots(2, 1, figsize=(13.0, 10.0), sharex=True)
    DB = panel(ax[0], B, 'BATS, 31.7$^\\circ$N 64.2$^\\circ$W', '(a)')
    DH = panel(ax[1], H, 'HOT (station ALOHA), 22.75$^\\circ$N 158$^\\circ$W', '(b)')
    ax[1].set_xlabel('Year', fontsize=19)
    ax[1].set_xticks(np.arange(2004, 2023, 3))
    ax[1].set_xticklabels([str(y) for y in np.arange(2004, 2023, 3)])
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, ncol=4, frameon=False, fontsize=19,
               loc='lower center', bbox_to_anchor=(.5, .008))
    fig.tight_layout(rect=[0, .062, 1, 1])
    png = os.path.join(OUT, 'fig_nitracline_mld_distance_bats_hot.png')
    fig.savefig(png, dpi=200)
    print('saved %s' % png)
    report('BATS', DB)
    report('HOT', DH)

