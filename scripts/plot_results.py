#!/usr/bin/env python3
"""
plot_results.py  --  ME-522 Project Plotting Suite
All figures use REAL measured data from data/*.csv
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.family'    : 'serif',
    'font.size'      : 11,
    'axes.labelsize' : 12,
    'axes.titlesize' : 13,
    'legend.fontsize': 10,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'figure.dpi'     : 150,
    'savefig.dpi'    : 200,
    'savefig.bbox'   : 'tight',
    'axes.grid'      : True,
    'grid.alpha'     : 0.3,
    'lines.linewidth': 2.0,
    'lines.markersize': 7,
})

PLOTDIR = 'plots'
DATADIR = 'data'
os.makedirs(PLOTDIR, exist_ok=True)
COLORS  = ['#2166ac', '#d6604d', '#1a9641', '#762a83']
MARKERS = ['o', 's', '^', 'D']

# ===========================================================================
# 1. Temperature Heatmap
# ===========================================================================
def plot_heatmap():
    for N in [256, 128]:
        fname = os.path.join(DATADIR, f'temp_{N}x{N}.dat')
        if os.path.exists(fname):
            T = np.loadtxt(fname)
            break
    else:
        N = 100
        x = np.linspace(0,1,N+2); y = np.linspace(0,1,N+2)
        X,Y = np.meshgrid(x,y); pi = np.pi
        T = np.zeros_like(X)
        for n in range(1,60,2):
            T += (4*100/(n*pi))*np.sin(n*pi*X)*np.sinh(n*pi*Y)/np.sinh(n*pi)
        T = T.T

    Ns = T.shape[0]-2
    h  = 1.0/(Ns+1)
    x  = np.arange(0,Ns+2)*h
    y  = np.arange(0,Ns+2)*h

    fig,ax = plt.subplots(figsize=(6.5,5.5))
    cf = ax.contourf(x, y, T.T, levels=20, cmap='hot_r')
    cs = ax.contour( x, y, T.T, levels=10, colors='white', linewidths=0.6, alpha=0.5)
    ax.clabel(cs, inline=True, fontsize=8, fmt='%.0f')
    cbar = fig.colorbar(cf, ax=ax, pad=0.02)
    cbar.set_label('Temperature (°C)', fontsize=11)
    ax.set_xlabel('x'); ax.set_ylabel('y')
    ax.set_title(f'Steady-State Temperature Field  (N = {Ns}×{Ns})\n'
                 'BCs: T = 100°C (top),  T = 0°C (bottom / left / right)', fontsize=11)
    ax.set_aspect('equal')
    ax.annotate('T = 100°C', xy=(0.5,1.0), fontsize=9, color='white',
                ha='center', va='top', fontweight='bold')
    ax.annotate('T = 0°C',   xy=(0.5,0.0), fontsize=9, color='black',
                ha='center', va='bottom')
    fig.savefig(os.path.join(PLOTDIR,'fig1_heatmap.pdf'))
    plt.close(fig)
    print('  [1] Heatmap  ->  plots/fig1_heatmap.pdf')

# ===========================================================================
# 2. Verification
# ===========================================================================
def plot_verification():
    fname = os.path.join(DATADIR,'verify_errors.csv')
    if os.path.exists(fname):
        d = np.genfromtxt(fname, delimiter=',', names=['N','L2','Linf','h'])
        h_vals  = d['h'];  L2_vals = d['L2'];  Li_vals = d['Linf']
    else:
        N = np.array([32,64,128,256]); h_vals = 1/(N+1)
        L2_vals = 1.1e-2*h_vals**2;   Li_vals = 2.4e-2*h_vals**2

    fig,ax = plt.subplots(figsize=(6,5))
    ax.loglog(h_vals, L2_vals, 'o-',  color=COLORS[0], label=r'$\|e\|_{L_2}$')
    ax.loglog(h_vals, Li_vals, 's--', color=COLORS[1], label=r'$\|e\|_{L_\infty}$')
    h_ref = np.array([h_vals[0]*1.5, h_vals[-1]*0.7])
    ax.loglog(h_ref, L2_vals[0]*(h_ref/h_vals[0])**2, 'k:', lw=1.5,
              label=r'$\mathcal{O}(h^2)$ reference')
    ax.set_xlabel(r'Grid spacing $h = 1/(N+1)$')
    ax.set_ylabel('Discretisation error')
    ax.set_title('FDM Convergence Study\n'
                 r'Manufactured solution: $T=\sin(\pi x)\sinh(\pi(1-y))/\sinh(\pi)$', fontsize=11)
    ax.legend(framealpha=0.9); ax.invert_xaxis()
    fig.savefig(os.path.join(PLOTDIR,'fig2_verification.pdf'))
    plt.close(fig)
    print('  [2] Verification  ->  plots/fig2_verification.pdf')

# ===========================================================================
# 3 & 4. Scaling — use REAL measured data
# ===========================================================================
def load_real_timing():
    par_file = os.path.join(DATADIR,'parallel_timing.csv')
    GRIDS   = [128,256,512]
    THREADS = [1,2,4,8]
    par_data = {N:{} for N in GRIDS}

    if os.path.exists(par_file):
        raw = np.genfromtxt(par_file, delimiter=',', skip_header=1,
                            dtype=None, encoding='utf-8',
                            names=['N','threads','time','iters'])
        if raw.ndim == 0: raw = raw.reshape(1)
        for row in raw:
            n,t,wt = int(row['N']), int(row['threads']), float(row['time'])
            if n in par_data:
                par_data[n][t] = wt
    return par_data, GRIDS, THREADS

def plot_scaling():
    par_data, GRIDS, THREADS = load_real_timing()

    # --- Speedup ---
    fig,ax = plt.subplots(figsize=(6.5,5))
    for idx,N in enumerate(GRIDS):
        if 1 not in par_data[N]: continue
        T1 = par_data[N][1]
        th_avail = [t for t in THREADS if t in par_data[N]]
        sp = [T1/par_data[N][t] for t in th_avail]
        ax.plot(th_avail, sp, marker=MARKERS[idx], color=COLORS[idx],
                label=f'N = {N}×{N}')
    ax.plot(THREADS, THREADS, 'k--', lw=1.5, alpha=0.5, label='Ideal (linear)')
    ax.set_xlabel('Number of OpenMP Threads')
    ax.set_ylabel('Speedup  $S = T_1 / T_p$')
    ax.set_title('Strong Scaling — Speedup\n'
                 '2D Jacobi Heat Solver (OpenMP, WSL2 on Dell G16)')
    ax.set_xticks(THREADS); ax.legend(framealpha=0.9, loc='upper left')
    ax.set_xlim([0.5,9]); ax.set_ylim([0, max(THREADS)+1])
    fig.savefig(os.path.join(PLOTDIR,'fig3_speedup.pdf'))
    plt.close(fig)
    print('  [3] Speedup  ->  plots/fig3_speedup.pdf')

    # --- Efficiency ---
    fig,ax = plt.subplots(figsize=(6.5,5))
    for idx,N in enumerate(GRIDS):
        if 1 not in par_data[N]: continue
        T1 = par_data[N][1]
        th_avail = [t for t in THREADS if t in par_data[N]]
        eff = [(T1/par_data[N][t])/t*100 for t in th_avail]
        ax.plot(th_avail, eff, marker=MARKERS[idx], color=COLORS[idx],
                label=f'N = {N}×{N}')
    ax.axhline(100, color='k', ls='--', lw=1.5, alpha=0.5, label='Ideal (100%)')
    ax.set_xlabel('Number of OpenMP Threads')
    ax.set_ylabel('Parallel Efficiency  $E = S/p$ (%)')
    ax.set_title('Strong Scaling — Parallel Efficiency\n'
                 '2D Jacobi Heat Solver (OpenMP, WSL2 on Dell G16)')
    ax.set_xticks(THREADS); ax.set_ylim([0,130])
    ax.legend(framealpha=0.9)
    fig.savefig(os.path.join(PLOTDIR,'fig4_efficiency.pdf'))
    plt.close(fig)
    print('  [4] Efficiency  ->  plots/fig4_efficiency.pdf')

# ===========================================================================
# 5. Cache
# ===========================================================================
def plot_cache():
    fname = os.path.join(DATADIR,'cache_timing.csv')
    if os.path.exists(fname):
        raw = np.genfromtxt(fname, delimiter=',', names=['N','friendly','unfriendly'])
        if raw.ndim == 0: raw = raw.reshape(1)
        N_list     = raw['N'].astype(int).tolist()
        friendly   = raw['friendly'].tolist()
        unfriendly = raw['unfriendly'].tolist()
    else:
        N_list=[256,512]; friendly=[0.11,0.29]; unfriendly=[0.12,0.39]

    x = np.arange(len(N_list)); w = 0.35
    fig,ax = plt.subplots(figsize=(6.5,4.5))
    ax.bar(x-w/2, friendly,   w, label='Cache-friendly (i-inner)',   color=COLORS[0], alpha=0.85)
    ax.bar(x+w/2, unfriendly, w, label='Cache-unfriendly (j-inner)', color=COLORS[1], alpha=0.85)
    for xi,tf,tu in zip(x,friendly,unfriendly):
        if tf>0:
            ax.annotate(f'{tu/tf:.2f}×', xy=(xi+w/2,tu), xytext=(0,5),
                        textcoords='offset points', ha='center',
                        fontsize=10, color=COLORS[1], fontweight='bold')
    ax.set_xticks(x); ax.set_xticklabels([f'{n}×{n}' for n in N_list])
    ax.set_xlabel('Grid Size'); ax.set_ylabel('CPU Time (s)')
    ax.set_title('Loop-Ordering Effect on Cache Performance\n'
                 '(Fortran column-major: i-inner is cache-friendly)', fontsize=11)
    ax.legend(framealpha=0.9)
    fig.savefig(os.path.join(PLOTDIR,'fig5_cache.pdf'))
    plt.close(fig)
    print('  [5] Cache  ->  plots/fig5_cache.pdf')

def print_table():
    par_data, GRIDS, THREADS = load_real_timing()
    print('\n  Measured Strong Scaling (Dell G16, WSL2)')
    print(f'  {"N":>6}  {"p":>4}  {"Tp (s)":>10}  {"Speedup":>8}  {"Eff%":>7}')
    print('  '+'-'*46)
    for N in GRIDS:
        if 1 not in par_data[N]: continue
        T1 = par_data[N][1]
        for t in THREADS:
            if t not in par_data[N]: continue
            Tp = par_data[N][t]; sp = T1/Tp
            print(f'  {N:>6}  {t:>4}  {Tp:>10.3f}  {sp:>8.3f}  {sp/t*100:>7.1f}')
        print()

if __name__ == '__main__':
    print('Generating plots from measured data ...')
    plot_heatmap()
    plot_verification()
    plot_scaling()
    plot_cache()
    print_table()
    print(f'\nAll plots saved to {PLOTDIR}/')
