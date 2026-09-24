import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit

# =========================
# SETTINGS
# =========================

file_nh3 = "output_RGC_NH3_asymmetries_fall22.csv"
file_nd3 = "output_RGC_ND3_asymmetries_fall22.csv"

outdir = "fit_results"
plotdir = "fit_results"
os.makedirs(plotdir, exist_ok=True)
os.makedirs(outdir, exist_ok=True)

titles_asymmetries = [
    r"$0.0 < P_{hT} < 0.25 \ GeV$",
    r"$0.25 < P_{hT} < 0.5 \ GeV$",
    r"$0.5 < P_{hT} < 0.8 \ GeV$",
    r"$0.8 < P_{hT} < 1.4 \ GeV$"
]

# gruppi PhT
pt_groups = {
    "PhT_bin1": range(1,8),
    "PhT_bin2": range(8,15),
    "PhT_bin3": range(15,22),
    "PhT_bin4": range(22,26)
}

# asimmetrie
asymmetries = [
    ("AUL_sinPhi", "AUL_sinPhi_err"),
    ("AUL_sin2Phi", "AUL_sin2Phi_err"),
    ("ALL", "ALL_err"),
    ("ALL_cosPhi", "ALL_cosPhi_err"),
    ("ALU_sinPhi", "ALU_sinPhi_err")
]

labels = {
    "AUL_sinPhi": r"$A_{UL}^{\sin\phi}$",
    "AUL_sin2Phi": r"$A_{UL}^{\sin2\phi}$",
    "ALL": r"$A_{LL}$",
    "ALL_cosPhi": r"$A_{LL}^{\cos\phi}$",
    "ALU_sinPhi": r"$A_{LU}^{\sin\phi}$"
}

# parametro D-state
w = 0.05
C = 2/(1-1.5*w)

# Plot

def plot_fit_1x4(z_list, A_list, Aerr_list, popt_list, pcov_list, titles, ylabel, filename, color_data='black', color_fit='red'):
    fig, axes = plt.subplots(1, 4, figsize=(18,4), sharex=True, sharey=True)

    for i in range(4):
        ax = axes[i]
        z_fit = np.linspace(0.2, 0.9, 100)
        popt = popt_list[i]
        pcov = pcov_list[i]

        if pcov is not None and not np.isnan(popt).any():
            try:
                # Campionamento per la banda (ora su 2 parametri)
                ps = np.random.multivariate_normal(popt, pcov, 500)
                y_variants = np.array([power_law_fixed_beta(z_fit, *p) for p in ps])
                
                lower = np.percentile(y_variants, 16, axis=0)
                upper = np.percentile(y_variants, 84, axis=0)
                
                ax.fill_between(z_fit, lower, upper, color=color_fit, alpha=0.2)
                ax.plot(z_fit, power_law_fixed_beta(z_fit, *popt), color=color_fit, lw=1.5)
            except Exception as e:
                print(f"Errore banda: {e}")

        ax.errorbar(z_list[i], A_list[i], yerr=Aerr_list[i], fmt="o", capsize=3, color=color_data)
        ax.axhline(0, color="gray", ls="--", alpha=0.5)
        ax.set_ylim(-0.3, 0.3) # Stringiamo il range per vedere meglio le asimmetrie piccole
        ax.set_title(titles[i])
        ax.grid(alpha=0.2)
        if i == 0: ax.set_ylabel(ylabel)
        ax.set_xlabel("z")
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()

# =========================
# FIT MODEL
# =========================

def linear(z, a, b):
    return a*z + b

beta_fixed = 2.0
def power_law_fixed_beta(z, N, alpha):
    return N * (z**alpha)

# =========================
# FUNZIONE FIT
# =========================

def do_fit(z, A, Aerr):
    # Rimuoviamo i punti con errore zero o NaN che rompono il fit
    mask = ~np.isnan(A) & ~np.isnan(Aerr) & (Aerr > 0) & (z>0.199)
    z, A, Aerr = z[mask], A[mask], Aerr[mask]
    if len(z) < 2: return [np.nan]*2, None, [np.nan]*2, np.nan, np.nan
    try:
        # Guadiamo l'altezza media dei dati per azzeccare N
        # A 0.3 (centro del range), z^1 * (1-z)^2 vale circa 0.15
        n_guess = np.mean(A) / 0.15
        p0 = [n_guess, 1.1]
        
        # Bounds: N tra -1 e 1 (asimmetria max 100%), alpha tra 0 e 5
        bounds = ([-1, 0.01], [1, 5.0])
        
        popt, pcov = curve_fit(
            power_law_fixed_beta, z, A,
            p0=p0, sigma=Aerr, bounds=bounds,
            absolute_sigma=True, maxfev=20000
        )
        
        perr = np.sqrt(np.diag(pcov))
        chi2 = np.sum(((A - power_law_fixed_beta(z, *popt)) / Aerr)**2)
        ndf = len(z) - len(popt)

        return popt, pcov, perr, chi2, ndf

    except Exception:
        # Se fallisce il fit curvo, proviamo il fit lineare N * z
        # (Opzionale, ma aiuta a non avere grafici vuoti)
        return [np.nan]*2, None, [np.nan]*2, np.nan, np.nan
    
def power_global(z, N, alpha):
    return N * z**alpha
    
def do_fit_global(z, A, Aerr):

    mask = ~np.isnan(A) & ~np.isnan(Aerr) & (Aerr > 0)
    z, A, Aerr = z[mask], A[mask], Aerr[mask]

    if len(z) < 3:
        return [np.nan, np.nan], None, [np.nan, np.nan], np.nan, np.nan

    try:
        popt, pcov = curve_fit(
            power_global,
            z, A,
            p0=[0.1, 1.0],
            sigma=Aerr,
            bounds=([-1, 0.01], [1, 5]),
            absolute_sigma=True,
            maxfev=10000
        )

        perr = np.sqrt(np.diag(pcov))
        chi2 = np.sum(((A - power_global(z, *popt)) / Aerr)**2)
        ndf = len(z) - len(popt)

        return popt, pcov, perr, chi2, ndf

    except Exception as e:
        print("Global fit failed:", e)
        return [np.nan, np.nan], None, [np.nan, np.nan], np.nan, np.nan
    
def plot_global(z, A, Aerr, popt, pcov, title, filename, color="black"):

    plt.figure()

    plt.errorbar(z, A, yerr=Aerr, fmt="o", capsize=3, color=color, label="Data")

    if pcov is not None and not np.isnan(popt).any():
        z_fit = np.linspace(0.2, 0.9, 200)
        A_fit = power_global(z_fit, *popt)

        plt.plot(z_fit, A_fit, color="red", lw=2, label="Global fit")

        # banda errore
        try:
            ps = np.random.multivariate_normal(popt, pcov, 300)
            y_var = np.array([power_global(z_fit, *p) for p in ps])

            low = np.percentile(y_var, 16, axis=0)
            high = np.percentile(y_var, 84, axis=0)

            plt.fill_between(z_fit, low, high, color="red", alpha=0.2)
        except:
            pass

    plt.axhline(0, color="gray", ls="--")

    plt.xlabel("z")
    plt.ylabel("Asymmetry")
    plt.title(title)

    plt.grid(alpha=0.3)
    plt.legend(frameon=False)

    plt.xlim(0.2, 0.9)
    plt.ylim(-0.4, 0.4)

    plt.savefig(filename, dpi=300)
    plt.close()

# =========================
# LOAD DATA
# =========================

nh3 = pd.read_csv(file_nh3, skipinitialspace=True)
nd3 = pd.read_csv(file_nd3, skipinitialspace=True)

results = []

# =========================
# MAIN LOOP
# =========================

for asym, err in asymmetries:
    z_NH3, A_NH3, Aerr_NH3, fit_NH3, cov_NH3 = [], [], [], [], []
    z_ND3, A_ND3, Aerr_ND3, fit_ND3, cov_ND3 = [], [], [], [], []
    z_DIFF, A_DIFF, Aerr_DIFF, fit_DIFF, cov_DIFF = [], [], [], [], []
    z_N, A_N, Aerr_N, fit_N, cov_N = [], [], [], [], []
    titles = []
    
    # ===== GLOBAL DATA =====

    print(f"\n=== {asym} ===")

    for label, bins in pt_groups.items():

        nh3_sel = nh3[nh3["bin_zPt"].isin(bins)].sort_values("mean_z")
        nd3_sel = nd3[nd3["bin_zPt"].isin(bins)].sort_values("mean_z")

        z = nh3_sel["mean_z"].values

        # =====================
        # NH3
        # =====================
        Ap = nh3_sel[asym].values
        Ap_err = nh3_sel[err].values

        popt, pcov,perr, chi2, ndf = do_fit(z, Ap, Ap_err)
        z_NH3.append(z)
        A_NH3.append(Ap)
        Aerr_NH3.append(Ap_err)
        fit_NH3.append(popt)
        cov_NH3.append(pcov)

        print(f"{label} NH3: N={popt[0]:.4f}±{perr[0]:.4f}, alpha={popt[1]:.4f}±{perr[1]:.4f}, beta={beta_fixed} (fixed), chi2/ndf={chi2/ndf:.2f}")

        results.append([asym, label, "NH3", *popt, beta_fixed,*perr, 0.0,chi2, ndf])

        # =====================
        # ND3
        # =====================
        Ad = nd3_sel[asym].values
        Ad_err = nd3_sel[err].values

        popt, pcov, perr, chi2, ndf = do_fit(z, Ad, Ad_err)
        z_ND3.append(z)
        A_ND3.append(Ad)
        Aerr_ND3.append(Ad_err)
        fit_ND3.append(popt)
        cov_ND3.append(pcov)

        # Esempio per ND3 (fai lo stesso per NH3, DIFF e NEUTRON)
        print(f"{label} ND3: N={popt[0]:.4f}±{perr[0]:.4f}, alpha={popt[1]:.4f}±{perr[1]:.4f}, beta={beta_fixed} (fixed), chi2/ndf={chi2/ndf:.2f}")

        results.append([asym, label, "ND3", *popt, beta_fixed,*perr, 0.0,chi2, ndf])

        # =====================
        # DIFFERENZA (ND3 - NH3)
        # =====================
        Adiff = Ad - Ap
        Adiff_err = np.sqrt(Ad_err**2 + Ap_err**2)

        popt, pcov, perr, chi2, ndf = do_fit(z, Adiff, Adiff_err)
        z_DIFF.append(z)
        A_DIFF.append(Adiff)
        Aerr_DIFF.append(Adiff_err)
        fit_DIFF.append(popt)
        cov_DIFF.append(pcov)

        print(f"{label} DIFF: N={popt[0]:.4f}±{perr[0]:.4f}, alpha={popt[1]:.4f}±{perr[1]:.4f}, beta={beta_fixed} (fixed), chi2/ndf={chi2/ndf:.2f}")

        results.append([asym, label, "ND3-NH3", *popt, beta_fixed,*perr, 0.0,chi2, ndf])

        # =====================
        # NEUTRONE (approx)
        # =====================
        An = C*Ad - Ap
        An_err = np.sqrt((C*Ad_err)**2 + Ap_err**2)

        popt, pcov, perr, chi2, ndf = do_fit(z, An, An_err)
        z_N.append(z)
        A_N.append(An)
        Aerr_N.append(An_err)
        fit_N.append(popt)
        cov_N.append(pcov)

        print(f"{label} NEUTRON: N={popt[0]:.4f}±{perr[0]:.4f}, alpha={popt[1]:.4f}±{perr[1]:.4f}, beta={beta_fixed} (fixed), chi2/ndf={chi2/ndf:.2f}")

        results.append([asym, label, "neutron", *popt, beta_fixed,*perr, 0.0,chi2, ndf])
        
        titles.append(label)
        
    z_all = np.concatenate(z_NH3)
    Ap_all = np.concatenate(A_NH3)
    Ap_err_all = np.concatenate(Aerr_NH3)
    Ad_all = np.concatenate(A_ND3)
    Ad_err_all = np.concatenate(Aerr_ND3)
    Adiff_all = np.concatenate(A_DIFF)
    Adiff_err_all = np.concatenate(Aerr_DIFF)
    An_all = np.concatenate(A_N)
    An_err_all = np.concatenate(Aerr_N)

    # plot
    # ===== NH3 =====
    popt, pcov, perr, chi2, ndf = do_fit_global(z_all, Ap_all, Ap_err_all)
    print(f"GLOBAL NH3: N={popt[0]:.4f}±{perr[0]:.4f}, alpha={popt[1]:.4f}±{perr[1]:.4f}, chi2/ndf={chi2/ndf:.2f}")

    plot_global(
        z_all, Ap_all, Ap_err_all,
        popt, pcov,
        f"{labels[asym]} GLOBAL (NH3)",
        f"{plotdir}/{asym}_GLOBAL_NH3.png",
        color="cornflowerblue"
    )
    plot_fit_1x4(
        z_NH3, A_NH3, Aerr_NH3, fit_NH3, cov_NH3,
        titles,
        f"{labels[asym]} (NH3)",
        f"{plotdir}/{asym}_NH3_1x4.png",
        color_data='cornflowerblue', color_fit='navy'
    )
    
    # ===== ND3 =====
    popt, pcov, perr, chi2, ndf = do_fit_global(z_all, Ad_all, Ad_err_all)

    print(f"GLOBAL ND3: N={popt[0]:.4f}±{perr[0]:.4f}, alpha={popt[1]:.4f}±{perr[1]:.4f}, chi2/ndf={chi2/ndf:.2f}")

    plot_global(
        z_all, Ad_all, Ad_err_all,
        popt, pcov,
        f"{labels[asym]} GLOBAL (ND3)",
        f"{plotdir}/{asym}_GLOBAL_ND3.png",
        color="orange"
    )
    plot_fit_1x4(
        z_ND3, A_ND3, Aerr_ND3, fit_ND3, cov_ND3,
        titles,
        f"{labels[asym]} (ND3)",
        f"{plotdir}/{asym}_ND3_1x4.png",
        color_data='orange', color_fit='darkred'
    )
    
    # ===== DIFF =====
    popt, pcov, perr, chi2, ndf = do_fit_global(z_all, Adiff_all, Adiff_err_all)

    print(f"GLOBAL DIFF: N={popt[0]:.4f}±{perr[0]:.4f}, alpha={popt[1]:.4f}±{perr[1]:.4f}, chi2/ndf={chi2/ndf:.2f}")

    plot_global(
        z_all, Adiff_all, Adiff_err_all,
        popt, pcov,
        f"{labels[asym]} GLOBAL (ND3-NH3)",
        f"{plotdir}/{asym}_GLOBAL_DIFF.png",
        color="crimson"
    )
    plot_fit_1x4(
        z_DIFF, A_DIFF, Aerr_DIFF, fit_DIFF, cov_DIFF,
        titles,
        f"{labels[asym]} (ND3 - NH3)",
        f"{plotdir}/{asym}_DIFF_1x4.png",
        color_data='crimson', color_fit='darkgreen'
    )
    
    # ===== NEUTRON =====
    popt, pcov, perr, chi2, ndf = do_fit_global(z_all, An_all, An_err_all)

    print(f"GLOBAL NEUTRON: N={popt[0]:.4f}±{perr[0]:.4f}, alpha={popt[1]:.4f}±{perr[1]:.4f}, chi2/ndf={chi2/ndf:.2f}")

    plot_global(
        z_all, An_all, An_err_all,
        popt, pcov,
        f"{labels[asym]} GLOBAL (neutron)",
        f"{plotdir}/{asym}_GLOBAL_neutron.png",
        color="seagreen"
    )
    plot_fit_1x4(
        z_N, A_N, Aerr_N, fit_N, cov_N,
        titles,
        f"{labels[asym]} (neutron approx)",
        f"{plotdir}/{asym}_neutron_1x4.png",
        color_data='seagreen', color_fit='indigo'
    )
    
    
# =========================
# SAVE RESULTS
# =========================

df_results = pd.DataFrame(results, columns=[
    "asymmetry", "PhT_bin", "target",
    "N", "alpha", "beta", 
    "err_N", "err_alpha", "err_beta",
    "chi2", "ndf"
])

outfile = f"{outdir}/fit_results.csv"
df_results.to_csv(outfile, index=False)

print(f"\nResults saved in {outfile}")