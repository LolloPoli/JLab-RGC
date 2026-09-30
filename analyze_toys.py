import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

# ------------------------------------------------------------------
# configurazione
# ------------------------------------------------------------------
period = "summer22"          # cambia secondo il tuo caso
toy_dir = "toy_model2"
outdir = "TOY_diagnostics2"
os.makedirs(outdir, exist_ok=True)

SYS_COMPONENTS = ["Acc_sys", "PID_sys", "Purity_sys", "Bin_mig_sys"]
ENABLED_COMPONENTS = {
    "Acc_sys":      True,
    "PID_sys":      True,
    "Purity_sys":   True,
    "Bin_mig_sys":  True,  
}
SYS_COMPONENTS = [c for c in SYS_COMPONENTS if ENABLED_COMPONENTS.get(c, True)]

modulations = ["AUL_sinPhi", "AUL_sin2Phi", "ALL_const", "ALL_cosPhi", "ALU_sinPhi"]
modulation_labels = {
    "AUL_sinPhi":  r"$A_{UL}^{\sin\phi}$",
    "AUL_sin2Phi": r"$A_{UL}^{\sin2\phi}$",
    "ALL_const":   r"$A_{LL}$",
    "ALL_cosPhi":  r"$A_{LL}^{\cos\phi}$",
    "ALU_sinPhi":  r"$A_{LU}^{\sin\phi}$",
}
sys_labels = {
    "Acc_sys": "Acceptance + Tracking",
    "PID_sys": "PID efficiency",
    "Purity_sys": "Purity and Contamination",
    "Bin_mig_sys": "Bin Migration",
}

comp_colors = {
    "Acc_sys": "tab:blue",
    "PID_sys": "tab:orange",
    "Purity_sys": "tab:green",
    "Bin_mig_sys": "tab:red",
}

# ------------------------------------------------------------------
# carica tutti i toy in un unico dataframe, con colonna 'toy'
# ------------------------------------------------------------------
pattern = os.path.join(toy_dir, f"table_RGC_MC_{period}_zPt_test_*.csv")
files = sorted(glob.glob(pattern))
if len(files) == 0:
    raise FileNotFoundError(f"Nessun file trovato con pattern: {pattern}")

print(f"Trovati {len(files)} file toy.")

dfs = []
for f in files:
    # estrae il numero del toy dal nome file
    toy_id = int(f.split("_")[-1].replace(".csv", ""))
    d = pd.read_csv(f, skipinitialspace=True)
    d.columns = d.columns.str.strip()
    d["modulation"] = d["modulation"].str.strip()
    d["toy"] = toy_id
    dfs.append(d)

all_toys = pd.concat(dfs, ignore_index=True)
n_toys = all_toys["toy"].nunique()
print(f"Numero di toy caricati: {n_toys}")

# ------------------------------------------------------------------
# calcolo bias (media), RMS (std), SEM e significativita' per ciascuna
# combinazione (modulation, bin_zPt, componente sistematica)
# ------------------------------------------------------------------
summary_rows = []

for modul in modulations:
    sub_mod = all_toys[all_toys["modulation"] == modul]

    for bin_zpt in sorted(sub_mod["bin_zPt"].unique()):
        sub = sub_mod[sub_mod["bin_zPt"] == bin_zpt]
        mean_z = sub["mean_z"].iloc[0]
        mean_pt = sub["mean_PhT"].iloc[0]

        row = {
            "modulation": modul,
            "bin_zPt": bin_zpt,
            "mean_z": mean_z,
            "mean_PhT": mean_pt,
            "n_toys": len(sub),
        }

        SIGNIF_THRESHOLD = 1.5   # soglia unica, usata sia nel calcolo che nella stampa

        for comp in SYS_COMPONENTS:
            vals = sub[comp].values
            bias = np.mean(vals)
            rms = np.std(vals, ddof=1) if len(vals) > 1 else np.nan
            sem = rms / np.sqrt(len(vals)) if len(vals) > 1 else np.nan
            signif = np.abs(bias) / sem if (sem and sem > 0) else np.nan

            # sistematica finale:
            # - se il bias e' significativo (>= soglia): uso il bias stesso come sistematica
            # - se NON e' significativo: non posso dire che sia zero, quindi uso il SEM
            #   come stima conservativa della sensibilita' del test (limite superiore)
            if not np.isnan(signif) and signif >= SIGNIF_THRESHOLD:
                final_sys = np.abs(bias)
            elif not np.isnan(sem):
                final_sys = sem
            else:
                final_sys = 0.0   # solo se non hai nemmeno un SEM valido (es. un solo toy)

            row[f"{comp}_bias"] = bias
            row[f"{comp}_rms"] = rms
            row[f"{comp}_sem"] = sem
            row[f"{comp}_signif"] = signif
            row[f"{comp}_raw_abs"] = np.abs(bias)
            row[f"{comp}_final"] = final_sys

        summary_rows.append(row)

summary = pd.DataFrame(summary_rows)
summary_path = os.path.join(outdir, f"toy_summary_{period}.csv")
summary.to_csv(summary_path, index=False)
print(f"Riepilogo salvato in: {summary_path}")

# ------------------------------------------------------------------
# stampa un avviso per i bin/componenti dove il bias e' significativo
# ------------------------------------------------------------------
print(f"\n=== Bin con bias sistematico significativo (>= {SIGNIF_THRESHOLD} sigma) ===")
for _, r in summary.iterrows():
    for comp in SYS_COMPONENTS:
        if r[f"{comp}_signif"] >= SIGNIF_THRESHOLD:
            print(f"  {r['modulation']:14s} bin_zPt={int(r['bin_zPt']):2d}  {comp:12s} "
                  f"bias={r[f'{comp}_bias']:+.4f}  sem={r[f'{comp}_sem']:.4f}  "
                  f"signif={r[f'{comp}_signif']:.1f}")
        else:
            print(f"  {r['modulation']:14s} bin_zPt={int(r['bin_zPt']):2d}  {comp:12s} "
                  f"signif={r[f'{comp}_signif']:.2f} (not significant → uso SEM={r[f'{comp}_sem']:.4f} come sistematica)")

# ------------------------------------------------------------------
# plot diagnostico: per ciascuna modulazione, 4 pannelli (uno per
# componente sistematica), bias vs bin_zPt con errorbar = SEM,
# banda ombreggiata = RMS (dispersione naturale tra toy)
# ------------------------------------------------------------------
for modul in modulations:
    sel = summary[summary["modulation"] == modul].sort_values("bin_zPt")

    fig, axes = plt.subplots(1, 4, figsize=(18, 4), sharex=True)

    for i, comp in enumerate(SYS_COMPONENTS):
        ax = axes[i]

        bias = -sel[f"{comp}_bias"].values
        rms = sel[f"{comp}_rms"].values
        sem = sel[f"{comp}_sem"].values
        x = sel["bin_zPt"].values

        # banda RMS (dispersione naturale toy-to-toy)
        #ax.fill_between(x, bias - rms, bias + rms, color=comp_colors[comp], alpha=0.15, label="±RMS (toy spread)")
        ax.fill_between(x, bias - sem, bias + sem, color=comp_colors[comp], alpha=0.15, label="±SEM (toy spread)")

        # bias con errorbar SEM
        ax.errorbar(x, bias, yerr=sem, fmt="o", color=comp_colors[comp], markersize=5,
                    capsize=3, linewidth=1.2, label="bias ± SEM")

        ax.axhline(0, color="black", linestyle="--", linewidth=1, alpha=0.6)

        ax.set_title(sys_labels[comp], fontsize=12)
        ax.set_xlabel("bin z-Pt")
        if i == 0:
            ax.set_ylabel(f"Δ = {modulation_labels[modul]}", fontsize=12)

        ax.tick_params(direction="in", top=True, right=True)
        ax.grid(alpha=0.25, linestyle=":")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", ncol=2, fontsize=10)

    fig.suptitle(f"Toy study ({n_toys} extractions) — {modulation_labels[modul]} — {period}", fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"toy_diag_{modul}_{period}.png"), dpi=300, bbox_inches="tight")
    plt.close()

print(f"\nPlot diagnostici salvati in: {outdir}/")

# ------------------------------------------------------------------
# BONUS: istogramma della distribuzione dei toy per un paio di bin
# scelti a mano (utile per controllare a occhio la forma: gaussiana? outlier?)
# ------------------------------------------------------------------
CHECK_BINS = [1, 9, 25]  # modifica secondo i bin che ti interessano di piu'

for modul in modulations:
    fig, axes = plt.subplots(len(CHECK_BINS), len(SYS_COMPONENTS), figsize=(16, 3.2 * len(CHECK_BINS)))

    for r, bin_zpt in enumerate(CHECK_BINS):
        sub = all_toys[(all_toys["modulation"] == modul) & (all_toys["bin_zPt"] == bin_zpt)]
        if sub.empty:
            continue
        for c, comp in enumerate(SYS_COMPONENTS):
            ax = axes[r, c] if len(CHECK_BINS) > 1 else axes[c]
            vals = sub[comp].values
            ax.hist(vals, bins=12, color="tab:blue", alpha=0.7, edgecolor="white")
            ax.axvline(0, color="black", linestyle="--", linewidth=1)
            ax.axvline(np.mean(vals), color="tab:red", linewidth=1.5, label=f"mean={np.mean(vals):.3f}")
            if r == 0:
                ax.set_title(sys_labels[comp], fontsize=11)
            if c == 0:
                ax.set_ylabel(f"bin_zPt={bin_zpt}", fontsize=10)
            ax.legend(fontsize=8)
            ax.tick_params(direction="in")

    fig.suptitle(f"Toy distribution — {modulation_labels[modul]} — {period}", fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"toy_hist_{modul}_{period}.png"), dpi=300, bbox_inches="tight")
    plt.close()

print("Istogrammi diagnostici (bin scelti) salvati.")

# ------------------------------------------------------------------
# SISTEMATICA TOTALE FINALE: somma in quadratura delle 4 componenti
# "final" (gia' passate per il test di significativita' >= 2 sigma)
# ------------------------------------------------------------------
final_cols = [f"{c}_final" for c in SYS_COMPONENTS]
sem_cols = [f"{c}_sem" for c in SYS_COMPONENTS]

summary["Total_sys_final"] = np.sqrt((summary[final_cols] ** 2).sum(axis=1))
# propagazione approssimata dell'incertezza sulla stima stessa del totale
# (somma in quadratura dei SEM delle componenti che contribuiscono davvero)
summary["Total_sys_final_err"] = np.sqrt((summary[sem_cols].fillna(0) ** 2).sum(axis=1))

summary.to_csv(summary_path, index=False)  # ri-salva con le colonne totali incluse

# ------------------------------------------------------------------
# REPORT FINALE: tabella leggibile, una riga per (modulation, bin_zPt),
# con il valore finale di ciascuna componente + il totale
# ------------------------------------------------------------------
report_lines = []
report_lines.append(f"{'='*100}")
report_lines.append(f"REPORT SISTEMATICHE FINALI — {period} — {n_toys} toy")
report_lines.append(f"{'='*100}\n")

for modul in modulations:
    sel = summary[summary["modulation"] == modul].sort_values("bin_zPt")
    report_lines.append(f"--- {modulation_labels[modul]}  ({modul}) ---")

    header_parts = [f"{'bin':>4}", f"{'mean_z':>7}"]
    for comp in SYS_COMPONENTS:
        short_name = comp.replace("_sys", "")
        header_parts.append(f"{short_name:>9}")
    header_parts.append(f"{'TOTAL':>10}")
    header_parts.append(f"{'±err':>8}")
    header = " ".join(header_parts)

    report_lines.append(header)
    report_lines.append("-" * len(header))

    for _, r in sel.iterrows():
        row_parts = [f"{int(r['bin_zPt']):>4}", f"{r['mean_z']:>7.3f}"]
        for comp in SYS_COMPONENTS:
            row_parts.append(f"{r[f'{comp}_final']:>9.4f}")
        row_parts.append(f"{r['Total_sys_final']:>10.4f}")
        row_parts.append(f"{r['Total_sys_final_err']:>8.4f}")
        report_lines.append(" ".join(row_parts))

    report_lines.append("")

report_text = "\n".join(report_lines)
report_path = os.path.join(outdir, f"final_report_{period}.txt")
with open(report_path, "w") as f:
    f.write(report_text)

print(report_text)
print(f"\nReport testuale salvato in: {report_path}")

# ------------------------------------------------------------------
# PLOT: breakdown a barre impilate -- quale componente domina in ogni bin
# ------------------------------------------------------------------

for modul in modulations:
    sel = summary[summary["modulation"] == modul].sort_values("bin_zPt")
    x = sel["bin_zPt"].values

    fig, ax = plt.subplots(figsize=(12, 5))

    bottom = np.zeros(len(sel))
    # impila i quadrati (cosi' l'altezza totale della barra e' coerente con la somma in quadratura)
    for comp in SYS_COMPONENTS:
        heights = sel[f"{comp}_final"].values #**2
        ax.bar(x, heights, bottom=bottom, color=comp_colors[comp], label=sys_labels[comp], width=0.8)
        bottom += heights

    ax.set_xlabel("bin z-Pt")
    ax.set_ylabel(r"$\sigma_{sys}$ (contributo in quadratura)")
    ax.set_title(f"Breakdown final systematics — {modulation_labels[modul]} — {period}", fontsize=13)
    ax.legend(fontsize=10)
    ax.tick_params(direction="in", top=True, right=True)
    ax.grid(alpha=0.25, linestyle=":", axis="y")

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"toy_breakdown_{modul}_{period}.png"), dpi=300, bbox_inches="tight")
    plt.close()

print(f"Plot breakdown salvati in: {outdir}/")

# ------------------------------------------------------------------
# PLOT: confronto della sistematica TOTALE tra tutte le modulazioni
# ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(12, 5))
markers = ["o", "s", "^", "D", "v"]
colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]

for modul, mk, col in zip(modulations, markers, colors):
    sel = summary[summary["modulation"] == modul].sort_values("bin_zPt")
    ax.errorbar(
        sel["bin_zPt"], sel["Total_sys_final"], yerr=sel["Total_sys_final_err"],
        fmt=mk, color=col, markersize=6, capsize=3, linewidth=1.2,
        label=modulation_labels[modul], markerfacecolor="white"
    )

ax.axhline(0, color="black", linestyle="--", linewidth=1, alpha=0.5)
ax.set_xlabel("bin z-Pt")
ax.set_ylabel(r"$\sigma_{sys}^{tot}$")
ax.set_title(f"Total systematics among the modulations — {period}", fontsize=13)
ax.legend(fontsize=10, ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.15))
ax.tick_params(direction="in", top=True, right=True)
ax.grid(alpha=0.25, linestyle=":")

plt.tight_layout()
plt.savefig(os.path.join(outdir, f"toy_total_comparison_{period}.png"), dpi=300, bbox_inches="tight")
plt.close()

print(f"Plot di confronto totale salvato in: {outdir}/toy_total_comparison_{period}.png")
print("\nFatto!")


fig, ax = plt.subplots(figsize=(12, 4))
markers = ["o", "s", "^", "D"]
colors = ["tab:blue", "tab:orange", "tab:green", "tab:red"]

for modul, mk, col in zip(modulations, markers, colors):
    sel = summary[summary["modulation"] == modul].sort_values("bin_zPt")
    ax.errorbar(
        sel["bin_zPt"], sel["Total_sys_final"], yerr=sel["Total_sys_final_err"],
        fmt=mk, color=col, markersize=6, capsize=3, linewidth=1.2,
        label=modulation_labels[modul], markerfacecolor="white"
    )

ax.axhline(0, color="black", linestyle="--", linewidth=1, alpha=0.5)
ax.set_ylim(-0.05, 0.06)
ax.set_xlabel("bin z-Pt")
ax.set_ylabel(r"$\sigma_{sys}^{tot}$")
ax.set_title(f"Total systematics among the modulations — {period}", fontsize=13)
ax.legend(fontsize=10, ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.15))
ax.tick_params(direction="in", top=True, right=True)
ax.grid(alpha=0.25, linestyle=":")

plt.tight_layout()
plt.savefig(os.path.join(outdir, f"toy_total_comparison_cut_{period}.png"), dpi=300, bbox_inches="tight")
plt.close()

print(f"Plot di confronto totale salvato in: {outdir}/toy_total_comparison_cut_{period}.png")
print("\nFatto!")