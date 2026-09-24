import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import os

# cartella output
outdir = "ASYMMETRIES_plot_w_sys"
os.makedirs(outdir, exist_ok=True)

# ------------------------------------------------------------------
# carica csv delle asimmetrie
# ------------------------------------------------------------------
sum22 = pd.read_csv("output_RGC_NH3_asymmetries_sum22.csv", skipinitialspace=True)
fall22 = pd.read_csv("output_RGC_NH3_asymmetries_fall22.csv", skipinitialspace=True)
# spring23 = pd.read_csv("output_RGC_NH3_asymmetries_spring23.csv", skipinitialspace=True)

sum22.columns = sum22.columns.str.strip()
fall22.columns = fall22.columns.str.strip()

# ------------------------------------------------------------------
# carica csv delle sistematiche, formato toy_summary (gia' prodotto da
# analyze_toys.py: contiene Total_sys_final, gia' testato per significativita'
# sui toy e sommato in quadratura -- qui non si ricalcola piu' nulla)
# ------------------------------------------------------------------
sys_sum22 = pd.read_csv("TOY_diagnostics3/toy_summary_summer22.csv", skipinitialspace=True)
sys_fall22 = pd.read_csv("TOY_diagnostics3/toy_summary_fall22.csv", skipinitialspace=True)

sys_sum22.columns = sys_sum22.columns.str.strip()
sys_fall22.columns = sys_fall22.columns.str.strip()
sys_sum22["modulation"] = sys_sum22["modulation"].str.strip()
sys_fall22["modulation"] = sys_fall22["modulation"].str.strip()

# ------------------------------------------------------------------
# mappatura nome asimmetria (csv asimmetrie) -> nome modulation (csv sistematiche)
# ------------------------------------------------------------------
asym_to_modulation = {
    "AUL_sinPhi":  "AUL_sinPhi",
    "AUL_sin2Phi": "AUL_sin2Phi",
    "ALL":         "ALL_const",
    "ALL_cosPhi":  "ALL_cosPhi",
    "ALU_sinPhi":  "ALU_sinPhi",
}


def merge_sys_for_asym(asym_df, sys_df, asym_name):
    """
    Filtra il csv sistematiche (toy_summary) sulla modulation corrispondente
    all'asimmetria richiesta, poi fa il merge su bin_zPt.
    """
    modulation = asym_to_modulation[asym_name]
    sel = sys_df[sys_df["modulation"] == modulation][["bin_zPt", "Total_sys_final", "Total_sys_final_err"]]
    sel = sel.rename(columns={
        "Total_sys_final": f"{asym_name}_totsys",
        "Total_sys_final_err": f"{asym_name}_totsys_err",
    })
    return asym_df.merge(sel, on="bin_zPt", how="left")


# ------------------------------------------------------------------
# gruppi di PhT
# ------------------------------------------------------------------
pt_groups = {
    "PhT_bin1": range(1, 8),
    "PhT_bin2": range(8, 15),
    "PhT_bin3": range(15, 22),
    "PhT_bin4": range(22, 26)
}

# lista asimmetrie
asymmetries = [
    ("AUL_sinPhi", "AUL_sinPhi_err"),
    ("AUL_sin2Phi", "AUL_sin2Phi_err"),
    ("ALL", "ALL_err"),
    ("ALL_cosPhi", "ALL_cosPhi_err"),
    ("ALU_sinPhi", "ALU_sinPhi_err")
]

titles_asymmetries = [
    r"$0.0 < P_{hT} < 0.25 \ GeV$",
    r"$0.25 < P_{hT} < 0.5 \ GeV$",
    r"$0.5 < P_{hT} < 0.8 \ GeV$",
    r"$0.8 < P_{hT} < 1.4 \ GeV$"
]

asym_labels = {
    "AUL_sinPhi": r"$F_{UL}^{\sin\phi}/F_{UU}$",
    "AUL_sin2Phi": r"$F_{UL}^{\sin2\phi}/F_{UU}$",
    "ALL": r"$F_{LL}/F_{UU}$",
    "ALL_cosPhi": r"$F_{LL}^{\cos\phi}/F_{UU}$",
    "ALU_sinPhi": r"$F_{LU}^{\sin\phi}/F_{UU}$"
}

SYS_BOX_WIDTH = 0.035


def plot_campaign_with_sys(df, sys_df, campaign_label, color, marker, fname_suffix):
    for asym, err in asymmetries:

        sel_full = merge_sys_for_asym(df, sys_df, asym)
        totsys_col = f"{asym}_totsys"

        fig, axes = plt.subplots(1, 4, figsize=(18, 3.5), sharey=True, sharex=True)

        for i, (label, bins) in enumerate(pt_groups.items()):

            sel = sel_full[sel_full["bin_zPt"].isin(bins)].sort_values("mean_z")
            ax = axes[i]

            ax.bar(
                sel["mean_z"],
                height=2 * sel[totsys_col],
                bottom=sel[asym] - sel[totsys_col],
                width=SYS_BOX_WIDTH,
                color="r",
                alpha=0.4,
                edgecolor="none",
                zorder=1
            )

            ax.errorbar(
                sel["mean_z"],
                sel[asym],
                yerr=sel[err],
                fmt=marker,
                markersize=6,
                capsize=3,
                linewidth=1.2,
                label=f"{campaign_label} asym.",
                color=color,
                markerfacecolor="white",
                zorder=2
            )

            ax.axhline(0, color="black", linestyle="--", linewidth=1, alpha=0.6)

            if asym == "ALU_sinPhi":
                ax.set_ylim(-0.15, 0.15)
            elif asym == "ALL":
                ax.set_ylim(-0.1, 0.8)
            elif asym == "ALL_cosPhi":
                ax.set_ylim(-0.5, 0.5)
            else:
                ax.set_ylim(-0.3, 0.3)

            ax.set_xlim(0.2, 1)

            ax.set_title(titles_asymmetries[i], fontsize=12)
            ax.set_xlabel(r"$z$")

            if i == 0:
                ax.set_ylabel(asym_labels[asym], fontsize=13)

            ax.tick_params(direction="in", top=True, right=True, length=5)
            for spine in ax.spines.values():
                spine.set_linewidth(1.2)
            ax.grid(alpha=0.25, linestyle=":")

        fig.text(0.5, 1.02, rf"$K^+$ NH$_3$ Run Group C — {campaign_label}", ha="center", fontsize=14)

        handles, labels = axes[0].get_legend_handles_labels()

        handles.append(Patch(facecolor="r", alpha=0.4, edgecolor="none"))
        labels.append("Systematic")

        legend_x = 0.155 if campaign_label == "Summer22" else 0.135
        fig.legend(handles, labels, bbox_to_anchor=(legend_x, 0.335))
        
        

        plt.tight_layout()
        plt.savefig(f"{outdir}/{asym}_{fname_suffix}.png", dpi=300, bbox_inches="tight")
        plt.close()


plot_campaign_with_sys(sum22, sys_sum22, "Summer22", "tab:blue", "o", "sum22_with_sys")
plot_campaign_with_sys(fall22, sys_fall22, "Fall22", "tab:orange", "s", "fall22_with_sys")

print("Plot sum22 e fall22 (con sistematiche, separati) creati!")