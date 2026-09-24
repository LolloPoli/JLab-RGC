import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

outdir = "plots_neutron_asymmetries"
os.makedirs(outdir, exist_ok=True)

# lettura file
nh3 = pd.read_csv("output_RGC_NH3_asymmetries_fall22.csv", skipinitialspace=True)
nd3 = pd.read_csv("output_RGC_ND3_asymmetries_fall22.csv", skipinitialspace=True)

# parametro D-state
w = 0.05
C = 2/(1-1.5*w)

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

titles_asymmetries = [
    r"$0.0 < P_{hT} < 0.25 \ GeV$",
    r"$0.25 < P_{hT} < 0.5 \ GeV$",
    r"$0.5 < P_{hT} < 0.8 \ GeV$",
    r"$0.8 < P_{hT} < 1.4 \ GeV$"
]

labels = {
    "AUL_sinPhi": r"$F_{UL}^{\sin\phi}/F_{UU}$",
    "AUL_sin2Phi": r"$F_{UL}^{\sin2\phi}/F_{UU}$",
    "ALL": r"$F_{LL}/F_{UU}$",
    "ALL_cosPhi": r"$F_{LL}^{\cos\phi}/F_{UU}$",
    "ALU_sinPhi": r"$F_{LU}^{\sin\phi}/F_{UU}$"
}

for asym, err in asymmetries:

    fig, axes = plt.subplots(1,4, figsize=(18,4), sharex=True, sharey=True)

    for i, (label, bins) in enumerate(pt_groups.items()):

        nh3_sel = nh3[nh3["bin_zPt"].isin(bins)].sort_values("mean_z")
        nd3_sel = nd3[nd3["bin_zPt"].isin(bins)].sort_values("mean_z")

        Ap = nh3_sel[asym].values
        Ap_err = nh3_sel[err].values

        Ad = nd3_sel[asym].values
        Ad_err = nd3_sel[err].values

        # neutron asymmetry
        An = C*Ad - Ap

        # errore
        An_err = np.sqrt((C*Ad_err)**2 + Ap_err**2)

        z = (nh3_sel["mean_z"] + nd3_sel["mean_z"]) / 2
        
        # differenza diretta
        Adiff = Ad - Ap
        Adiff_err = np.sqrt(Ad_err**2 + Ap_err**2)

        ax = axes[i]

        ax.errorbar(
            z,
            An,
            yerr=An_err,
            fmt="o",
            color="seagreen",
            capsize=3,
            markersize=6,
            label="Neutron (approx.)"
        )
        ax.errorbar(
            z,
            Adiff,
            yerr=Adiff_err,
            fmt="s",
            color="crimson",
            capsize=3,
            markersize=6,
            label="ND3 - NH3"
        )

        ax.axhline(0, color="black", linestyle="--")

        if asym == "ALU_sinPhi":
            ax.set_ylim(-0.2,0.2)
        elif asym == "ALL":
            ax.set_ylim(-0.8,0.2)
        else:
            ax.set_ylim(-0.5,0.5)

        ax.set_xlim(0.2,1)

        ax.set_title(titles_asymmetries[i])

        ax.tick_params(direction='in', top=True, right=True)

        for spine in ax.spines.values():
            spine.set_linewidth(1.2)

        ax.grid(alpha=0.3)

        if i == 0:
            ax.set_ylabel(labels[asym] + r" (neutron)")
            ax.legend(loc="upper left")

        ax.set_xlabel("z")

    fig.suptitle(f"{labels[asym]}  |  $K^+$  (neutron approx. vs ND3 - NH3)", fontsize=14)

    plt.tight_layout()

    plt.savefig(f"{outdir}/{asym}_neutron_all_PhT.png", dpi=300)
    plt.close()

print("Neutron asymmetry plots created!")