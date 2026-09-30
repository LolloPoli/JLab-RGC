import pandas as pd
import matplotlib.pyplot as plt
import os

# ============================================================
# FLAG
# True  -> plot Summer22 + Fall22 + Spring23
# False -> plot Summer22 + Fall22
# ============================================================

include_spring23 = True


# cartella output
outdir = "overlayed_asymmetries"
os.makedirs(outdir, exist_ok=True)


# ============================================================
# Carica CSV
# ============================================================

sum22 = pd.read_csv(
    "output_RGC_NH3_asymmetries_sum22.csv",
    skipinitialspace=True
)

fall22 = pd.read_csv(
    "output_RGC_NH3_asymmetries_fall22.csv",
    skipinitialspace=True
)

# Spring23 viene caricato solo se richiesto
if include_spring23:
    spring23 = pd.read_csv(
        "output_RGC_NH3_asymmetries_spring23.csv",
        skipinitialspace=True
    )


# ============================================================
# Gruppi di PhT
# ============================================================

pt_groups = {
    "PhT_bin1": range(1, 8),
    "PhT_bin2": range(8, 15),
    "PhT_bin3": range(15, 22),
    "PhT_bin4": range(22, 26)
}


# ============================================================
# Lista asimmetrie
# ============================================================

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


# ============================================================
# Loop sulle asimmetrie
# ============================================================

for asym, err in asymmetries:

    fig, axes = plt.subplots(
        1, 4,
        figsize=(18, 3.5),
        sharey=True,
        sharex=True
    )

    for i, (label, bins) in enumerate(pt_groups.items()):

        # -------------------------
        # Selezione bin
        # -------------------------

        sum22_sel = (
            sum22[sum22["bin_zPt"].isin(bins)]
            .sort_values("mean_z")
        )

        fall22_sel = (
            fall22[fall22["bin_zPt"].isin(bins)]
            .sort_values("mean_z")
        )

        if include_spring23:
            spring23_sel = (
                spring23[spring23["bin_zPt"].isin(bins)]
                .sort_values("mean_z")
            )

        ax = axes[i]


    
        # ====================================================
        # Spring23
        # Solo se include_spring23 == True
        # ====================================================

        if include_spring23:

            ax.errorbar(
                spring23_sel["mean_z"],
                spring23_sel[asym],
                yerr=spring23_sel[err],
                fmt="^",
                markersize=6,
                capsize=3,
                linewidth=1.2,
                label="Spring23",
                color="tab:green",
                markerfacecolor="white"
            )
        
        # ====================================================
        # Summer22
        # ====================================================

        ax.errorbar(
            sum22_sel["mean_z"],
            sum22_sel[asym],
            yerr=sum22_sel[err],
            fmt="o",
            markersize=6,
            capsize=3,
            linewidth=1.2,
            label="Summer22",
            color="tab:blue",
            markerfacecolor="white"
        )


        # ====================================================
        # Fall22
        # ====================================================

        ax.errorbar(
            fall22_sel["mean_z"],
            fall22_sel[asym],
            yerr=fall22_sel[err],
            fmt="s",
            markersize=6,
            capsize=3,
            linewidth=1.2,
            label="Fall22",
            color="tab:orange",
            markerfacecolor="white"
        )

        # ====================================================
        # Zero line
        # ====================================================

        ax.axhline(
            0,
            color="black",
            linestyle="--",
            linewidth=1,
            alpha=0.6
        )


        # ====================================================
        # Y limits
        # ====================================================

        if asym == "ALU_sinPhi":
            ax.set_ylim(-0.2, 0.2)

        elif asym == "ALL":
            ax.set_ylim(-0.1, 0.8)

        elif asym == "ALL_cosPhi":
            ax.set_ylim(-0.6, 0.6)

        elif asym == "AUL_sin2Phi":
            ax.set_ylim(-0.4, 0.5)

        else:
            ax.set_ylim(-0.3, 0.3)


        ax.set_xlim(0.2, 1)


        # ====================================================
        # Titoli / assi
        # ====================================================

        ax.set_title(
            titles_asymmetries[i],
            fontsize=12
        )

        ax.set_xlabel(r"$z$")

        if i == 0:
            ax.set_ylabel(
                asym_labels[asym],
                fontsize=13
            )


        # ====================================================
        # Stile
        # ====================================================

        ax.tick_params(
            direction="in",
            top=True,
            right=True,
            length=5
        )

        for spine in ax.spines.values():
            spine.set_linewidth(1.2)

        ax.grid(
            alpha=0.25,
            linestyle=":"
        )


    # ========================================================
    # Titolo generale
    # ========================================================

    fig.text(
        0.5,
        1.02,
        r"$K^+$ NH$_3$ Run Group C",
        ha="center",
        fontsize=14
    )


    # ========================================================
    # Legenda unica
    # ========================================================

    handles, labels = axes[0].get_legend_handles_labels()

    fig.legend(
        handles,
        labels,
        bbox_to_anchor=(0.135, 0.385)
    )


    plt.tight_layout()


    # ========================================================
    # Nome output diverso a seconda della flag
    # ========================================================

    if include_spring23:
        output_name = f"{asym}_campaign_comparison_all.png"
    else:
        output_name = f"{asym}_campaign_comparison_sum_fall.png"


    plt.savefig(
        f"{outdir}/{output_name}",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()


print("All asymmetry plots created!")