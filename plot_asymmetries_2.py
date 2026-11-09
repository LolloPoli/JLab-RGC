import pandas as pd
import matplotlib.pyplot as plt
import os

# ============================================================
# FLAG
# True  -> plot Summer22 + Fall22 + Spring23
# False -> plot Summer22 + Fall22
# ============================================================

include_spring23 = True

# Mostra anche le modulazioni non polarizzate (richiede le colonne AUU nei CSV).
include_AUU = True


# cartella output
outdir = "overlayed_asymmetries"
os.makedirs(outdir, exist_ok=True)


# ============================================================
# Carica CSV
# ============================================================

sum22 = pd.read_csv("RECO_CSV/output_RGC_NH3_asymmetries_sum22.csv",skipinitialspace=True)
fall22 = pd.read_csv("RECO_CSV/output_RGC_NH3_asymmetries_fall22.csv",skipinitialspace=True)

# Spring23 viene caricato solo se richiesto
if include_spring23:
    spring23 = pd.read_csv("RECO_CSV/output_RGC_NH3_asymmetries_spring23.csv",skipinitialspace=True)


sum22_xB = pd.read_csv("RECO_CSV/output_RGC_NH3_asymmetries_sum22_xQ2.csv",skipinitialspace=True)
fall22_xB = pd.read_csv("RECO_CSV/output_RGC_NH3_asymmetries_fall22_xQ2.csv",skipinitialspace=True)

if include_spring23:
    spring23_xB = pd.read_csv("RECO_CSV/output_RGC_NH3_asymmetries_spring23_xQ2.csv",skipinitialspace=True)

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


if include_AUU:
    asymmetries += [
        ("AUU_cosPhi", "AUU_cosPhi_err"),
        ("AUU_cos2Phi", "AUU_cos2Phi_err"),
    ]

# Controlla i CSV prima di produrre i grafici.
input_tables = [
    ("Summer22 zPt", sum22),
    ("Fall22 zPt", fall22),
    ("Summer22 xQ2", sum22_xB),
    ("Fall22 xQ2", fall22_xB),
]
if include_spring23:
    input_tables += [
        ("Spring23 zPt", spring23),
        ("Spring23 xQ2", spring23_xB),
    ]

for table_name, table in input_tables:
    table.columns = table.columns.str.strip()
    required = {col for pair in asymmetries for col in pair}
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(
            f"{table_name}: colonne mancanti: {', '.join(missing)}. "
            "Esporta valori ed errori AUU dal fit aggiornato oppure imposta "
            "include_AUU = False per mostrare solo le ampiezze polarizzate."
        )


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
    "ALU_sinPhi": r"$F_{LU}^{\sin\phi}/F_{UU}$",
    "AUU_cosPhi": r"$F_{UU}^{\cos\phi}/F_{UU}$",
    "AUU_cos2Phi": r"$F_{UU}^{\cos 2\phi}/F_{UU}$"
}


# ============================================================
# Loop sulle asimmetrie
# ============================================================

for asym, err in asymmetries:
    # Piccoli offset orizzontali per separare le campagne
    Z_OFFSET = 0.01

    offset_summer22 = -Z_OFFSET
    offset_fall22 = 0.0
    offset_spring23 = +Z_OFFSET

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
                spring23_sel["mean_z"] + offset_spring23,
                spring23_sel[asym],
                yerr=spring23_sel[err],
                fmt="^",
                markersize=5,
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
            sum22_sel["mean_z"] + offset_summer22,
            sum22_sel[asym],
            yerr=sum22_sel[err],
            fmt="o",
            markersize=5,
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
            fall22_sel["mean_z"] + offset_fall22,
            fall22_sel[asym],
            yerr=fall22_sel[err],
            fmt="s",
            markersize=5,
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

        if asym.startswith("AUU_"):
            ax.margins(y=0.15)

        elif asym == "ALU_sinPhi":
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

    if asym.startswith("AUU_"):
        fig.legend(handles,labels,bbox_to_anchor=(0.135, 0.865))
    else:
        fig.legend(handles,labels,bbox_to_anchor=(0.135, 0.385))


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
    


# ============================================================
# Plot xB-Q2
# ============================================================

for asym, err in asymmetries:

    fig, ax = plt.subplots(
        figsize=(7, 5)
    )

    # --------------------------------------------------------
    # Spring23 - prima, così rimane sotto
    # --------------------------------------------------------

    if include_spring23:

        spring23_sel = spring23_xB.sort_values("mean_xB")

        ax.errorbar(
            spring23_sel["mean_xB"],
            spring23_sel[asym],
            yerr=spring23_sel[err],
            fmt="^",
            markersize=5,
            capsize=3,
            linewidth=1.2,
            label="Spring23",
            color="tab:green",
            markerfacecolor="white",
            zorder=1
        )


    # --------------------------------------------------------
    # Summer22
    # --------------------------------------------------------

    sum22_sel = sum22_xB.sort_values("mean_xB")

    ax.errorbar(
        sum22_sel["mean_xB"],
        sum22_sel[asym],
        yerr=sum22_sel[err],
        fmt="o",
        markersize=5,
        capsize=3,
        linewidth=1.2,
        label="Summer22",
        color="tab:blue",
        markerfacecolor="white",
        zorder=2
    )


    # --------------------------------------------------------
    # Fall22
    # --------------------------------------------------------

    fall22_sel = fall22_xB.sort_values("mean_xB")

    ax.errorbar(
        fall22_sel["mean_xB"],
        fall22_sel[asym],
        yerr=fall22_sel[err],
        fmt="s",
        markersize=5,
        capsize=3,
        linewidth=1.2,
        label="Fall22",
        color="tab:orange",
        markerfacecolor="white",
        zorder=3
    )


    # zero line
    ax.axhline(
        0,
        color="black",
        linestyle="--",
        linewidth=1,
        alpha=0.6
    )


    # stessi limiti dei plot precedenti
    if asym.startswith("AUU_"):
        ax.margins(y=0.15)

    elif asym == "ALU_sinPhi":
        ax.set_ylim(-0.2, 0.2)

    elif asym == "ALL":
        ax.set_ylim(-0.1, 0.8)

    elif asym == "ALL_cosPhi":
        ax.set_ylim(-0.6, 0.6)

    elif asym == "AUL_sin2Phi":
        ax.set_ylim(-0.4, 0.5)

    else:
        ax.set_ylim(-0.3, 0.3)


    # assi
    ax.set_xlabel(r"$x_B$", fontsize=13)
    ax.set_ylabel(asym_labels[asym], fontsize=13)

    ax.set_title(
        r"$K^+$ NH$_3$ Run Group C",
        fontsize=14
    )


    # stile
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

    ax.legend()


    plt.tight_layout()


    # output
    if include_spring23:
        output_name = f"{asym}_xB_Q2_campaign_comparison_all.png"
    else:
        output_name = f"{asym}_xB_Q2_campaign_comparison_sum_fall.png"


    plt.savefig(
        f"{outdir}/{output_name}",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()


print("All asymmetry plots created in:", outdir)

# ============================================================
# Plot xB-Q2: tre campagne affiancate per ogni modulazione
# Solo errori statistici, colori per intervallo Q2
# ============================================================

comparison_outdir = os.path.join(outdir, "xB_Q2_side_by_side")
os.makedirs(comparison_outdir, exist_ok=True)

asym_bin_column = "bin_xQ2"
q2_groups = [
    (range(1, 6), "tab:blue", r"$1 < Q^2 < 3\ \mathrm{GeV}^2$"),
    (range(6, 11), "tab:purple", r"$3 < Q^2 < 5\ \mathrm{GeV}^2$"),
    (range(11, 13), "deeppink", r"$5 < Q^2 < 7\ \mathrm{GeV}^2$"),
    (range(13, 15), "tab:orange", r"$7 < Q^2 < 11\ \mathrm{GeV}^2$"),
]

campaigns_xq2 = [
    ("Summer22", sum22_xB, "o"),
    ("Fall22", fall22_xB, "s"),
]
if include_spring23:
    campaigns_xq2.append(("Spring23", spring23_xB, "^"))

for campaign_label, campaign_df, _ in campaigns_xq2:
    for col in ("mean_xB", asym_bin_column):
        if col not in campaign_df.columns:
            raise KeyError(f"Colonna '{col}' mancante nei dati {campaign_label}")

for asym, err in asymmetries:
    fig, axes = plt.subplots(
        1, len(campaigns_xq2),
        figsize=(6 * len(campaigns_xq2), 5),
        sharex=True, sharey=True,
        squeeze=False,
    )
    axes = axes[0]

    for ax, (campaign_label, campaign_df, marker) in zip(axes, campaigns_xq2):
        sel = campaign_df.sort_values("mean_xB")

        for bins, q2_color, q2_label in q2_groups:
            sel_q2 = sel[sel[asym_bin_column].isin(bins)]
            if sel_q2.empty:
                continue
            ax.errorbar(
                sel_q2["mean_xB"],
                sel_q2[asym],
                yerr=sel_q2[err],
                fmt=marker,
                markersize=5,
                capsize=3,
                linewidth=1,
                color=q2_color,
                markerfacecolor="white",
                zorder=2,
            )

        ax.axhline(0, color="black", linestyle="--", linewidth=1, alpha=0.6)
        ax.set_xlabel(r"$x_B$", fontsize=13)
        ax.set_title(campaign_label, fontsize=14)
        ax.tick_params(direction="in", top=True, right=True, length=5)
        for spine in ax.spines.values():
            spine.set_linewidth(1.2)
        ax.grid(alpha=0.25, linestyle=":")

    # Stessi limiti y su tutti e tre i pannelli
    if asym.startswith("AUU_"):
        axes[0].margins(y=0.15)
    elif asym == "ALU_sinPhi":
        axes[0].set_ylim(-0.2, 0.2)
    elif asym == "ALL":
        axes[0].set_ylim(-0.1, 0.8)
    elif asym == "ALL_cosPhi":
        axes[0].set_ylim(-0.6, 0.6)
    elif asym == "AUL_sin2Phi":
        axes[0].set_ylim(-0.4, 0.5)
    else:
        axes[0].set_ylim(-0.3, 0.3)

    axes[0].set_ylabel(asym_labels[asym], fontsize=13)
    fig.suptitle(r"$K^+$ NH$_3$ Run Group C", fontsize=15, y=0.92)

    # Legenda unica, anche se un intervallo Q2 non ha punti
    from matplotlib.lines import Line2D
    legend_handles = [
        Line2D([0], [0], color=color, marker="o", linestyle="none",
               markerfacecolor="white", markersize=7, label=label)
        for _, color, label in q2_groups
    ]
    # ========================================================
    # Legenda dentro il pannello centrale
    # ========================================================

    legend_ax = axes[1] if len(axes) > 1 else axes[0]
    legend_ax.legend(
        handles=legend_handles,
        loc="upper left",
        fontsize=11,
        ncol=1,
        frameon=True,
        framealpha=0.85,
        edgecolor="lightgray",
        handletextpad=0.4,
        columnspacing=0.8
    )

    # Più spazio dedicato ai grafici
    fig.tight_layout(rect=(0, 0, 1, 0.93))

    output_file = os.path.join(comparison_outdir, f"{asym}_xB_Q2_side_by_side.png")
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)

print("Side-by-side xB-Q2 plots created in:", comparison_outdir)
