import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os


# ============================================================
# CONFIGURAZIONE
# ============================================================

# Scegli:
#   "zPt" -> binning z-P_hT, con 4 pannelli
#   "xQ2" -> binning xB-Q2, con un singolo pannello
binning = "xQ2"

# ============================================================
# CONFIGURAZIONE DIPENDENTE DAL BINNING
# ============================================================

if binning == "zPt":
    
    outdir = "ASYMMETRIES_zPt_plot_w_sys"
    # CSV asimmetrie
    asym_sum22_file = "RECO_CSV/output_RGC_NH3_asymmetries_sum22.csv"
    asym_fall22_file = "RECO_CSV/output_RGC_NH3_asymmetries_fall22.csv"

    # CSV sistematiche
    sys_sum22_file = "TOY_diagnostics/toy_summary_summer22_zPt.csv"
    sys_fall22_file = "TOY_diagnostics/toy_summary_fall22_zPt.csv"

    # colonne CSV asimmetrie
    asym_bin_column = "bin_zPt"
    x_column = "mean_z"

    # label asse x
    x_label = r"$z$"

    # larghezza box sistematica
    SYS_BOX_WIDTH = 0.035


elif binning == "xQ2":

    outdir = "ASYMMETRIES_xQ2_plot_w_sys"
    # CSV asimmetrie
    asym_sum22_file = "RECO_CSV/output_RGC_NH3_asymmetries_sum22_xQ2.csv"
    asym_fall22_file = "RECO_CSV/output_RGC_NH3_asymmetries_fall22_xQ2.csv"

    # CSV sistematiche
    sys_sum22_file = "TOY_diagnostics/toy_summary_summer22_xQ2.csv"
    sys_fall22_file = "TOY_diagnostics/toy_summary_fall22_xQ2.csv"

    # colonne CSV asimmetrie
    asym_bin_column = "bin_xQ2"
    x_column = "mean_xB"

    # label asse x
    x_label = r"$x_B$"

    # larghezza box sistematica
    SYS_BOX_WIDTH = 0.008


else:
    raise ValueError(
        f"Binning non riconosciuto: {binning}. "
        'Usare "zPt" oppure "xQ2".'
    )

os.makedirs(outdir, exist_ok=True)

# ============================================================
# CARICA CSV DELLE ASIMMETRIE
# ============================================================

sum22 = pd.read_csv(
    asym_sum22_file,
    skipinitialspace=True
)

fall22 = pd.read_csv(
    asym_fall22_file,
    skipinitialspace=True
)

sum22.columns = sum22.columns.str.strip()
fall22.columns = fall22.columns.str.strip()


# ============================================================
# CARICA CSV DELLE SISTEMATICHE
# ============================================================

sys_sum22 = pd.read_csv(
    sys_sum22_file,
    skipinitialspace=True
)

sys_fall22 = pd.read_csv(
    sys_fall22_file,
    skipinitialspace=True
)

sys_sum22.columns = sys_sum22.columns.str.strip()
sys_fall22.columns = sys_fall22.columns.str.strip()

sys_sum22["modulation"] = sys_sum22["modulation"].str.strip()
sys_fall22["modulation"] = sys_fall22["modulation"].str.strip()


# ============================================================
# MAPPA:
# nome asimmetria CSV RECO -> modulation CSV sistematiche
# ============================================================

asym_to_modulation = {
    "AUL_sinPhi":  "AUL_sinPhi",
    "AUL_sin2Phi": "AUL_sin2Phi",
    "ALL":          "ALL_const",
    "ALL_cosPhi":   "ALL_cosPhi",
    "ALU_sinPhi":   "ALU_sinPhi",
}


# ============================================================
# LISTA ASIMMETRIE
# ============================================================

asymmetries = [
    ("AUL_sinPhi", "AUL_sinPhi_err"),
    ("AUL_sin2Phi", "AUL_sin2Phi_err"),
    ("ALL", "ALL_err"),
    ("ALL_cosPhi", "ALL_cosPhi_err"),
    ("ALU_sinPhi", "ALU_sinPhi_err")
]


# ============================================================
# LABEL ASIMMETRIE
# ============================================================

asym_labels = {
    "AUL_sinPhi": r"$F_{UL}^{\sin\phi}/F_{UU}$",
    "AUL_sin2Phi": r"$F_{UL}^{\sin2\phi}/F_{UU}$",
    "ALL": r"$F_{LL}/F_{UU}$",
    "ALL_cosPhi": r"$F_{LL}^{\cos\phi}/F_{UU}$",
    "ALU_sinPhi": r"$F_{LU}^{\sin\phi}/F_{UU}$"
}


# ============================================================
# PhT GROUPS
# usati solamente per binning == "zPt"
# ============================================================

pt_groups = {
    "PhT_bin1": range(1, 8),
    "PhT_bin2": range(8, 15),
    "PhT_bin3": range(15, 22),
    "PhT_bin4": range(22, 26)
}


titles_asymmetries = [
    r"$0.0 < P_{hT} < 0.25 \ \mathrm{GeV}$",
    r"$0.25 < P_{hT} < 0.5 \ \mathrm{GeV}$",
    r"$0.5 < P_{hT} < 0.8 \ \mathrm{GeV}$",
    r"$0.8 < P_{hT} < 1.4 \ \mathrm{GeV}$"
]

# COLORI
q2_groups = [
    (range(1, 6),  "tab:blue",   r"$1 < Q^2 < 3\ \mathrm{GeV}^2$"),
    (range(6, 11), "tab:purple", r"$3 < Q^2 < 5\ \mathrm{GeV}^2$"),
    (range(11, 13), "deeppink",  r"$5 < Q^2 < 7\ \mathrm{GeV}^2$"),
    (range(13, 15), "tab:orange", r"$7 < Q^2 < 11\ \mathrm{GeV}^2$")
]

# ============================================================
# FUNZIONE PER MERGE ASIMMETRIE + SISTEMATICHE
# ============================================================

def merge_sys_for_asym(asym_df, sys_df, asym_name):

    modulation = asym_to_modulation[asym_name]

    # Nei nuovi toy_summary il bin si chiama semplicemente "bin"
    sel = sys_df[
        sys_df["modulation"] == modulation
    ][
        [
            "bin",
            "Total_sys_final",
            "Total_sys_final_err"
        ]
    ].copy()

    # Rinomina "bin" in modo che corrisponda al CSV delle asimmetrie
    #
    # zPt -> bin_zPt
    # xQ2 -> bin_xQ2
    sel = sel.rename(
        columns={
            "bin": asym_bin_column,
            "Total_sys_final": f"{asym_name}_totsys",
            "Total_sys_final_err": f"{asym_name}_totsys_err",
        }
    )

    merged = asym_df.merge(
        sel,
        on=asym_bin_column,
        how="left"
    )

    return merged


# ============================================================
# FUNZIONE PER I LIMITI Y
# ============================================================

def set_y_limits(ax, asym):

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


# ============================================================
# STILE COMUNE DEGLI ASSI
# ============================================================

def style_axis(ax):

    ax.axhline(
        0,
        color="black",
        linestyle="--",
        linewidth=1,
        alpha=0.6
    )

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


# ============================================================
# FUNZIONE PRINCIPALE DI PLOT
# ============================================================

def plot_campaign_with_sys(
    df,
    sys_df,
    campaign_label,
    color,
    marker,
    fname_suffix
):

    # loop sulle 5 asimmetrie
    for asym, err in asymmetries:

        # merge RECO + sistematiche
        sel_full = merge_sys_for_asym(
            df,
            sys_df,
            asym
        )

        totsys_col = f"{asym}_totsys"


        # ====================================================
        # zPt BINNING
        # 4 pannelli, uno per intervallo PhT
        # ====================================================

        if binning == "zPt":

            fig, axes = plt.subplots(
                1,
                4,
                figsize=(18, 3.5),
                sharey=True,
                sharex=True
            )

            for i, (label, bins) in enumerate(pt_groups.items()):

                sel = (
                    sel_full[
                        sel_full[asym_bin_column].isin(bins)
                    ]
                    .sort_values(x_column)
                )

                ax = axes[i]


                # --------------------------------------------
                # BOX SISTEMATICA
                # --------------------------------------------

                ax.bar(
                    sel[x_column],
                    height=2 * sel[totsys_col],
                    bottom=sel[asym] - sel[totsys_col],
                    width=2 * SYS_BOX_WIDTH,
                    color="r",
                    alpha=0.4,
                    edgecolor="none",
                    zorder=1
                )


                # --------------------------------------------
                # ASIMMETRIA + ERRORE STATISTICO
                # --------------------------------------------

                ax.errorbar(
                    sel[x_column],
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


                # --------------------------------------------
                # stile
                # --------------------------------------------

                style_axis(ax)
                set_y_limits(ax, asym)

                ax.set_xlim(0.2, 1)

                ax.set_title(
                    titles_asymmetries[i],
                    fontsize=12
                )

                ax.set_xlabel(x_label)

                if i == 0:
                    ax.set_ylabel(
                        asym_labels[asym],
                        fontsize=13
                    )


            # titolo generale
            fig.text(
                0.5,
                1.02,
                rf"$K^+$ NH$_3$ Run Group C — {campaign_label}",
                ha="center",
                fontsize=14
            )


            # --------------------------------------------
            # legenda
            # --------------------------------------------

            handles, labels = axes[0].get_legend_handles_labels()

            handles.append(
                Patch(
                    facecolor="r",
                    alpha=0.4,
                    edgecolor="none"
                )
            )

            labels.append("Systematic")

            legend_x = (
                0.155
                if campaign_label == "Summer22"
                else 0.135
            )

            fig.legend(
                handles,
                labels,
                bbox_to_anchor=(legend_x, 0.335)
            )


            plt.tight_layout()


            # --------------------------------------------
            # salva
            # --------------------------------------------

            output_file = os.path.join(
                outdir,
                f"{asym}_{fname_suffix}_{binning}.png"
            )

            plt.savefig(
                output_file,
                dpi=300,
                bbox_inches="tight"
            )

            plt.close()


        # ====================================================
        # xB-Q2 BINNING
        # un singolo pannello per asimmetria
        # ====================================================

        elif binning == "xQ2":

            sel = sel_full.sort_values(x_column)

            fig, ax = plt.subplots(figsize=(7, 5))


            # --------------------------------------------
            # BOX SISTEMATICA
            # tutti rossi
            # --------------------------------------------

            ax.bar(
                sel[x_column],
                height=2 * sel[totsys_col],
                bottom=sel[asym] - sel[totsys_col],
                width=2 * SYS_BOX_WIDTH,
                color="r",
                alpha=0.4,
                edgecolor="none",
                zorder=1
            )


            # --------------------------------------------
            # ASIMMETRIA + ERRORE STATISTICO
            # colore diverso per ogni intervallo Q2
            # --------------------------------------------

            for bins, q2_color, q2_label in q2_groups:

                sel_q2 = sel[
                    sel[asym_bin_column].isin(bins)
                ]

                if sel_q2.empty:
                    continue

                ax.errorbar(
                    sel_q2[x_column],
                    sel_q2[asym],
                    yerr=sel_q2[err],
                    fmt=marker,
                    markersize=6,
                    capsize=3,
                    linewidth=1.2,
                    color=q2_color,
                    markerfacecolor="white",
                    zorder=2,
                    label=q2_label
                )


            # --------------------------------------------
            # stile
            # --------------------------------------------

            style_axis(ax)
            set_y_limits(ax, asym)

            ax.set_xlabel(
                x_label,
                fontsize=13
            )

            ax.set_ylabel(
                asym_labels[asym],
                fontsize=13
            )


            # --------------------------------------------
            # titolo
            # --------------------------------------------

            ax.set_title(
                rf"$K^+$ NH$_3$ Run Group C — {campaign_label}",
                fontsize=14
            )


            # --------------------------------------------
            # legenda
            # --------------------------------------------

            handles, labels = ax.get_legend_handles_labels()

            handles.append(
                Patch(
                    facecolor="r",
                    alpha=0.4,
                    edgecolor="none"
                )
            )

            labels.append("Systematic")

            ax.legend(
                handles,
                labels,
                loc="best"
            )


            plt.tight_layout()


            # --------------------------------------------
            # salva
            # --------------------------------------------

            output_file = os.path.join(
                outdir,
                f"{asym}_{fname_suffix}_{binning}.png"
            )

            plt.savefig(
                output_file,
                dpi=300,
                bbox_inches="tight"
            )

            plt.close()


# ============================================================
# PRODUCI I PLOT
# ============================================================

print("\n============================================")
print(f" Binning selezionato: {binning}")
print("============================================\n")


# Summer22
plot_campaign_with_sys(
    sum22,
    sys_sum22,
    "Summer22",
    "tab:blue",
    "o",
    "sum22_with_sys"
)


# Fall22
plot_campaign_with_sys(
    fall22,
    sys_fall22,
    "Fall22",
    "tab:orange",
    "s",
    "fall22_with_sys"
)


# ============================================================
# FINE
# ============================================================

print(
    f"\nPlot Summer22 e Fall22 con sistematiche creati "
    f"per binning = {binning}"
)

print(
    f"Output directory: {os.path.abspath(outdir)}"
)