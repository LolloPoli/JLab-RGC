
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ============================================================
# CONFIGURAZIONE
# ============================================================

target = "NH3"
period = "fall22"   # "sum22", "fall22", "spring23"

csv_filename = (
    f"RECO_CSV/output_RGC_{target}_asymmetries_4D_{period}.csv"
)

outdir = f"asymmetries_4D_grid/{period}"
os.makedirs(outdir, exist_ok=True)

# Limiti degli assi
XLIM = (0.0, 0.7)

# Modificabili separatamente per ogni asimmetria
YLIMS = {
    "ALU_sinPhi":  (-0.3, 0.3),
    "ALL":         (-0.1, 0.8),
    "ALL_cosPhi":  (-0.5, 0.6),
    "AUL_sinPhi":  (-0.4, 0.4),
    "AUL_sin2Phi": (-0.5, 0.5),
}

# Riproduce il limite sugli errori usato nel codice ROOT
MAX_ERROR = 0.4

# ============================================================
# ASIMMETRIE
# ============================================================

asymmetries = {
    "ALU_sinPhi": (
        "ALU_sinPhi_err",
        r"$A_{LU}^{\sin\phi_h}$"
    ),
    "ALL": (
        "ALL_err",
        r"$A_{LL}$"
    ),
    "ALL_cosPhi": (
        "ALL_cosPhi_err",
        r"$A_{LL}^{\cos\phi_h}$"
    ),
    "AUL_sinPhi": (
        "AUL_sinPhi_err",
        r"$A_{UL}^{\sin\phi_h}$"
    ),
    "AUL_sin2Phi": (
        "AUL_sin2Phi_err",
        r"$A_{UL}^{\sin2\phi_h}$"
    ),
}

# ============================================================
# INTERVALLI CINEMATICI
# ============================================================

# Colonne: da sinistra a destra
z_ranges = [
    (0.20, 0.35),
    (0.35, 0.50),
    (0.50, 0.70),
    (0.70, 1.00),
]

# Righe: dall'alto verso il basso
pt_ranges = [
    (0.60, 1.40),
    (0.40, 0.60),
    (0.20, 0.40),
    (0.00, 0.20),
]

# Gruppi Q2 basati sul bin_xB (1-based)
q2_groups = [
    (1.0, 3.0, "tab:blue",   r"$1<Q^2<3$"),
    (3.0, 5.0, "tab:purple", r"$3<Q^2<5$"),
    (5.0, 7.0, "deeppink",   r"$5<Q^2<7$"),
    (7.0, 11.0, "tab:orange", r"$7<Q^2<11$"),
]

# ============================================================
# CARICAMENTO CSV
# ============================================================

df = pd.read_csv(csv_filename, skipinitialspace=True)
df.columns = df.columns.str.strip()

required = {"bin_xB", "bin_zPt", "mean_xB", "mean_Q2"}

for asym, (err, _) in asymmetries.items():
    required.update([asym, err])

missing = sorted(required - set(df.columns))

if missing:
    raise ValueError(
        f"Colonne mancanti nel CSV: {', '.join(missing)}"
    )

# Conversione a numerico; eventuali valori non validi -> NaN
for col in required:
    df[col] = pd.to_numeric(df[col], errors="coerce")

# ============================================================
# STILE GENERALE
# ============================================================

plt.rcParams.update({
    "font.size": 11,
    "axes.linewidth": 1.0,
    "mathtext.fontset": "stix",
    "font.family": "STIXGeneral",
    "savefig.dpi": 300,
})

# ============================================================
# LOOP SULLE CINQUE ASIMMETRIE
# ============================================================

for asym, (err, asym_label) in asymmetries.items():

    fig, axes = plt.subplots(
        4, 4,
        figsize=(18, 12),
        sharex=True,
        sharey=True,
        gridspec_kw={
            "wspace": 0.0,
            "hspace": 0.0,
        }
    )

    for row in range(4):
        for col in range(4):

            ax = axes[row, col]

            # Bin zPt:
            # riga superiore: 13,14,15,16
            # riga inferiore: 1,2,3,4
            bin_zPt = (3 - row) * 4 + col + 1

            selected = df[
                df["bin_zPt"] == bin_zPt
            ].copy()

            # Valori validi
            selected = selected.replace(
                [np.inf, -np.inf], np.nan
            ).dropna(
                subset=["bin_xB", "mean_xB", asym, err]
            )

            # Stesso taglio sulle ampiezze del codice ROOT
            selected = selected[
                selected[asym].between(-1.0, 1.0)
                & (selected[err] >= 0)
            ]

            # Limite sugli errori
            selected[err] = selected[err].clip(
                upper=MAX_ERROR
            )

            # --------------------------------------------
            # Punti colorati per intervallo Q2
            # --------------------------------------------

            for q2min, q2max, color, q2_label in q2_groups:

                subset = selected[
                    (selected["mean_Q2"] >= q2min) &
                    (selected["mean_Q2"] < q2max)
                ].sort_values("mean_xB")

                if subset.empty:
                    continue

                ax.errorbar(
                    subset["mean_xB"],
                    subset[asym],
                    yerr=subset[err],
                    fmt="o",
                    markersize=3.5,
                    capsize=2,
                    elinewidth=0.9,
                    color=color,
                    markerfacecolor="white",
                    zorder=3,
                )

            # --------------------------------------------
            # Assi e linea a zero
            # --------------------------------------------

            ax.set_xlim(*XLIM)
            ax.set_ylim(*YLIMS[asym])

            ax.axhline(
                0,
                color="gray",
                linestyle="--",
                linewidth=0.8,
                alpha=0.7,
                zorder=1,
            )

            ax.tick_params(
                direction="in",
                top=True,
                right=True,
                length=4,
                width=0.9,
                labelsize=10,
            )

            # Numeri dell'asse X solo in basso
            ax.tick_params(
                axis="x",
                labelbottom=(row == 3)
            )

            # Numeri dell'asse Y solo a sinistra
            ax.tick_params(
                axis="y",
                labelleft=(col == 0)
            )

            # Titoli X solo nella riga inferiore
            if row == 3:
                ax.set_xlabel(
                    r"$x_B$",
                    fontsize=14,
                    labelpad=7,
                )

            # Titoli Y solo nella prima colonna
            if col == 0:
                ax.set_ylabel(
                    asym_label,
                    fontsize=14,
                    labelpad=9,
                )

            # Nessun titolo dentro i pannelli

            # Griglia leggera
            ax.grid(
                alpha=0.15,
                linestyle=":",
                linewidth=0.7,
            )

            # Bordi dei pannelli
            for spine in ax.spines.values():
                spine.set_linewidth(1.0)

    # ========================================================
    # LAYOUT: MATRICE COMPATTA
    # ========================================================

    fig.subplots_adjust(
        left=0.095,
        right=0.91,
        bottom=0.085,
        top=0.88, 
        wspace=0.0,
        hspace=0.0,
    )

    # ========================================================
    # ETICHETTE z SOPRA LE COLONNE
    # ========================================================

    for col, (zmin, zmax) in enumerate(z_ranges):

        pos = axes[0, col].get_position()
        xcenter = (pos.x0 + pos.x1) / 2

        # Etichette z sopra le colonne
        fig.text(
            xcenter,
            0.885,
            rf"${zmin:.2f}<z<{zmax:.2f}$",
            ha="center",
            va="bottom",
            fontsize=13,
        )

    # ========================================================
    # ETICHETTE Pt A DESTRA, RUOTATE
    # ========================================================

    for row, (ptmin, ptmax) in enumerate(pt_ranges):

        pos = axes[row, 3].get_position()
        ycenter = (pos.y0 + pos.y1) / 2

        fig.text(
            0.925,
            ycenter,
            rf"${ptmin:.2f}<P_{{hT}}<{ptmax:.2f}"
            r"\ \mathrm{GeV}$",
            rotation=-90,
            ha="center",
            va="center",
            fontsize=13,
        )

    # ========================================================
    # TITOLO GENERALE
    # ========================================================

    # Titolo generale
    fig.suptitle(
        rf"$K^+$ {target} Run Group C"
        + "  -  "
        + asym_label + r" vs $x_B$ in 4D bins of $(x_B, Q^2, z, P_{hT})$",
        fontsize=20,
        y=0.98,
    )

    # ========================================================
    # LEGENDA UNICA
    # ========================================================

    legend_handles = [
        axes[0, 0].errorbar(
            [np.nan], [np.nan],
            yerr=[0.1],
            fmt="o",
            color=color,
            markerfacecolor="white",
            markersize=6,
            capsize=3,
            elinewidth=1.2,
            linestyle="none",
            label=label + r"$\ \mathrm{GeV}^2$",
        )
        for _, _, color, label in q2_groups
    ]

    fig.legend(
        handles=legend_handles,
        loc="center",
        bbox_to_anchor=(0.50, 0.925),
        ncol=4,
        fontsize=14,
        frameon=False,
        columnspacing=2.0,
        handletextpad=0.4,
    )

    # ========================================================
    # SALVATAGGIO
    # ========================================================

    output_base = os.path.join(
        outdir,
        f"{asym}_4D_grid_{period}"
    )

    fig.savefig(output_base + ".png",dpi=300,bbox_inches="tight",)

    #fig.savefig(output_base + ".pdf",bbox_inches="tight",)

    plt.close(fig)

    print(f"Created: {output_base}.png / .pdf")

print("All 4D grids created in:", outdir)
