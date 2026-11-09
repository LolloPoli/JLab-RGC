#!/usr/bin/env python3
"""Fall22: combinazioni di elicita' dalle ampiezze di spin e, se presenti, UU.

Richiede numpy, pandas, matplotlib. Le A in ingresso devono essere gia'
corrette per polarizzazioni, diluizione e depolarizzazione, con denominatore
comune U0 = F_UU,T + epsilon F_UU,L e convenzioni di segno coerenti.
Non si applicano ulteriori fattori di depolarizzazione alle A dei CSV.

Esempio, eseguito dalla cartella rgc:
    python3 plot_helicity_sf.py --input-dir RECO_CSV

Le colonne canoniche di covarianza vengono riconosciute automaticamente:
    cov_ALU_AUL_sinPhi           = Cov(ALU_sinPhi, AUL_sinPhi)
    cov_ALL_cosPhi_AUL_sinPhi    = Cov(ALL_cosPhi, AUL_sinPhi)
Se una colonna manca, per quella coppia si assume rho=0 con avviso.
Per specificare nomi diversi:
    python3 plot_helicity_sf.py --input-dir RECO_CSV \
        --cov-column mia_cov_LU_UL --complex-cov-column mia_cov_LLcos_UL
In alternativa --rho e --complex-rho impostano correlazioni costanti.
Gli errori propagati sono quelli delle colonne *_err, senza aggiungere
sistematiche non presenti nei file. Modulo e fase: propagazione lineare;
la fase non viene mostrata quando l'ellisse gaussiana congiunta al 68.3%
contiene l'origine. Le combinazioni non applicano la correzione fascio-fotone.
Per mostrare tutte le fasi definite, senza il criterio sull'ellisse:
    python3 plot_helicity_sf.py --input-dir RECO_CSV --phase-mode all \
        --outdir helicity_sf_fall22_all_phases
"""

import argparse
from pathlib import Path
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd


PT_GROUPS = [
    (range(1, 8), r"$0 < P_{hT} < 0.25\ \mathrm{GeV}$"),
    (range(8, 15), r"$0.25 < P_{hT} < 0.5\ \mathrm{GeV}$"),
    (range(15, 22), r"$0.5 < P_{hT} < 0.8\ \mathrm{GeV}$"),
    (range(22, 26), r"$0.8 < P_{hT} < 1.4\ \mathrm{GeV}$"),
]
Q2_GROUPS = [
    (range(1, 6), "tab:blue", r"$1 < Q^2 < 3\ \mathrm{GeV}^2$"),
    (range(6, 11), "tab:purple", r"$3 < Q^2 < 5\ \mathrm{GeV}^2$"),
    (range(11, 13), "deeppink", r"$5 < Q^2 < 7\ \mathrm{GeV}^2$"),
    (range(13, 15), "tab:orange", r"$7 < Q^2 < 11\ \mathrm{GeV}^2$"),
]
OBSERVABLES = [
    ("ImGpp_over_U0", r"$\mathrm{Im}\,\mathcal{G}_{+1,0}^{++}/U_0$"),
    ("ImGmm_over_U0", r"$\mathrm{Im}\,\mathcal{G}_{+1,0}^{--}/U_0$"),
    ("DeltaG_T_over_U0", r"$(\mathcal{G}_{+1,+1}^{++}-\mathcal{G}_{+1,+1}^{--})/U_0$"),
    ("ImG_TT_over_U0", r"$\mathrm{Im}\,\mathcal{G}_{+1,-1}^{++}/U_0$"),
    ("ReDeltaG_LT_over_U0", r"$\mathrm{Re}\,\Delta\mathcal{G}_{LT}/U_0$"),
    ("ImDeltaG_LT_over_U0", r"$\mathrm{Im}\,\Delta\mathcal{G}_{LT}/U_0$"),
    ("AbsDeltaG_LT_over_U0", r"$|\Delta\mathcal{G}_{LT}|/U_0$"),
    ("PhaseDeltaG_LT_deg", r"$\arg\Delta\mathcal{G}_{LT}\ [\mathrm{deg}]$"),
    ("ReGpp_over_U0", r"$\mathrm{Re}\,\mathcal{G}_{+1,0}^{++}/U_0$"),
    ("ReGmm_over_U0", r"$\mathrm{Re}\,\mathcal{G}_{+1,0}^{--}/U_0$"),
    ("ReG_TT_over_U0", r"$\mathrm{Re}\,\mathcal{G}_{+1,-1}^{++}/U_0$"),
    ("AbsGpp_over_U0", r"$|\mathcal{G}_{+1,0}^{++}|/U_0$"),
    ("PhaseGpp_deg", r"$\arg\mathcal{G}_{+1,0}^{++}\ [\mathrm{deg}]$"),
    ("AbsGmm_over_U0", r"$|\mathcal{G}_{+1,0}^{--}|/U_0$"),
    ("PhaseGmm_deg", r"$\arg\mathcal{G}_{+1,0}^{--}\ [\mathrm{deg}]$"),
    ("AbsG_TT_over_U0", r"$|\mathcal{G}_{+1,-1}^{++}|/U_0$"),
    ("PhaseG_TT_deg", r"$\arg\mathcal{G}_{+1,-1}^{++}\ [\mathrm{deg}]$"),
]
TITLE = r"$K^+$ NH$_3$ Run Group C - Fall22"
FORMULA_TITLES = {
    "ImGpp_over_U0": (
        r"$\mathrm{Im}\,\mathcal{G}_{+1,0}^{++}/U_0"
        r" = -(A_{UL}^{\sin\phi_h}+A_{LU}^{\sin\phi_h})/\sqrt{2}$"
    ),
    "ImGmm_over_U0": (
        r"$\mathrm{Im}\,\mathcal{G}_{+1,0}^{--}/U_0"
        r" = (A_{UL}^{\sin\phi_h}-A_{LU}^{\sin\phi_h})/\sqrt{2}$"
    ),
    "DeltaG_T_over_U0": r"$(\mathcal{G}_{+1,+1}^{++}-\mathcal{G}_{+1,+1}^{--})/U_0=2A_{LL}$",
    "ImG_TT_over_U0": r"$\mathrm{Im}\,\mathcal{G}_{+1,-1}^{++}/U_0=-A_{UL}^{\sin2\phi_h}$",
    "ReDeltaG_LT_over_U0": r"$\mathrm{Re}\,\Delta\mathcal{G}_{LT}/U_0=-\sqrt{2}\,A_{LL}^{\cos\phi_h}$",
    "ImDeltaG_LT_over_U0": r"$\mathrm{Im}\,\Delta\mathcal{G}_{LT}/U_0=-\sqrt{2}\,A_{UL}^{\sin\phi_h}$",
    "AbsDeltaG_LT_over_U0": r"$|\Delta\mathcal{G}_{LT}|/U_0=\sqrt{2}\sqrt{(A_{LL}^{\cos\phi_h})^2+(A_{UL}^{\sin\phi_h})^2}$",
    "PhaseDeltaG_LT_deg": r"$\arg\Delta\mathcal{G}_{LT}=\mathrm{atan2}(-A_{UL}^{\sin\phi_h},-A_{LL}^{\cos\phi_h})$",
    "ReGpp_over_U0": r"$\mathrm{Re}\,\mathcal{G}_{+1,0}^{++}/U_0=-(A_{UU}^{\cos\phi_h}+A_{LL}^{\cos\phi_h})/\sqrt{2}$",
    "ReGmm_over_U0": r"$\mathrm{Re}\,\mathcal{G}_{+1,0}^{--}/U_0=-(A_{UU}^{\cos\phi_h}-A_{LL}^{\cos\phi_h})/\sqrt{2}$",
    "ReG_TT_over_U0": r"$\mathrm{Re}\,\mathcal{G}_{+1,-1}^{++}/U_0=-A_{UU}^{\cos2\phi_h}$",
}
COMPLEX_TITLE = (
    r"$\Delta\mathcal{G}_{LT}/U_0"
    r"=(\mathcal{G}_{+1,0}^{++}-\mathcal{G}_{+1,0}^{--})/U_0"
    r"=-\sqrt{2}(A_{LL}^{\cos\phi_h}+iA_{UL}^{\sin\phi_h})$"
)
COMPLEX_FAMILIES = {
    "DeltaG_LT": COMPLEX_TITLE,
    "Gpp": (r"$\mathcal{G}_{+1,0}^{++}/U_0=-[(A_{UU}^{\cos\phi_h}+A_{LL}^{\cos\phi_h})"
            r"+i(A_{LU}^{\sin\phi_h}+A_{UL}^{\sin\phi_h})]/\sqrt{2}$"),
    "Gmm": (r"$\mathcal{G}_{+1,0}^{--}/U_0=-[(A_{UU}^{\cos\phi_h}-A_{LL}^{\cos\phi_h})"
            r"+i(A_{LU}^{\sin\phi_h}-A_{UL}^{\sin\phi_h})]/\sqrt{2}$"),
    "G_TT": r"$\mathcal{G}_{+1,-1}^{++}/U_0=-A_{UU}^{\cos2\phi_h}-iA_{UL}^{\sin2\phi_h}$",
}
for _stem, _symbol in (("Gpp", r"\mathcal{G}_{+1,0}^{++}"),
                        ("Gmm", r"\mathcal{G}_{+1,0}^{--}"),
                        ("G_TT", r"\mathcal{G}_{+1,-1}^{++}")):
    FORMULA_TITLES[f"Abs{_stem}_over_U0"] = (
        r"$|" + _symbol + r"|/U_0=\sqrt{(\mathrm{Re}\," + _symbol
        + r"/U_0)^2+(\mathrm{Im}\," + _symbol + r"/U_0)^2}$")
    FORMULA_TITLES[f"Phase{_stem}_deg"] = (
        r"$\arg " + _symbol + r"=\mathrm{atan2}(\mathrm{Im}\,"
        + _symbol + r",\mathrm{Re}\," + _symbol + r")$")
# Quantile 68.27% del chi-quadro a 2 gradi di liberta'.
ELLIPSE_LEVEL = -2 * np.log(1 - 0.682689492)


def covariance_values(df, column, rho, canonical, a, b, path):
    if column is None and rho is None and canonical in df:
        column = canonical
    if column:
        if column not in df:
            raise ValueError(f"{path}: colonna di covarianza assente: {column}")
        cov = pd.to_numeric(df[column], errors="raise").to_numpy(dtype=float)
        source = column
    else:
        rho = 0.0 if rho is None else rho
        cov = rho * df[a + "_err"].to_numpy() * df[b + "_err"].to_numpy()
        source = f"assumed_rho={rho:g}"
        warnings.warn(f"{path.name}: {canonical} non usata; rho={rho:g} per {a}, {b}.")
    bound = df[a + "_err"].to_numpy() * df[b + "_err"].to_numpy()
    invalid = ~np.isfinite(cov) | (np.abs(cov) > bound * (1 + 1e-10))
    if invalid.any():
        details = []
        for index in np.flatnonzero(invalid)[:5]:
            row = df.iloc[index]
            bin_id = ", ".join(f"{col}={row[col]:g}" for col in ("bin_zPt", "bin_xQ2") if col in df)
            correlation = cov[index] / bound[index] if bound[index] > 0 else np.nan
            details.append(
                f"{bin_id or f'riga {index + 1}'}: cov={cov[index]:.9g}, "
                f"sigma({a})={row[a + '_err']:.9g}, sigma({b})={row[b + '_err']:.9g}, "
                f"sigma_a*sigma_b={bound[index]:.9g}, rho={correlation:.9g}"
            )
        raise ValueError(f"{path}: covarianza {source} non valida (|rho| > 1 o non finita).\n"
                         + "\n".join(details))
    return np.clip(cov, -bound, bound), source


def optional_covariance(df, a, b, path, missing):
    """Accetta cov_<nome parametro>_<nome parametro>, in entrambi gli ordini."""
    canonical = f"cov_{a}_{b}"
    column = next((col for col in (canonical, f"cov_{b}_{a}") if col in df), None)
    if column:
        cov, source = covariance_values(df, column, None, canonical, a, b, path)
    else:
        cov, source = np.zeros(len(df)), "assumed_rho=0"
        missing.append(canonical)
    df[canonical + "_used"] = cov
    df[canonical + "_source"] = source
    return cov


def add_polar_quantities(df, stem):
    """Propagazione locale della covarianza cartesiana in modulo e fase."""
    re = df[f"Re{stem}_over_U0"].to_numpy()
    im = df[f"Im{stem}_over_U0"].to_numpy()
    vr = df[f"Re{stem}_over_U0_err"].to_numpy()**2
    vi = df[f"Im{stem}_over_U0_err"].to_numpy()**2
    cri = df[f"cov_Re_Im_{stem}_over_U0"].to_numpy()
    if (np.abs(cri) > np.sqrt(vr * vi) * (1 + 1e-10)).any():
        raise ValueError(f"{stem}: covarianza Re/Im non positiva semidefinita; "
                         "controllare coerenza delle covarianze fornite e di quelle mancanti.")
    radius = np.hypot(re, im)
    var_radius = np.full(len(df), np.nan)
    var_phase = np.full(len(df), np.nan)
    nonzero = radius > 0
    var_radius[nonzero] = (
        re[nonzero]**2 * vr[nonzero] + im[nonzero]**2 * vi[nonzero]
        + 2 * re[nonzero] * im[nonzero] * cri[nonzero]
    ) / radius[nonzero]**2
    var_phase[nonzero] = (
        im[nonzero]**2 * vr[nonzero] + re[nonzero]**2 * vi[nonzero]
        - 2 * re[nonzero] * im[nonzero] * cri[nonzero]
    ) / radius[nonzero]**4
    df[f"Abs{stem}_over_U0"] = radius
    df[f"Abs{stem}_over_U0_err"] = np.sqrt(np.maximum(0, var_radius))
    df[f"Phase{stem}_deg"] = np.where(nonzero, np.degrees(np.arctan2(im, re)), np.nan)
    df[f"Phase{stem}_deg_err"] = np.degrees(np.sqrt(np.maximum(0, var_phase)))
    determinant = vr * vi - cri**2
    distance2 = np.full(len(df), np.nan)
    invertible = determinant > 1e-12 * vr * vi
    distance2[invertible] = (
        re[invertible]**2 * vi[invertible] + im[invertible]**2 * vr[invertible]
        - 2 * re[invertible] * im[invertible] * cri[invertible]
    ) / determinant[invertible]
    df[f"{stem}_origin_mahalanobis2"] = distance2
    df[f"Phase{stem}_displayed"] = nonzero & invertible & (distance2 > ELLIPSE_LEVEL)
    df[f"Phase{stem}_valid"] = df[f"Phase{stem}_displayed"]


def add_uu_studies(df, path):
    """Completa LT++ e LT-- e TT; nessuna separazione dei diagonali T/L."""
    missing = []
    c_u_l = optional_covariance(df, "AUU_cosPhi", "ALL_cosPhi", path, missing)
    c_u_a = optional_covariance(df, "AUU_cosPhi", "ALU_sinPhi", path, missing)
    c_u_b = optional_covariance(df, "AUU_cosPhi", "AUL_sinPhi", path, missing)
    c_l_a = optional_covariance(df, "ALL_cosPhi", "ALU_sinPhi", path, missing)
    c_tt = optional_covariance(df, "AUU_cos2Phi", "AUL_sin2Phi", path, missing)
    if missing:
        warnings.warn(f"{path.name}: covarianze aggiuntive assunte nulle: " + ", ".join(missing))
    df["uu_covariances_assumed_zero"] = ";".join(missing)
    u = df["AUU_cosPhi"].to_numpy()
    l = df["ALL_cosPhi"].to_numpy()
    vu = df["AUU_cosPhi_err"].to_numpy()**2
    vl = df["ALL_cosPhi_err"].to_numpy()**2
    c_l_b = df["cov_ALLcos_AUL_used"].to_numpy()
    for stem, sign in (("Gpp", 1), ("Gmm", -1)):
        df[f"Re{stem}_over_U0"] = -(u + sign*l) / np.sqrt(2)
        df[f"Re{stem}_over_U0_err"] = np.sqrt(np.maximum(0, (vu + vl + 2*sign*c_u_l)/2))
        # Cov(U +/- L, A +/- B)/2, con A=LU, B=UL.
        df[f"cov_Re_Im_{stem}_over_U0"] = (
            c_u_a + sign*c_u_b + sign*c_l_a + c_l_b) / 2
        add_polar_quantities(df, stem)
    df["cov_ReGpp_ReGmm_over_U0"] = (vu - vl) / 2
    df["ReG_TT_over_U0"] = -df["AUU_cos2Phi"]
    df["ReG_TT_over_U0_err"] = df["AUU_cos2Phi_err"]
    df["cov_Re_Im_G_TT_over_U0"] = c_tt
    add_polar_quantities(df, "G_TT")


def calculate(path, bin_column, x_column, allowed_bins, cov_column=None, rho=None,
              complex_cov_column=None, complex_rho=None):
    """Combina ampiezze della stessa riga/bin e propaga la covarianza."""
    df = pd.read_csv(path, skipinitialspace=True)
    df.columns = df.columns.str.strip()
    amplitudes = ["ALU_sinPhi", "AUL_sinPhi", "AUL_sin2Phi", "ALL", "ALL_cosPhi"]
    uu_columns = {"AUU_cosPhi", "AUU_cosPhi_err", "AUU_cos2Phi", "AUU_cos2Phi_err"}
    has_uu = bool(uu_columns & set(df.columns))
    if has_uu:
        if not uu_columns <= set(df.columns):
            raise ValueError(f"{path}: gruppo UU incompleto: {sorted(uu_columns - set(df.columns))}")
        amplitudes += ["AUU_cosPhi", "AUU_cos2Phi"]
    required = [bin_column, x_column] + amplitudes + [a + "_err" for a in amplitudes]
    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f"{path}: colonne mancanti: {sorted(missing)}")
    if df.empty:
        raise ValueError(f"{path}: nessuna riga.")
    for col in required:
        df[col] = pd.to_numeric(df[col], errors="raise")
    if not np.isfinite(df[required].to_numpy(dtype=float)).all():
        raise ValueError(f"{path}: valori mancanti o non finiti nelle colonne richieste.")
    if df[bin_column].duplicated().any():
        raise ValueError(f"{path}: identificatori di bin duplicati.")
    if not df[bin_column].isin(allowed_bins).all():
        raise ValueError(f"{path}: bin non contemplati nei gruppi dello script.")

    errors = df[[a + "_err" for a in amplitudes]].to_numpy()
    if (errors < 0).any():
        raise ValueError(f"{path}: errori negativi.")
    valid = (errors > 0).all(axis=1)
    if not valid.any():
        raise ValueError(
            f"{path}: nessun bin con tutti gli errori richiesti > 0. "
            "Verificare che il CSV contenga risultati del fit, non segnaposto."
        )
    if not valid.all():
        warnings.warn(f"{path}: esclusi bin con errore nullo: "
                      f"{df.loc[~valid, bin_column].tolist()}")
    df = df.loc[valid].copy()
    lu = df["ALU_sinPhi"].to_numpy()
    ul = df["AUL_sinPhi"].to_numpy()
    elu = df["ALU_sinPhi_err"].to_numpy()
    eul = df["AUL_sinPhi_err"].to_numpy()
    cov, source = covariance_values(df, cov_column, rho, "cov_ALU_AUL_sinPhi",
                                    "ALU_sinPhi", "AUL_sinPhi", path)

    df["cov_ALU_AUL_used"] = cov
    df["covariance_source"] = source
    df["ImGpp_over_U0"] = -(lu + ul) / np.sqrt(2)
    df["ImGmm_over_U0"] = -(lu - ul) / np.sqrt(2)
    df["ImGpp_over_U0_err"] = np.sqrt(np.maximum(0, (elu**2 + eul**2 + 2*cov)/2))
    df["ImGmm_over_U0_err"] = np.sqrt(np.maximum(0, (elu**2 + eul**2 - 2*cov)/2))
    # Le due quantita' derivate rimangono generalmente correlate.
    df["cov_ImGpp_ImGmm_over_U0"] = (elu**2 - eul**2) / 2

    df["DeltaG_T_over_U0"] = 2 * df["ALL"]
    df["DeltaG_T_over_U0_err"] = 2 * df["ALL_err"]
    df["ImG_TT_over_U0"] = -df["AUL_sin2Phi"]
    df["ImG_TT_over_U0_err"] = df["AUL_sin2Phi_err"]
    cov_lt, source_lt = covariance_values(
        df, complex_cov_column, complex_rho, "cov_ALL_cosPhi_AUL_sinPhi",
        "ALL_cosPhi", "AUL_sinPhi", path)
    df["cov_ALLcos_AUL_used"] = cov_lt
    df["complex_covariance_source"] = source_lt
    re = -np.sqrt(2) * df["ALL_cosPhi"].to_numpy()
    im = -np.sqrt(2) * ul
    vr = 2 * df["ALL_cosPhi_err"].to_numpy()**2
    vi = 2 * eul**2
    cri = 2 * cov_lt
    df["ReDeltaG_LT_over_U0"] = re
    df["ImDeltaG_LT_over_U0"] = im
    df["ReDeltaG_LT_over_U0_err"] = np.sqrt(vr)
    df["ImDeltaG_LT_over_U0_err"] = np.sqrt(vi)
    df["cov_Re_Im_DeltaG_LT_over_U0"] = cri
    add_polar_quantities(df, "DeltaG_LT")
    if has_uu:
        add_uu_studies(df, path)
    return df


def style(ax):
    ax.axhline(0, color="black", linestyle="--", linewidth=1, alpha=0.6)
    ax.tick_params(direction="in", top=True, right=True, length=5)
    ax.grid(alpha=0.25, linestyle=":")
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)


def configure_phase_display(df, mode):
    """Cambia soltanto i punti visualizzati; valori ed errori restano invariati."""
    if mode not in ("all", "resolved"):
        raise ValueError(f"Modalita' fase sconosciuta: {mode}")
    for stem in COMPLEX_FAMILIES:
        if f"Phase{stem}_deg" not in df:
            continue
        if mode == "all":
            df[f"Phase{stem}_displayed"] = (
                np.isfinite(df[f"Phase{stem}_deg"])
                & np.isfinite(df[f"Phase{stem}_deg_err"])
            )
        else:
            df[f"Phase{stem}_displayed"] = df[f"Phase{stem}_valid"]
    df["phase_display_mode"] = mode


def points(ax, df, x, observable, color, label):
    if observable.startswith("Phase") and observable.endswith("_deg"):
        df = df[df[observable.removesuffix("_deg") + "_displayed"]]
        if not df.empty:
            ax.plot(df[x], df[observable], "o", markersize=6, color=color,
                    markerfacecolor="white", label=label)
            # Intervalli circolari: disegna anche la porzione oltre il taglio +/-180.
            for _, row in df.iterrows():
                angle, err = row[observable], min(row[observable + "_err"], 180)
                for shift in (-360, 0, 360):
                    lo = max(-180, angle - err + shift)
                    hi = min(180, angle + err + shift)
                    if lo <= hi:
                        ax.vlines(row[x], lo, hi, color=color, linewidth=1.2)
        return
    if not df.empty:
        df = df.sort_values(x)
        ax.errorbar(df[x], df[observable], yerr=df[observable + "_err"],
                    fmt="o", markersize=6, capsize=3, linewidth=1.2,
                    color=color, markerfacecolor="white", label=label)


def save(fig, outdir, stem):
    fig.savefig(outdir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def phase_display_note(fig, df, mode, observable):
    shown = int(df[observable.removesuffix("_deg") + "_displayed"].sum())
    rule = (
        "All defined phases; local linear uncertainties with periodic intervals."
        if mode == "all" else
        "Display rule: joint 68.3% ellipse excludes the origin; not a fit-quality cut."
    )
    fig.supxlabel(
        f"Additional representation: {shown}/{len(df)} bins shown.\n"
        + rule,
        fontsize=8,
    )


def plot(zpt, xq2, outdir, phase_mode="resolved"):
    for observable, ylabel in OBSERVABLES:
        if observable not in zpt or observable not in xq2:
            continue
        is_phase = observable.startswith("Phase") and observable.endswith("_deg")
        suffix = "_all_bins" if is_phase and phase_mode == "all" else ""
        fig, axes = plt.subplots(1, 4, figsize=(18, 4), sharex=True, sharey=True,
                                 layout="constrained")
        for ax, (bins, title) in zip(axes, PT_GROUPS):
            selected = zpt[zpt["bin_zPt"].isin(bins)]
            points(ax, selected, "mean_z", observable, "tab:orange", "Fall22")
            if selected.empty:
                ax.text(0.5, 0.5, "Nessun bin valido", transform=ax.transAxes,
                        ha="center")
            style(ax)
            ax.set(title=title, xlabel=r"$z$", xlim=(0.2, 1.0))
            if is_phase:
                ax.set_ylim(-180, 180)
                ax.set_yticks([-180, -90, 0, 90, 180])
                if not selected[observable.removesuffix("_deg") + "_displayed"].any():
                    ax.text(.5, .7, "No phase shown\nwith this display rule",
                            transform=ax.transAxes, ha="center", fontsize=9)
        axes[0].set_ylabel(ylabel, fontsize=13)
        fig.suptitle(TITLE + "\n" + FORMULA_TITLES[observable], fontsize=12)
        if is_phase:
            phase_display_note(fig, zpt, phase_mode, observable)
        save(fig, outdir, f"{observable}_fall22_zPt{suffix}")

        fig, ax = plt.subplots(figsize=(7.5, 5.5), layout="constrained")
        for bins, color, label in Q2_GROUPS:
            points(ax, xq2[xq2["bin_xQ2"].isin(bins)], "mean_xB",
                   observable, color, label)
        style(ax)
        ax.set(xlabel=r"$x_B$", ylabel=ylabel)
        ax.set_title(TITLE + "\n" + FORMULA_TITLES[observable], fontsize=12)
        if is_phase:
            ax.set_ylim(-180, 180)
            ax.set_yticks([-180, -90, 0, 90, 180])
            phase_display_note(fig, xq2, phase_mode, observable)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(loc="best")
        save(fig, outdir, f"{observable}_fall22_xB_Q2{suffix}")
    plot_complex(zpt, xq2, outdir)


def complex_points(ax, df, bin_column, color, label=None, stem="DeltaG_LT"):
    ax.plot(df[f"Re{stem}_over_U0"], df[f"Im{stem}_over_U0"], "o",
            color=color, markerfacecolor="white", markersize=5, label=label)
    for _, row in df.iterrows():
        re, im = row[f"Re{stem}_over_U0"], row[f"Im{stem}_over_U0"]
        covariance = np.array([
            [row[f"Re{stem}_over_U0_err"]**2, row[f"cov_Re_Im_{stem}_over_U0"]],
            [row[f"cov_Re_Im_{stem}_over_U0"], row[f"Im{stem}_over_U0_err"]**2],
        ])
        eigenvalues, eigenvectors = np.linalg.eigh(covariance)
        width, height = 2 * np.sqrt(ELLIPSE_LEVEL * np.maximum(0, eigenvalues))
        angle = np.degrees(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0]))
        ax.add_patch(Ellipse((re, im), width, height, angle=angle,
                             facecolor=color, edgecolor=color, alpha=.18))
        ax.annotate(str(int(row[bin_column])), (re, im), xytext=(4, 4),
                    textcoords="offset points", fontsize=7, color=color)
    ax.autoscale_view()


def complex_style(ax, stem):
    style(ax)
    ax.axvline(0, color="black", linestyle="--", linewidth=1, alpha=.6)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel(dict(OBSERVABLES)[f"Re{stem}_over_U0"])


def plot_complex(zpt, xq2, outdir):
    for stem, title in COMPLEX_FAMILIES.items():
        if f"Re{stem}_over_U0" in zpt and f"Re{stem}_over_U0" in xq2:
            plot_complex_family(zpt, xq2, outdir, stem, title)


def plot_complex_family(zpt, xq2, outdir, stem, formula):
    fig, axes = plt.subplots(1, 4, figsize=(18, 5), layout="constrained")
    for ax, (bins, title) in zip(axes, PT_GROUPS):
        complex_points(ax, zpt[zpt["bin_zPt"].isin(bins)], "bin_zPt", "tab:orange", stem=stem)
        complex_style(ax, stem)
        ax.set_title(title, fontsize=11)
    axes[0].set_ylabel(dict(OBSERVABLES)[f"Im{stem}_over_U0"])
    fig.suptitle(TITLE + "\n" + formula, fontsize=12)
    fig.supxlabel("All input bins shown; joint 68.3% ellipses. Numbers identify bins.",
                  fontsize=9)
    save(fig, outdir, f"Complex{stem}_fall22_zPt")
    fig, ax = plt.subplots(figsize=(8, 7), layout="constrained")
    for bins, color, label in Q2_GROUPS:
        complex_points(ax, xq2[xq2["bin_xQ2"].isin(bins)], "bin_xQ2", color, label, stem)
    complex_style(ax, stem)
    ax.set_ylabel(dict(OBSERVABLES)[f"Im{stem}_over_U0"])
    ax.set_title(TITLE + "\n" + formula, fontsize=10)
    ax.legend(loc="best", fontsize=9)
    fig.supxlabel("All input bins shown; joint 68.3% ellipses. Numbers identify bins.",
                  fontsize=9)
    save(fig, outdir, f"Complex{stem}_fall22_xB_Q2")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-dir", type=Path, default=Path("RECO_CSV"))
    parser.add_argument("--outdir", type=Path, default=Path("helicity_sf_fall22"))
    parser.add_argument("--phase-mode", choices=("resolved", "all"), default="resolved",
                        help="resolved: ellisse esclude origine (default); all: tutte le fasi definite.")
    covariance = parser.add_mutually_exclusive_group()
    covariance.add_argument("--cov-column", help="Colonna Cov(ALU_sinPhi, AUL_sinPhi).")
    covariance.add_argument("--rho", type=float,
                            help="Forza correlazione ALU/AUL costante, invece della colonna.")
    complex_covariance = parser.add_mutually_exclusive_group()
    complex_covariance.add_argument("--complex-cov-column",
                                    help="Colonna Cov(ALL_cosPhi, AUL_sinPhi).")
    complex_covariance.add_argument("--complex-rho", type=float,
                                    help="Forza correlazione ALL_cosPhi/AUL_sinPhi costante.")
    args = parser.parse_args()
    for name in ("rho", "complex_rho"):
        value = getattr(args, name)
        if value is not None and (not np.isfinite(value) or not -1 <= value <= 1):
            parser.error(f"--{name.replace('_', '-')} deve essere compreso tra -1 e 1.")
    try:
        zpt = calculate(args.input_dir / "output_RGC_NH3_asymmetries_fall22.csv",
                        "bin_zPt", "mean_z", range(1, 26), args.cov_column, args.rho,
                        args.complex_cov_column, args.complex_rho)
        xq2 = calculate(args.input_dir / "output_RGC_NH3_asymmetries_fall22_xQ2.csv",
                        "bin_xQ2", "mean_xB", range(1, 15), args.cov_column, args.rho,
                        args.complex_cov_column, args.complex_rho)
    except (ValueError, OSError) as exc:
        parser.exit(1, f"Errore: {exc}\n")
    if ("ReGpp_over_U0" in zpt) != ("ReGpp_over_U0" in xq2):
        parser.exit(1, "Errore: fornire le colonne UU in entrambe le griglie, oppure in nessuna.\n")
    for df in (zpt, xq2):
        configure_phase_display(df, args.phase_mode)
    args.outdir.mkdir(parents=True, exist_ok=True)
    zpt.to_csv(args.outdir / "helicity_sf_fall22_zPt.csv", index=False)
    xq2.to_csv(args.outdir / "helicity_sf_fall22_xQ2.csv", index=False)
    plot(zpt, xq2, args.outdir, args.phase_mode)
    for name, df in (("zPt", zpt), ("xQ2", xq2)):
        reason = ("fase o errore non definito" if args.phase_mode == "all" else
                  "origine inclusa nell'ellisse al 68.3% o covarianza singolare")
        for stem in COMPLEX_FAMILIES:
            column = f"Phase{stem}_displayed"
            if column in df:
                omitted = len(df) - int(df[column].sum())
                print(f"{name}, {stem}: fase non mostrata in {omitted}/{len(df)} bin ({reason}).")
    count = 2 * sum(key in zpt and key in xq2 for key, _ in OBSERVABLES)
    count += 2 * sum(f"Re{stem}_over_U0" in zpt for stem in COMPLEX_FAMILIES)
    print(f"Creati {count} grafici PDF e 2 CSV in: {args.outdir.resolve()}")


if __name__ == "__main__":
    main()
