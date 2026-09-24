import pandas as pd
import matplotlib.pyplot as plt

# carica i file
nh3 = pd.read_csv("output_RGC_NH3_multi_fall22.csv", skipinitialspace=True)
nd3 = pd.read_csv("output_RGC_ND3_multi_fall22.csv", skipinitialspace=True)

titles_multi_xQ2bin = [
    r"$x_{B} < 0.12$  &  $1 < Q^{2} < 3$ GeV$^{2}$",
    r"$0.12 < x_{B} < 0.18$  &  $1 < Q^{2} < 3$ GeV$^{2}$",
    r"$0.18 < x_{B} < 0.24$  &  $1 < Q^{2} < 3$ GeV$^{2}$",
    r"$0.24 < x_{B} < 0.30$  &  $1 < Q^{2} < 3$ GeV$^{2}$",
    r"$x_{B} > 0.30$  &  $1 < Q^{2} < 3$ GeV$^{2}$",

    r"$x_{B} < 0.24$  &  $3 < Q^{2} < 5$ GeV$^{2}$",
    r"$0.24 < x_{B} < 0.30$  &  $3 < Q^{2} < 5$ GeV$^{2}$",
    r"$0.30 < x_{B} < 0.36$  &  $3 < Q^{2} < 5$ GeV$^{2}$",
    r"$0.36 < x_{B} < 0.44$  &  $3 < Q^{2} < 5$ GeV$^{2}$",
    r"$x_{B} > 0.44$  &  $3 < Q^{2} < 5$ GeV$^{2}$",

    r"$x_{B} < 0.44$  &  $5 < Q^{2} < 7$ GeV$^{2}$",
    r"$x_{B} > 0.44$  &  $5 < Q^{2} < 7$ GeV$^{2}$",

    r"$x_{B} < 0.55$  &  $7 < Q^{2} < 10$ GeV$^{2}$",
    r"$x_{B} > 0.55$  &  $7 < Q^{2} < 10$ GeV$^{2}$"
]

# loop sui bin di xB
for xb in sorted(nh3["bin_xB"].unique()):

    nh3_bin = nh3[nh3["bin_xB"] == xb]
    nd3_bin = nd3[nd3["bin_xB"] == xb]

    # ordina per z
    nh3_bin = nh3_bin.sort_values("mean_z")
    nd3_bin = nd3_bin.sort_values("mean_z")

    plt.figure()

    plt.errorbar(
        nh3_bin["mean_z"],
        nh3_bin["Multi"],
        yerr=nh3_bin["Multi_err"],
        fmt="o",
        label="NH3",
        color="cornflowerblue",
        capsize=3
    )

    plt.errorbar(
        nd3_bin["mean_z"],
        nd3_bin["Multi"],
        yerr=nd3_bin["Multi_err"],
        fmt="s",
        label="ND3",
        color="orange",
        capsize=3
    )

    plt.xlabel("z")
    plt.ylabel("$M^{K+}(x_{B}, Q^{2}, z, P_{hT})$")
    plt.title(titles_multi_xQ2bin[xb-1])
    plt.legend()
    plt.grid(True)

    plt.savefig(f"plots_multiplicities/multiplicity_xBbin_{xb}.png", dpi=300)

    plt.close()

print("Plots created in plots_multiplicities/ directory.")