from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import randompack
import varmapack

OUTPUT_DIR = Path(__file__).parents[1] / "docs" / "figures"
FONT_SIZE = 10
CROSS_FONT_SIZE = 15

def compute_results():
    n = 80
    nrep = 5
    A1 = np.array([[0.75, 0.05], [0, 0.5]])
    A2 = np.array([[0.13, 0], [0, 0.05]])
    B1 = np.array([[0.4, 0.15], [0.05, 0.2]])
    Sig = np.array([[1, 0.99], [0.99, 1]])
    model = varmapack.Model(A=[A1, A2], B=B1, Sig=Sig)

    rng = randompack.Rng()
    rng.seed(123)
    startup = np.linspace(0, 5, 11)
    initial_values = np.column_stack((startup, startup))
    paths = model.sim(n, nrep=nrep, X0=initial_values, rng=rng)

    Gamma = model.acvf(16)
    scale = np.sqrt(Gamma[0, 0, 0]*Gamma[0, 1, 1])
    lags = np.arange(-4, 17)
    indices = np.abs(lags)
    rho_xy = np.where(
        lags >= 0, Gamma[indices, 0, 1], Gamma[indices, 1, 0]
    )/scale
    return paths, len(initial_values), lags, rho_xy

def plot_cross_correlation(lags, rho_xy):
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.plot(lags, rho_xy, color="tab:blue", linewidth=1.5)
    ax.axhline(0, color="0.3", linewidth=0.8)
    ax.set_xlabel(r"Lag $k$ (positive when $y$ leads $x$)",
                  fontsize=CROSS_FONT_SIZE)
    ax.set_ylabel(r"Cross-correlation $\rho_{xy}(k)$", fontsize=CROSS_FONT_SIZE)
    ax.set_xlim(-4, 16)
    ax.set_ylim(0, 0.8)
    ax.set_xticks(np.arange(-4, 17, 2))
    ax.tick_params(axis="both", labelsize=CROSS_FONT_SIZE)
    ax.grid(axis="y", color="0.88", linewidth=0.6)
    fig.tight_layout()
    output = OUTPUT_DIR / "bivariate_varma_cross_correlation"
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output.with_suffix(".svg"), bbox_inches="tight")

def plot_paths(paths, h):
    nrep, n, _ = paths.shape
    t = np.arange(n)
    colors = plt.get_cmap("tab10").colors[:nrep]
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.2), sharex=True)
    labels = (r"$x_t$ (asymptotic decay rate 0.90)",
              r"$y_t$ (asymptotic decay rate 0.59)")
    for component, ax in enumerate(axes):
        ax.axvspan(0, h - 1, color="0.92")
        ax.axvline(h - 1, color="0.35", linestyle="--", linewidth=1)
        for replicate, color in enumerate(colors):
            ax.plot(t, paths[replicate, :, component], color=color,
                    alpha=0.8, linewidth=1.5)
        ax.set_ylabel(labels[component], fontsize=FONT_SIZE)
        ax.set_ylim(-8, 8)
        ax.set_yticks(np.arange(-8, 9, 2))
        ax.tick_params(axis="both", labelsize=FONT_SIZE)
        ax.grid(axis="both", color="0.88", linewidth=0.6)
    axes[-1].set_xlabel("Time $t$", fontsize=FONT_SIZE)
    axes[-1].set_xlim(0, n)
    axes[-1].set_xticks(np.arange(0, n + 1, 10))
    fig.tight_layout()
    output = OUTPUT_DIR / "bivariate_varma_paths"
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output.with_suffix(".svg"), bbox_inches="tight")

if __name__ == "__main__":
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    paths, h, lags, rho_xy = compute_results()
    plot_cross_correlation(lags, rho_xy)
    plot_paths(paths, h)
