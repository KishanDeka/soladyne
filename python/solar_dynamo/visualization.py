"""Static and animated model diagnostics."""
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np
from .simulation import RunResult


FIELDS = {"Toroidal field": "toroidal", "Poloidal potential": "poloidal", "Rotation": "omega"}


def make_figure(result: RunResult, frame: int = -1, field: str = "Toroidal field"):
    if field not in FIELDS: raise ValueError(f"field must be one of {tuple(FIELDS)}")
    values = getattr(result, FIELDS[field])[frame]
    theta, radius = np.meshgrid(result.theta, result.radius)
    fig, ax = plt.subplots(figsize=(7, 5), subplot_kw={"projection": "polar"})
    contour = ax.contourf(theta, radius, values, levels=40, cmap="RdBu_r")
    ax.set_thetamin(0); ax.set_thetamax(180); ax.set_rorigin(0.35)
    ax.set_theta_zero_location("N"); ax.grid(alpha=.2)
    ax.set_title(f"{field} · t={result.time[frame]:.3e}")
    fig.colorbar(contour, ax=ax, pad=.08, shrink=.75)
    fig.tight_layout()
    return fig


def save_animation(result: RunResult, path: str | Path, field="Toroidal field", fps=12) -> Path:
    if field not in FIELDS: raise ValueError(f"field must be one of {tuple(FIELDS)}")
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    theta, radius = np.meshgrid(result.theta, result.radius)
    series = getattr(result, FIELDS[field]); bound = max(float(np.max(np.abs(series))), 1e-12)
    fig, ax = plt.subplots(figsize=(7, 5), subplot_kw={"projection": "polar"})
    def draw(i):
        ax.clear(); ax.contourf(theta, radius, series[i], levels=np.linspace(-bound, bound, 41), cmap="RdBu_r")
        ax.set_thetamin(0); ax.set_thetamax(180); ax.set_rorigin(0.35)
        ax.set_theta_zero_location("N"); ax.set_title(f"{field} · t={result.time[i]:.3e}")
        return ()
    animation = FuncAnimation(fig, draw, frames=len(result.time), interval=1000/fps)
    animation.save(path, writer=PillowWriter(fps=fps)); plt.close(fig)
    return path
