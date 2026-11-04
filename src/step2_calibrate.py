"""Step 2 - what the 8-bit numbers mean.

1. IR images: fit BT = offset + scale * DN through scenes of known temperature.
2. Histograms: the quantisation "combs" that show how the images were made.
3. file4: decode which colour is which channel.

Output: outputs/step2_calibration/calibration.json and three figures
"""

import json

import matplotlib.pyplot as plt
import numpy as np
from pyproj import Transformer

from src.common import ANCHORS, BT_UNCERTAINTY, FILES, box_values, out_dir, pretty, read


def fit_calibration(files):
    """Weighted least-squares line through the 'fit' anchors; residuals for all anchors."""
    rows = []
    for name, file, box, pct, expected, sigma, role in ANCHORS:
        dn = float(np.percentile(box_values(files[file], box), pct))
        rows.append({"name": name, "file": file, "dn": round(dn, 1), "expected": expected, "sigma": sigma, "role": role})
    fit = [r for r in rows if r["role"] == "fit"]
    dn = np.array([r["dn"] for r in fit])
    w = 1 / np.array([r["sigma"] for r in fit])
    A = np.vstack([np.ones_like(dn), dn]).T
    (offset, scale), *_ = np.linalg.lstsq(A * w[:, None], np.array([r["expected"] for r in fit]) * w, rcond=None)
    for r in rows:
        r["fitted"] = round(offset + scale * r["dn"], 1)
        r["residual"] = round(r["fitted"] - r["expected"], 1)
    return {"offset": round(float(offset), 3), "scale": round(float(scale), 5), "anchors": rows}


LABEL_OFFSETS = [(8, -4, "left"), (8, 6, "left"), (-8, -16, "right"), (8, 4, "left"), (8, -12, "left"), (-8, 10, "right")]


def plot_calibration(cal, path):
    fig, (ax, res) = plt.subplots(2, 1, figsize=(7.5, 6.5), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    d = np.linspace(0, 255, 2)
    ax.plot(d, cal["offset"] + cal["scale"] * d, "k-", lw=1, label="fitted line")
    ax.fill_between(d, cal["offset"] + cal["scale"] * d - BT_UNCERTAINTY, cal["offset"] + cal["scale"] * d + BT_UNCERTAINTY,
                    color="0.85", label=f"±{BT_UNCERTAINTY:.0f} K uncertainty")
    for i, r in enumerate(sorted(cal["anchors"], key=lambda r: r["dn"])):
        marker, color = ("o", "#0072B2") if r["role"] == "fit" else ("s", "#E69F00")
        ax.errorbar(r["dn"], r["expected"], yerr=r["sigma"], fmt=marker, color=color, capsize=3)
        res.errorbar(r["dn"], r["residual"], yerr=r["sigma"], fmt=marker, color=color, capsize=3)
        dx, dy, ha = LABEL_OFFSETS[i]  # hand-placed so the two close ocean anchors do not collide
        ax.annotate(r["name"], (r["dn"], r["expected"]), xytext=(dx, dy), textcoords="offset points", fontsize=7.5, ha=ha)
    ax.plot([], [], "o", color="#0072B2", label="anchor used for the fit")
    ax.plot([], [], "s", color="#E69F00", label="independent check")
    ax.set(ylabel="brightness temperature (K)",
           title=f"IR calibration recovered from known scenes: BT = {cal['offset']:.1f} − {-cal['scale']:.3f} · DN")
    ax.secondary_yaxis("right", functions=(lambda k: k - 273.15, lambda c: c + 273.15)).set_ylabel("°C")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.3)
    res.axhline(0, color="k", lw=0.8)
    res.set(xlabel="DN (8-bit pixel value; cold cloud = high DN)", ylabel="fit − expected (K)", xlim=(-5, 260))
    res.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_histograms(files, path):
    fig, axes = plt.subplots(2, 2, figsize=(12, 7))
    for ax, name in zip(axes.flat, FILES):
        f = files[name]
        labels = ["all pixels"] if len(f["bands"]) == 1 else ["R = VIS 0.6 µm", "G = NIR 0.8 µm", "B = inverted IR 10.8 µm"]
        colors = ["0.25"] if len(f["bands"]) == 1 else ["#D55E00", "#009E73", "#0072B2"]
        for band, color, label in zip(f["bands"], colors, labels):
            ax.plot(np.bincount(band[f["valid"]], minlength=256), color=color, lw=0.8, label=label)
        ax.set(yscale="log", xlabel="DN (8-bit value)", ylabel="number of pixels (log)", xlim=(0, 255),
               title=f"{name}: {pretty(FILES[name]['satellite'])}, {pretty(FILES[name]['channel'])}")
        if len(f["bands"]) > 1:
            ax.legend(fontsize=7.5, frameon=False)
    fig.suptitle("Histograms of valid pixels. The regular \"comb\" in file1 and file3: the instrument's radiance counts were\n"
                 "re-scaled linearly in temperature (not in radiance) onto 0-255, so some DN values receive more counts than others",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def solar_elevation(lon, lat, t):
    """Approximate solar elevation in degrees (good to ~1 deg)."""
    dec = np.radians(-23.44 * np.cos(np.radians(360 / 365 * (t.timetuple().tm_yday + 10))))
    hour_angle = np.radians((t.hour + t.minute / 60 - 12) * 15 + lon)
    lat = np.radians(lat)
    return np.degrees(np.arcsin(np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(hour_angle)))


def decode_rgb(f, path, step=4):
    """Show that R and G follow the sun (reflected light) and B does not (thermal IR)."""
    r, g, b = (band[::step, ::step].astype(float) for band in f["bands"])
    ok = f["valid"][::step, ::step]
    rows, cols = np.mgrid[0:f["valid"].shape[0]:step, 0:f["valid"].shape[1]:step]
    x, y = f["transform"] * (cols + 0.5, rows + 0.5)
    lon, lat = Transformer.from_crs(f["crs"], "EPSG:4326", always_xy=True).transform(x, y)
    elev = solar_elevation(lon, lat, f["time"])
    edges = np.arange(-4, 22, 3)
    means = {c: [v[ok & (elev >= lo) & (elev < hi)].mean() for lo, hi in zip(edges[:-1], edges[1:])]
             for c, v in zip("RGB", (r, g, b))}
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    names = {"R": "R: reflected sunlight (VIS 0.6 µm)", "G": "G: reflected sunlight (NIR 0.8 µm)",
             "B": "B: emitted heat (IR 10.8 µm, inverted)"}
    for c, color in zip("RGB", ["#D55E00", "#009E73", "#0072B2"]):
        ax.plot(edges[:-1] + 1.5, means[c], "o-", color=color, label=names[c])
    ax.axvline(0, color="0.5", lw=0.8, ls="--")
    ax.text(0.3, 5, "sunset line", fontsize=8, color="0.4")
    ax.set(xlabel="solar elevation at the pixel (°)", ylabel="mean DN in that band of elevation",
           title="file4: R and G dim with the sun, B does not\n→ R, G are solar channels, B is thermal")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    corr = np.corrcoef([r[ok], g[ok], b[ok]])
    return {"corr_R_G": round(corr[0, 1], 3), "corr_R_B": round(corr[0, 2], 3), "corr_G_B": round(corr[1, 2], 3)}


def main():
    out = out_dir("step2_calibration")
    files = {n: read(n) for n in FILES}
    cal = fit_calibration(files)
    plot_calibration(cal, out / "calibration_fit.png")
    plot_histograms(files, out / "histograms.png")
    cal["file4_channel_correlations"] = decode_rgb(files["file4"], out / "file4_channels_vs_sun.png")
    (out / "calibration.json").write_text(json.dumps(cal, indent=2))
    print(f"BT = {cal['offset']} {cal['scale']:+} * DN")
    for r in cal["anchors"]:
        print(f"  {r['role']:5s} {r['name']:40s} DN {r['dn']:6.1f}  expected {r['expected']}  residual {r['residual']:+}")
    print(cal["file4_channel_correlations"])


if __name__ == "__main__":
    main()
