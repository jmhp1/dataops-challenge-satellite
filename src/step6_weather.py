"""Step 6 - weather insights.

Everything is read from the step 5 product, except the open-cell measurement,
which needs file4's native ~1 km pixels (the 0.03 deg grid is too coarse for
~20 km cells). Temperature thresholds are shifted by +/-5 K to show how robust
each threshold-based result is.

Output: outputs/step6_weather/weather.json and four figures
"""

import json

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from pyproj import Transformer
from rasterio.features import rasterize
from scipy import ndimage

from matplotlib.ticker import ScalarFormatter

from src.common import BT_RANGE, BT_UNCERTAINTY, GRID_TRANSFORM, OUT, RES, VZA_LIMIT, ir_colormap, out_dir, read

KM_PER_DEG = 111.32
SHIELD, CORE, COLD_CLOUD = 241.0, 221.0, 235.0   # K: Maddox (1980) MCC thresholds; ice-topped cloud
MCC_SHIELD_KM2, MCC_CORE_KM2 = 100_000, 50_000
ITCZ_SECTORS = {"Atlantic 40W-10W": (-40, -10), "Africa 10E-40E": (10, 40), "Indian Ocean 55E-90E": (55, 90),
                "Maritime Continent 100E-150E": (100, 150), "West Pacific 150E-180": (150, 180)}
LAND_SEA = {  # land box and nearby ocean box: (west, east, south, north)
    "Sahara (07:00 local)": ((0, 25, 20, 30), (-25, -15, 20, 30)),
    "Arabia (09:30 local)": ((42, 55, 18, 28), (60, 68, 10, 20)),
    "Southern Africa (07:45 local)": ((15, 30, -30, -20), (0, 10, -30, -20)),
    "Australia (17:00 local)": ((120, 145, -30, -20), (95, 110, -30, -20)),
}
OPEN_CELLS_BOX = (-16, -9, 61, 63)
# Remnant of Severe Tropical Cyclone Sean at 06 UTC 22 Jan 2025 (Bureau of Meteorology bulletin)
EX_SEAN = {"lon": 107.1, "lat": -27.8, "box": (90, 130, -45, -10)}
COLORS = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7"]  # Okabe-Ito, colour-blind safe


def rows_cols(lat, lon, box):
    w, e, s, n = box
    return np.ix_((lat >= s) & (lat <= n), (lon >= w) & (lon <= e))


def land_sea(lat, lon, bt):
    """Warmest 10 % (clearest) of land vs nearby ocean: the diurnal cycle as a check."""
    land = rasterize(((g, 1) for g in cfeature.NaturalEarthFeature("physical", "land", "50m").geometries()),
                     out_shape=bt.shape, transform=GRID_TRANSFORM, fill=0, dtype="uint8").astype(bool)
    out = {}
    for name, (land_box, ocean_box) in LAND_SEA.items():
        t_land = float(np.nanpercentile(bt[rows_cols(lat, lon, land_box)][land[rows_cols(lat, lon, land_box)]], 90))
        t_sea = float(np.nanpercentile(bt[rows_cols(lat, lon, ocean_box)][~land[rows_cols(lat, lon, ocean_box)]], 90))
        out[name] = {"land_K": round(t_land, 1), "ocean_K": round(t_sea, 1), "land_minus_ocean_K": round(t_land - t_sea, 1)}
    return out


def itcz(lat, lon, bt, path):
    """Latitude where cold (ice-topped) cloud is most frequent, per longitude sector."""
    tropics = np.abs(lat) <= 30
    out = {}
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for (name, (w, e)), color in zip(ITCZ_SECTORS.items(), COLORS):
        b = bt[np.ix_(tropics, (lon >= w) & (lon <= e))]
        n_valid = np.maximum(np.isfinite(b).sum(1), 1)
        peaks = {}
        for shift in (-BT_UNCERTAINTY, 0, BT_UNCERTAINTY):
            frac = ndimage.uniform_filter1d((b < COLD_CLOUD + shift).sum(1) / n_valid, size=round(2 / RES))  # 2 deg mean
            peaks[shift] = round(float(lat[tropics][np.argmax(frac)]), 1)
            if shift == 0:
                ax.plot(lat[tropics], frac, color=color, lw=2, label=f"{name}: peak at {abs(peaks[0]):.1f}°{'N' if peaks[0] >= 0 else 'S'}")
        out[name] = {"peak_lat": peaks[0], "peak_lat_at_-5K": peaks[-BT_UNCERTAINTY], "peak_lat_at_+5K": peaks[BT_UNCERTAINTY]}
    ax.axvline(0, color="0.6", lw=0.8)
    ax.text(0.3, ax.get_ylim()[1] * 0.97, "equator", fontsize=8, color="0.4", va="top")
    ax.set(xlim=(-30, 30), xlabel="latitude (°, negative = south)",
           ylabel=f"share of pixels colder than {COLD_CLOUD:.0f} K (−38 °C)",
           title="Where is the deep tropical convection? Ice-topped cloud by latitude, per longitude sector\n"
                 "(22 Jan 2025, 2° running mean. In January the belt shifts south over land and the Indian Ocean)")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return out


def cold_shields(lat, bt, vza, shift=0.0):
    """Connected areas colder than 241 K, with their core area colder than 221 K."""
    cold = np.nan_to_num(bt, nan=999.0)
    labels, n = ndimage.label(cold <= SHIELD + shift, structure=np.ones((3, 3)))
    ids = np.arange(1, n + 1)
    area = np.broadcast_to(((RES * KM_PER_DEG) ** 2 * np.cos(np.radians(lat)))[:, None], bt.shape)
    shield = ndimage.sum(area, labels, ids)
    core = ndimage.sum(area * (cold <= CORE + shift), labels, ids)
    big = ids[shield >= MCC_SHIELD_KM2 / 5]
    centres = ndimage.center_of_mass(np.ones(bt.shape), labels, big)
    mean_vza = ndimage.mean(np.nan_to_num(vza, nan=90.0), labels, big)
    min_bt = ndimage.minimum(cold, labels, big)
    systems = []
    for i, (r, c), v, t in zip(big, centres, mean_vza, min_bt):
        clat = float(lat[int(r)])
        systems.append({"id": int(i), "row_col": [int(r), int(c)], "lat": round(clat, 1),
                        "shield_km2": int(shield[i - 1]), "core_km2": int(core[i - 1]), "min_bt_K": round(float(t), 1),
                        "mean_vza": round(float(v), 1),
                        "trustworthy": bool(abs(clat) <= 25 and v <= VZA_LIMIT),  # tropical and well viewed
                        "mcc_size": bool(shield[i - 1] >= MCC_SHIELD_KM2 and core[i - 1] >= MCC_CORE_KM2)})
    return sorted(systems, key=lambda s: -s["shield_km2"]), labels


def plot_shields(lat, lon, bt, systems, labels, path, top=10):
    k = 4
    fig = plt.figure(figsize=(15, 8))
    ax = fig.add_subplot(projection=ccrs.Robinson())
    ax.set_global()
    im = ax.imshow(bt[::k, ::k], extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree(), origin="upper",
                   cmap=ir_colormap(), vmin=BT_RANGE[0], vmax=BT_RANGE[1])
    fig.colorbar(im, ax=ax, orientation="horizontal", shrink=0.5, pad=0.03, label="brightness temperature (K)")
    ax.add_feature(cfeature.COASTLINE.with_scale("110m"), linewidth=0.4, edgecolor="goldenrod")
    good = [s for s in systems if s["trustworthy"]][:top]
    false = [s for s in systems[:12] if s["mean_vza"] > 80]
    for group, color, style in ((good, COLORS[3], "-"), (false, COLORS[0], "--")):
        for s in group:
            ax.contour(lon[::k], lat[::k], (labels[::k, ::k] == s["id"]).astype(float), levels=[0.5],
                       colors=[color], linewidths=1.2, linestyles=style, transform=ccrs.PlateCarree())
    for rank, s in enumerate(good, 1):
        ax.text(float(lon[s["row_col"][1]]), s["lat"], str(rank), transform=ccrs.PlateCarree(), fontsize=9,
                fontweight="bold", ha="center", va="center", bbox=dict(boxstyle="circle,pad=0.2", fc="white", ec=COLORS[3]))
    ax.plot([], [], color=COLORS[3], label=f"largest tropical shields, ranked by area (|lat| ≤ 25°, viewing angle ≤ {VZA_LIMIT:.0f}°)")
    ax.plot([], [], color=COLORS[0], ls="--", label="false 'shields' seen at viewing angle > 80°: cold ground / limb cooling, not storms")
    ax.legend(loc="lower left", fontsize=9)
    ax.set_title(f"Cold cloud shields: connected areas colder than {SHIELD:.0f} K (−32 °C), 22 Jan 2025")
    fig.savefig(path)
    plt.close(fig)


def open_cells(path):
    """Typical size of the open cells in file4, from the 2-D power spectrum of the visible channel."""
    f = read("file4")
    w, e, s, n = OPEN_CELLS_BOX
    xs, ys = Transformer.from_crs("EPSG:4326", f["crs"], always_xy=True).transform([w, e, w, e], [s, s, n, n])
    inv = ~f["transform"]
    cols, rows = zip(*(inv * (x, y) for x, y in zip(xs, ys)))
    vis = f["bands"][0][int(min(rows)):int(max(rows)), int(min(cols)):int(max(cols))].astype(float)
    dx, dy = f["transform"].a / 1000, -f["transform"].e / 1000  # km per pixel (1.15 x 1.0)
    v = (vis - vis.mean()) * np.outer(np.hanning(vis.shape[0]), np.hanning(vis.shape[1]))  # taper the edges
    power = np.abs(np.fft.fftshift(np.fft.fft2(v))) ** 2
    ky, kx = np.fft.fftshift(np.fft.fftfreq(v.shape[0], dy)), np.fft.fftshift(np.fft.fftfreq(v.shape[1], dx))
    k = np.hypot(*np.meshgrid(kx, ky)).ravel()     # cycles per km
    bins = np.linspace(1 / 200, 1 / 4, 120)        # wavelengths 4-200 km
    which = np.digitize(k, bins)
    spectrum = np.array([power.ravel()[which == i].mean() for i in range(1, len(bins))])
    k_mid = (bins[1:] + bins[:-1]) / 2
    peak = float(1 / k_mid[np.argmax(k_mid * spectrum)])  # k * P(k) removes the background slope
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5))
    h, w_px = vis.shape
    a1.imshow(vis, cmap="gray", extent=[0, w_px * dx, 0, h * dy])
    a1.set(title=f"file4 VIS 0.6 µm, {s}-{n}°N, {-e}-{-w}°W (between Iceland and the Faroes)\n"
                 "open cells: dark centres, bright cloud rims", xlabel="km", ylabel="km")
    a2.plot(1 / k_mid, k_mid * spectrum / np.max(k_mid * spectrum), color=COLORS[0], lw=2)
    a2.axvline(peak, color=COLORS[3], ls="--")
    a2.set(xscale="log", xlabel="wavelength (km)", ylabel="k · P(k), normalised (energy per log-wavelength)",
           title=f"Spatial power spectrum: dominant cell size ≈ {peak:.0f} km")
    a2.set_xticks([5, 10, 20, 50, 100, 200])
    a2.xaxis.set_major_formatter(ScalarFormatter())
    a2.grid(alpha=0.25, which="both")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return {"box": OPEN_CELLS_BOX, "pixel_km": [round(dx, 2), round(dy, 2)], "cell_size_km": round(peak, 1)}


def plot_ex_sean(lat, lon, bt, path):
    """Himawari view around the last reported position of ex-Tropical Cyclone Sean."""
    w, e, s, n = EX_SEAN["box"]
    rows, cols = rows_cols(lat, lon, EX_SEAN["box"])
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(projection=ccrs.PlateCarree())
    im = ax.imshow(bt[rows, cols], extent=[w, e, s, n], transform=ccrs.PlateCarree(), origin="upper",
                   cmap=ir_colormap(), vmin=BT_RANGE[0], vmax=BT_RANGE[1])
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.6, edgecolor="goldenrod")
    ax.plot(EX_SEAN["lon"], EX_SEAN["lat"], "+", color=COLORS[3], ms=16, mew=2.5,
            label="ex-TC Sean, BoM position 06 UTC")
    ax.legend(loc="upper left", fontsize=8)
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color="0.5", linestyle=":")
    gl.top_labels = gl.right_labels = False
    fig.colorbar(im, ax=ax, orientation="horizontal", shrink=0.7, pad=0.08, label="brightness temperature (K)")
    ax.set_title("Remnant of Cyclone Sean, Himawari IR, 22 Jan 2025 08:00 UTC\n"
                 "low-level swirl without deep convection, cloud band joining a front to the south", fontsize=10)
    fig.savefig(path)
    plt.close(fig)


def main():
    out = out_dir("step6_weather")
    with xr.open_dataset(OUT / "step5_merged" / "satellite_mosaic.nc", chunks={"lat": 600, "lon": 1200}) as ds:
        lat, lon = ds.lat.values, ds.lon.values
        bt = ds["toa_brightness_temperature"].astype("float32").values
        vza = ds["sensor_zenith_angle"].astype("float32").values
    systems, labels = cold_shields(lat, bt, vza)
    counts = {f"{s:+.0f} K": int(sum(x["mcc_size"] and x["trustworthy"] for x in cold_shields(lat, bt, vza, s)[0]))
              for s in (-BT_UNCERTAINTY, BT_UNCERTAINTY)}
    counts["+0 K"] = int(sum(x["mcc_size"] and x["trustworthy"] for x in systems))
    plot_shields(lat, lon, bt, systems, labels, out / "cold_shields.png")
    plot_ex_sean(lat, lon, bt, out / "ex_sean.png")
    result = {"land_minus_ocean": land_sea(lat, lon, bt), "itcz": itcz(lat, lon, bt, out / "itcz.png"),
              "mcc_size_tropical_systems": counts,
              "largest_tropical_systems": [s for s in systems if s["trustworthy"]][:10],
              "largest_shields_anywhere": systems[:3], "open_cells": open_cells(out / "open_cells.png")}
    (out / "weather.json").write_text(json.dumps(result, indent=2))
    print(json.dumps({k: result[k] for k in ("land_minus_ocean", "itcz", "mcc_size_tropical_systems", "open_cells")}, indent=1))


if __name__ == "__main__":
    main()
