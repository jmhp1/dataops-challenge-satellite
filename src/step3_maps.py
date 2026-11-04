"""Step 3 - a map of each file in its own projection (nothing resampled yet).

IR images are shown as brightness temperature (step 2 calibration), the water
vapour image as relative levels, file4 as delivered (RGB).

Output: outputs/step3_maps/map_file1.png ... map_file4.png
"""

import json

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np

from src.common import BT_RANGE, FILES, OUT, ir_colormap, map_gridlines, out_dir, pretty, read, wv_colormap


def map_projection(f):
    """Cartopy projection and image extent (left, right, bottom, top) for one file."""
    crs, t = f["crs"], f["transform"]
    h, w = f["valid"].shape
    extent = [t.c, t.c + w * t.a, t.f + h * t.e, t.f]
    p = crs.to_dict()
    globe = ccrs.Globe(ellipse=p.get("ellps"), datum=p.get("datum"))
    method = crs.coordinate_operation.method_name
    if method == "Equidistant Cylindrical":  # lat/lon stored in metres: convert to degrees around the data
        m_per_deg = np.pi * crs.ellipsoid.semi_major_metre / 180
        centre = round((extent[0] + extent[1]) / 2 / m_per_deg)
        return ccrs.PlateCarree(central_longitude=centre), [extent[0] / m_per_deg - centre, extent[1] / m_per_deg - centre,
                                                            extent[2] / m_per_deg, extent[3] / m_per_deg]
    if method.startswith("Geostationary"):
        return ccrs.Geostationary(central_longitude=p["lon_0"], satellite_height=p["h"],
                                  sweep_axis=p.get("sweep", "y"), globe=globe), extent
    return ccrs.Stereographic(central_latitude=p["lat_0"], central_longitude=p["lon_0"],
                              true_scale_latitude=p["lat_ts"], globe=globe), extent


def draw(name, cal, path, max_px=3000):
    f = read(name)
    k = max(1, max(f["valid"].shape) // max_px)  # thin out for display
    valid = f["valid"][::k, ::k]
    proj, extent = map_projection(f)

    kind = FILES[name]["kind"]
    if kind == "rgb":
        img = np.dstack([b[::k, ::k] / 255 for b in f["bands"]] + [valid.astype(float)])
        style, label = {}, None
    else:
        img = np.where(valid, f["bands"][0][::k, ::k].astype(float), np.nan)
        if kind == "ir":
            img = cal["offset"] + cal["scale"] * img
            style, label = dict(cmap=ir_colormap(), vmin=BT_RANGE[0], vmax=BT_RANGE[1]), "brightness temperature (K)"
        else:
            style, label = dict(cmap=wv_colormap(), vmin=0, vmax=255), \
                "uncalibrated level (DN): dark = dry upper troposphere, white = moist air / high cloud"

    fig = plt.figure(figsize=(10, 9))
    ax = fig.add_subplot(projection=proj)
    im = ax.imshow(img, extent=extent, transform=proj, origin="upper", interpolation="nearest", **style)
    if not isinstance(proj, ccrs.Geostationary):  # a full disk is the default view (its corners are in space)
        ax.set_extent(extent, crs=proj)
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.4, edgecolor="goldenrod")
    map_gridlines(ax)
    if label:
        fig.colorbar(im, ax=ax, orientation="horizontal", shrink=0.7, pad=0.04, label=label)
    info = FILES[name]
    when = f["time"].strftime("%Y-%m-%d %H:%M UTC") if f["time"] else "time unknown"
    ax.set_title(f"{name}: {pretty(info['satellite'])} {info['instrument']}, {pretty(info['channel'])}\n"
                 f"{when}, shown in the file's own projection", fontsize=11)
    if kind == "ir":
        fig.text(0.5, 0.01, "Grey: warm surfaces and low cloud.  Purple: cloud tops colder than 235 K (−38 °C), "
                 "ice clouds, mostly deep convection or frontal cloud.", ha="center", fontsize=8, color="0.3")
    if kind == "rgb":
        fig.text(0.5, 0.01, "White: thick, cold cloud.  Yellow: low, warm cloud (bright but not cold).  Green: land.  "
                 "Dark: open sea.\nBlue: no sunlight left (sun near or below the horizon), only the IR channel shows.",
                 ha="center", fontsize=8, color="0.3")
    fig.savefig(path)
    plt.close(fig)


def main():
    cal = json.loads((OUT / "step2_calibration" / "calibration.json").read_text())
    out = out_dir("step3_maps")
    for name in FILES:
        draw(name, cal, out / f"map_{name}.png")
        print("drew", name)


if __name__ == "__main__":
    main()
