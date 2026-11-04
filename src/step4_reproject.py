"""Step 4 - put every file on the same lat/lon grid (EPSG:4326, 0.03 deg).

- Resampling = average: every file is being made coarser, and averaging keeps
  fine detail from turning into noise (aliasing).
- "No data" pixels go in as NaN, so space never mixes into the edge of a disk.
- file2 runs past 180 E, so it is warped twice (once shifted by -360 deg).

Output: outputs/step4_reprojected/<file>.nc (values still raw DN) and an overview figure
"""

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from affine import Affine
from rasterio.enums import Resampling
from rasterio.warp import reproject

import json

from src.common import (BT_RANGE, FILES, GRID_LAT, GRID_LON, GRID_SHAPE, GRID_TRANSFORM, OUT, ir_colormap,
                        out_dir, pretty, read, wv_colormap)

LAYERS = {"ir": ["dn"], "wv": ["dn"], "rgb": ["red_vis06", "green_nir08", "blue_ir108_inverted"]}


def to_grid(f, band):
    """Average-resample one band of file `f` onto the common grid; NaN where there is no data."""
    src = np.where(f["valid"], band, np.nan).astype("float32")
    shifts = [0.0]
    if f["crs"].coordinate_operation.method_name == "Equidistant Cylindrical":
        half_world = np.pi * f["crs"].ellipsoid.semi_major_metre
        if f["transform"].c + f["valid"].shape[1] * f["transform"].a > half_world:  # past 180 E
            shifts.append(-2 * half_world)
    grid = np.full(GRID_SHAPE, np.nan, "float32")
    for dx in shifts:
        part = np.full(GRID_SHAPE, np.nan, "float32")
        reproject(src, part, src_transform=Affine.translation(dx, 0) * f["transform"], src_crs=f["crs"],
                  src_nodata=np.nan, dst_transform=GRID_TRANSFORM, dst_crs="EPSG:4326", dst_nodata=np.nan,
                  resampling=Resampling.average, num_threads=4)
        grid = np.where(np.isnan(grid), part, grid)
    return grid


def main():
    out = out_dir("step4_reprojected")
    cal = json.loads((OUT / "step2_calibration" / "calibration.json").read_text())  # for display only
    fig, axes = plt.subplots(2, 2, figsize=(16, 9), subplot_kw={"projection": ccrs.PlateCarree()})
    for ax, name in zip(axes.flat, FILES):
        f = read(name)
        info = FILES[name]
        ds = xr.Dataset({layer: (("lat", "lon"), to_grid(f, band)) for layer, band in zip(LAYERS[info["kind"]], f["bands"])},
                        coords={"lat": GRID_LAT, "lon": GRID_LON},
                        attrs={**{k: v for k, v in info.items() if isinstance(v, (str, float))},
                               "time": f["time"].isoformat() if f["time"] else "unknown"})
        ds.to_netcdf(out / f"{name}.nc", encoding={v: {"zlib": True, "complevel": 1} for v in ds.data_vars})
        cover = float(np.isfinite(ds[LAYERS[info["kind"]][0]]).mean())
        print(f"{name}: covers {cover:.1%} of the grid")

        thin = ds.isel(lat=slice(None, None, 4), lon=slice(None, None, 4))
        if info["kind"] == "rgb":
            rgb = np.dstack([thin[v].values / 255 for v in LAYERS["rgb"]])
            img = np.dstack([np.nan_to_num(rgb), np.isfinite(rgb[..., 0]).astype(float)])
            ax.imshow(img, extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree(), origin="upper")
        elif info["kind"] == "ir":  # shown as brightness temperature, same colours as the step 3 maps
            im_ir = ax.imshow(cal["offset"] + cal["scale"] * thin["dn"].values, extent=[-180, 180, -90, 90],
                              transform=ccrs.PlateCarree(), origin="upper", cmap=ir_colormap(),
                              vmin=BT_RANGE[0], vmax=BT_RANGE[1], interpolation="nearest")
        else:
            ax.imshow(thin["dn"].values, extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree(), origin="upper",
                      cmap=wv_colormap(), vmin=0, vmax=255, interpolation="nearest")
        ax.add_feature(cfeature.COASTLINE.with_scale("110m"), linewidth=0.4, edgecolor="goldenrod")
        ax.set_global()
        ax.set_title(f"{name}: {pretty(info['satellite'])}, {pretty(info['channel'])}, covers {cover:.0%} of the grid",
                     fontsize=9)
    fig.suptitle("All four files on the common grid (lat/lon, EPSG:4326, 0.03° ≈ 3.3 km)", fontweight="bold")
    fig.tight_layout(rect=(0, 0.06, 1, 0.95))
    fig.colorbar(im_ir, cax=fig.add_axes([0.3, 0.04, 0.4, 0.015]), orientation="horizontal",
                 label="IR brightness temperature (K), files 1 and 2")
    fig.savefig(out / "all_files_on_grid.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
