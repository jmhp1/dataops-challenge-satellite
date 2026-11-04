"""Step 5 - merge everything into one documented NetCDF file (CF-1.10 / ACDD-1.3).

The file contains:
  - every input as its own variable on the common grid (IR in kelvin, others relative)
  - the IR mosaic: per pixel, the satellite that sees it most directly (smallest
    viewing zenith angle), plus which satellite, when, and at what angle

Arrays are opened lazily in chunks (dask), so the 6000 x 12000 grid never has to
fit in memory at once.

Output: outputs/step5_merged/satellite_mosaic.nc and composite.png
"""

import json
from datetime import datetime, timezone

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from src.common import BT_RANGE, BT_UNCERTAINTY, FILES, OUT, RES, VZA_LIMIT, ir_colormap, out_dir, wv_colormap

CHUNKS = {"lat": 600, "lon": 1200}
R_EARTH, H_GEO = 6378137.0, 35786023.0  # m: Earth's equatorial radius, geostationary height
PACK_K = {"dtype": "int16", "scale_factor": 0.01, "add_offset": 250.0, "_FillValue": -32768}  # 0.01 K steps
PACK_LEVEL = {"dtype": "int16", "scale_factor": 0.01, "add_offset": 0.0, "_FillValue": -32768}


def viewing_zenith(lat, lon, sub_lon):
    """Angle between the local vertical and the line to a geostationary satellite."""
    r_s = R_EARTH + H_GEO
    cos_g = np.cos(np.radians(lat)) * np.cos(np.radians(lon - sub_lon))
    d = np.sqrt(R_EARTH**2 + r_s**2 - 2 * R_EARTH * r_s * cos_g)
    vza = np.degrees(np.arcsin(np.clip(r_s * np.sqrt(1 - cos_g**2) / d, -1, 1)))
    return vza.where(cos_g > R_EARTH / r_s)  # NaN where the satellite is below the horizon


def build_product(cal):
    src = {n: xr.open_dataset(OUT / "step4_reprojected" / f"{n}.nc", chunks=CHUNKS) for n in FILES}
    lat, lon = src["file1"].lat.chunk(CHUNKS["lat"]), src["file1"].lon.chunk(CHUNKS["lon"])
    ds = xr.Dataset()
    encoding = {}

    # 1. every input as its own variable
    ir_inputs = []
    for n, info in FILES.items():
        common = {"platform": info["satellite"], "instrument": info["instrument"], "source_file": f"files/{n}",
                  "observation_time": src[n].attrs["time"], "grid_mapping": "crs"}
        for layer in src[n].data_vars:
            if info["kind"] == "ir":
                name = f"bt_{info['instrument'].lower()}"
                ds[name] = (cal["offset"] + cal["scale"] * src[n][layer]).assign_attrs(
                    standard_name="toa_brightness_temperature", units="K",
                    long_name=f"{info['satellite']} {info['instrument']} {info['channel']} brightness temperature",
                    comment=f"From 8-bit display levels: BT = {cal['offset']} + ({cal['scale']}) * DN, "
                            f"+/-{BT_UNCERTAINTY} K.", **common)
                encoding[name] = PACK_K
                ir_inputs.append((name, n))
            else:
                name = f"{info['instrument'].split('-')[0].lower()}_{layer}" if info["kind"] == "rgb" else "wv_abi"
                ds[name] = src[n][layer].assign_attrs(
                    units="1", long_name=f"{info['satellite']} {info['channel']}, layer {layer}, relative level 0-255",
                    comment="Uncalibrated display level, regridded by area average.", **common)
                encoding[name] = PACK_LEVEL

    # 2. the IR mosaic: keep, per pixel, the input with the smallest viewing zenith angle
    times = {n: datetime.fromisoformat(src[n].attrs["time"]) for _, n in ir_inputs}
    day = min(times.values()).replace(hour=0, minute=0)
    best = best_vza = source = minutes = None
    for i, (name, n) in enumerate(ir_inputs, start=1):
        vza = viewing_zenith(lat, lon, FILES[n]["sub_lon"]).where(ds[name].notnull())
        t = (times[n] - day).total_seconds() / 60
        if best is None:
            best, best_vza = ds[name], vza
            source, minutes = xr.where(ds[name].notnull(), i, 0), xr.where(ds[name].notnull(), t, np.nan)
        else:
            take = vza.fillna(np.inf) < best_vza.fillna(np.inf)
            best, best_vza = xr.where(take, ds[name], best), xr.where(take, vza, best_vza)
            source, minutes = xr.where(take, i, source), xr.where(take, t, minutes)
    flags = " ".join(f"{n}_{FILES[n]['instrument'].lower()}" for _, n in ir_inputs)
    ds["toa_brightness_temperature"] = xr.DataArray(best.data, dims=best.dims, attrs={
        "standard_name": "toa_brightness_temperature", "units": "K", "grid_mapping": "crs",
        "long_name": "IR brightness temperature mosaic: per pixel, the satellite with the smallest viewing angle",
        "ancillary_variables": "source_satellite sensor_zenith_angle observation_time"})
    ds["source_satellite"] = xr.DataArray(source.astype("uint8").data, dims=best.dims, attrs={
        "long_name": "input used for each mosaic pixel", "flag_values": np.arange(1, len(ir_inputs) + 1, dtype="uint8"),
        "flag_meanings": flags, "grid_mapping": "crs"})
    ds["sensor_zenith_angle"] = xr.DataArray(best_vza.data, dims=best.dims, attrs={
        "standard_name": "sensor_zenith_angle", "units": "degree", "grid_mapping": "crs",
        "comment": f"Above {VZA_LIMIT} deg the view is very oblique: use the temperatures with care."})
    ds["observation_time"] = xr.DataArray(minutes.data, dims=best.dims, attrs={
        "standard_name": "time", "units": f"minutes since {day:%Y-%m-%d %H:%M:%S}", "calendar": "standard",
        "long_name": "scan time of the satellite image used for each mosaic pixel", "grid_mapping": "crs"})
    encoding.update({"toa_brightness_temperature": PACK_K, "sensor_zenith_angle": PACK_LEVEL,
                     "source_satellite": {"dtype": "uint8", "_FillValue": 0},
                     "observation_time": {"dtype": "int16", "_FillValue": -1}})

    # 3. coordinates, grid mapping and file-level metadata
    ds = ds.assign_coords(lat=src["file1"].lat, lon=src["file1"].lon)
    ds["lat"].attrs = {"standard_name": "latitude", "long_name": "latitude", "units": "degrees_north", "axis": "Y"}
    ds["lon"].attrs = {"standard_name": "longitude", "long_name": "longitude", "units": "degrees_east", "axis": "X"}
    ds["crs"] = xr.DataArray(np.int32(0), attrs={"grid_mapping_name": "latitude_longitude",
                                                 "semi_major_axis": 6378137.0, "inverse_flattening": 298.257223563})
    known = sorted(src[n].attrs["time"] for n in FILES if src[n].attrs["time"] != "unknown")
    ds.attrs = {
        "Conventions": "CF-1.10, ACDD-1.3",
        "title": "Global geostationary IR mosaic and multi-sensor layers, 22 January 2025",
        "summary": "Four satellite display images on a common 0.03 degree lat/lon grid. IR temperatures are recovered "
                   f"from 8-bit levels (+/-{BT_UNCERTAINTY} K) and merged by smallest viewing angle.",
        "keywords": "satellite, geostationary, infrared, brightness temperature, water vapour, mosaic",
        "source": "StormGeo dataops-challenge-satellite files file1-file4",
        "processing_level": "L3 (regridded, merged)",
        "creator_name": "Joan Mateu Horrach",
        "time_coverage_start": known[0] + "Z", "time_coverage_end": known[-1] + "Z",
        "geospatial_lat_min": -90.0, "geospatial_lat_max": 90.0, "geospatial_lon_min": -180.0, "geospatial_lon_max": 180.0,
        "geospatial_lat_resolution": f"{RES} degree", "geospatial_lon_resolution": f"{RES} degree",
        "date_created": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "history": "\n".join([
            "step1: files identified (src/step1_describe.py)",
            f"step2: IR calibration BT = {cal['offset']} + ({cal['scale']}) * DN from physical anchors (src/step2_calibrate.py)",
            "step4: all files averaged onto EPSG:4326 0.03 deg; no-data masked; file3 CRS sweep=x (src/step4_reproject.py)",
            "step5: IR mosaic by smallest viewing zenith angle, written as CF/ACDD NetCDF (src/step5_merge.py)"]),
    }
    for v in encoding:
        encoding[v] = {**encoding[v], "zlib": True, "complevel": 4, "chunksizes": tuple(CHUNKS.values())}
    encoding.update({"lat": {"_FillValue": None}, "lon": {"_FillValue": None}})  # CF: coordinates have no fill
    return ds, encoding


def plot_composite(path, fig_path):
    """IR mosaic, water vapour only where IR is missing, file4 as an inset, seam in red."""
    with xr.open_dataset(path) as full:
        ds = full.isel(lat=slice(None, None, 4), lon=slice(None, None, 4)).load()
    pc, extent = ccrs.PlateCarree(), [-180, 180, -90, 90]
    fig = plt.figure(figsize=(16, 9.5))
    ax = fig.add_axes([0.02, 0.12, 0.96, 0.8], projection=ccrs.Robinson())
    ax.set_global()
    ax.add_feature(cfeature.LAND.with_scale("110m"), facecolor="0.9")
    bt = ds["toa_brightness_temperature"]
    im_wv = ax.imshow(ds["wv_abi"].where(bt.isnull()).values, extent=extent, transform=pc, origin="upper",
                      cmap=wv_colormap(), vmin=0, vmax=255)
    im_ir = ax.imshow(bt.values, extent=extent, transform=pc, origin="upper", cmap=ir_colormap(),
                      vmin=BT_RANGE[0], vmax=BT_RANGE[1])
    ax.add_feature(cfeature.COASTLINE.with_scale("110m"), linewidth=0.4, edgecolor="goldenrod")
    s = ds["source_satellite"].fillna(0).values
    rows, cols = np.nonzero((s[:, :-1] > 0) & (s[:, 1:] > 0) & (s[:, :-1] != s[:, 1:]))
    ax.scatter(ds.lon.values[cols], ds.lat.values[rows], s=0.3, color="crimson", transform=pc,
               label="seam: Meteosat | Himawari switch (equal viewing angle, ~79°)")
    ax.legend(loc="lower left", markerscale=10, fontsize=9)
    ax.text(0.99, 0.01, "light grey land: no geostationary IR data", transform=ax.transAxes, ha="right", fontsize=8,
            color="0.35")
    for pos, im, label in (([0.08, 0.07, 0.38, 0.02], im_ir, "IR mosaic: brightness temperature (K)"),
                           ([0.56, 0.07, 0.36, 0.02], im_wv, "gap fill (Americas): GOES water vapour, dry (dark) to moist (white)")):
        fig.colorbar(im, cax=fig.add_axes(pos), orientation="horizontal").set_label(label)
    rgb_vars = [v for v in ds.data_vars if v.startswith("avhrr_")]
    region = ds[rgb_vars].dropna("lat", how="all").dropna("lon", how="all")
    rgb = np.dstack([region[v].values / 255 for v in rgb_vars])
    ins = fig.add_axes([0.80, 0.70, 0.19, 0.22], projection=pc)
    ins.imshow(np.dstack([np.nan_to_num(rgb), np.isfinite(rgb[..., 0])]), transform=pc, origin="upper",
               extent=[region.lon[0], region.lon[-1], region.lat[-1], region.lat[0]])
    ins.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.3, edgecolor="white")
    ins.set_title("file4: polar-orbiter RGB over Scandinavia, 20 Jan 12:19 UTC", fontsize=8)
    ax.set_title("Global satellite mosaic, 22 January 2025: Meteosat IR (06:15 UTC) + Himawari IR (08:00 UTC), "
                 "GOES water vapour where no IR")
    fig.savefig(fig_path)
    plt.close(fig)


def main():
    cal = json.loads((OUT / "step2_calibration" / "calibration.json").read_text())
    out = out_dir("step5_merged")
    path = out / "satellite_mosaic.nc"
    ds, encoding = build_product(cal)
    ds.to_netcdf(path, encoding=encoding)

    with xr.open_dataset(path, chunks=CHUNKS) as check:  # reopen as a user would
        bt, w = check["toa_brightness_temperature"], np.cos(np.radians(check.lat))
        print(f"wrote {path.name}: {path.stat().st_size / 1e6:.0f} MB, {len(check.data_vars) - 1} variables")
        print(f"mosaic covers {float((bt.notnull() * w).sum() / (w.sum() * check.sizes['lon'])):.1%} of the globe, "
              f"BT {float(bt.min()):.0f}-{float(bt.max()):.0f} K")
    plot_composite(path, out / "composite.png")


if __name__ == "__main__":
    main()
