"""Step 1 - describe the four files: format, projection, size, time, coverage.

Output: outputs/step1_files/files.md
"""

import numpy as np
import rasterio
from pyproj import Transformer

from src.common import FILES, DATA, out_dir, read


def subsatellite_lon(f):
    """Longitude a geostationary disk is centred on, from the symmetry of its edge.

    Only for disks resampled to lat/lon (eqc): on rows at +-60/70 deg the first and
    last valid pixels sit symmetrically around the satellite's longitude.
    """
    if f["crs"].coordinate_operation.method_name != "Equidistant Cylindrical":
        return None
    m_per_deg = np.pi * f["crs"].ellipsoid.semi_major_metre / 180
    mids = []
    for lat in (70, 60, -60, -70):
        row = rasterio.transform.rowcol(f["transform"], 0, lat * m_per_deg)[0]
        cols = np.nonzero(f["valid"][row])[0]
        if 0 < cols[0] and cols[-1] < f["valid"].shape[1] - 1:  # edge not cut by the image border
            x_west = (f["transform"] * (cols[0] + 0.5, row))[0]
            x_east = (f["transform"] * (cols[-1] + 0.5, row))[0]
            mids.append((x_west + x_east) / 2 / m_per_deg)
    return round(float(np.median(mids)), 2)


def describe(name):
    f = read(name)
    with rasterio.open(DATA / name) as src:
        size, dtype, compression = (src.width, src.height), src.dtypes[0], src.compression
        declared_nodata, res = src.nodata, src.res
    rows, cols = np.nonzero(f["valid"])
    cx, cy = f["transform"] * ((cols.min() + cols.max()) / 2, (rows.min() + rows.max()) / 2)
    lon, lat = Transformer.from_crs(f["crs"], "EPSG:4326", always_xy=True).transform(cx, cy)
    return {**FILES[name], "file": name, "size_px": size, "bands": f["bands"].shape[0], "dtype": dtype,
            "compression": str(compression), "projection": f["crs"].coordinate_operation.method_name,
            "pixel_size_m": [round(r, 1) for r in res], "declared_nodata": declared_nodata,
            "time": f["time"].isoformat() if f["time"] else None,
            "valid_fraction": round(float(f["valid"].mean()), 3),
            "centre_lon_lat": [round(lon, 1), round(lat, 1)], "subsatellite_lon": subsatellite_lon(f)}


def main():
    rows = [describe(n) for n in FILES]
    out = out_dir("step1_files")
    keys = ["satellite", "instrument", "channel", "size_px", "bands", "projection", "pixel_size_m",
            "time", "declared_nodata", "valid_fraction", "centre_lon_lat", "subsatellite_lon"]
    table = ["| | " + " | ".join(r["file"] for r in rows) + " |", "|---" * (len(rows) + 1) + "|"]
    table += [f"| {k} | " + " | ".join(str(r[k]) for r in rows) + " |" for k in keys]
    (out / "files.md").write_text("# Step 1: the four files\n\n" + "\n".join(table) + "\n")
    print("\n".join(table))


if __name__ == "__main__":
    main()
