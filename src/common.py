"""Settings and helpers shared by all steps.

Every decision (what each file is, the grid, the thresholds) is written here
once, so the step scripts only contain processing.
"""

from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold", "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.dpi": 130, "savefig.bbox": "tight", "figure.facecolor": "white"})
import numpy as np
import rasterio
from affine import Affine
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from pyproj import CRS, Transformer
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "files"
OUT = ROOT / "outputs"

# --------------------------------------------------------------------------- the four files
# What each file is (found in step 1) and how to tell real data from "no data".
FILES = {
    "file1": {"satellite": "Meteosat 0 deg (MSG)", "instrument": "SEVIRI", "channel": "IR 10.8 um",
              "kind": "ir", "fill": ("value", 0), "sub_lon": 0.0},
    "file2": {"satellite": "Himawari-9", "instrument": "AHI", "channel": "IR 10.4 um",
              "kind": "ir", "fill": ("value", 0), "sub_lon": 140.7},
    "file3": {"satellite": "GOES-East", "instrument": "ABI", "channel": "water vapour 6-7 um",
              "kind": "wv", "fill": ("border", 255), "sub_lon": -75.0,
              # GeoTIFF cannot store the GOES scan "sweep" axis (x), so it is set here.
              "crs": "+proj=geos +lon_0=-75 +h=35786023 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +sweep=x +no_defs"},
    "file4": {"satellite": "polar orbiter", "instrument": "AVHRR-class",
              "channel": "RGB: VIS 0.6 um, NIR 0.8 um, inverted IR 10.8 um", "kind": "rgb", "fill": ("alpha", 4)},
}

# --------------------------------------------------------------------------- common grid (step 4)
RES = 0.03                                       # degrees (~3.3 km), about SEVIRI's 3 km
GRID_TRANSFORM = Affine(RES, 0, -180, 0, -RES, 90)
GRID_SHAPE = (round(180 / RES), round(360 / RES))  # 6000 rows x 12000 columns
GRID_LAT = 90 - (np.arange(GRID_SHAPE[0]) + 0.5) * RES
GRID_LON = -180 + (np.arange(GRID_SHAPE[1]) + 0.5) * RES

# --------------------------------------------------------------------------- calibration (step 2)
# Scenes with a known 10.8 um brightness temperature.
# "fit" anchors determine the line, "check" anchors test it.
ANCHORS = [
    # name, file, box [west, east, south, north], DN percentile, expected BT (K), uncertainty (K), role
    ("Clear ocean, S Atlantic", "file1", (-10, 0, -25, -15), 5, 294.0, 2.5, "fit"),
    ("Clear ocean, W Pacific", "file2", (150, 170, 15, 25), 5, 296.0, 2.5, "fit"),
    ("Coldest storm tops, Maritime Continent", "file2", (100, 150, -10, 5), 99.9, 192.0, 5.0, "fit"),
    ("Antarctic plateau", "file1", (0, 30, -80, -75), 50, 240.0, 8.0, "check"),
    ("Sahara before sunrise", "file1", (0, 20, 20, 28), 50, 281.0, 5.0, "check"),
    ("Central Australia, late afternoon", "file2", (125, 140, -28, -22), 50, 315.0, 8.0, "check"),
]
BT_UNCERTAINTY = 5.0  # K

# --------------------------------------------------------------------------- display + analysis
BT_RANGE = (184.0, 330.0)  # K, the recovered stretch
GLACIATED = 235.0          # K, colder cloud tops are ice (about -38 C)
VZA_LIMIT = 70.0           # degrees, beyond this the view is too oblique for quantitative use


def read(name):
    """One input file: bands, valid-pixel mask, transform, CRS and observation time."""
    info = FILES[name]
    with rasterio.open(DATA / name) as src:
        bands = src.read()
        transform, tag = src.transform, src.tags().get("TIFFTAG_DATETIME")
        crs = CRS.from_user_input(info.get("crs") or src.crs)
    rule, value = info["fill"]
    if rule == "alpha":    # transparent pixels are outside the swath
        valid, bands = bands[value - 1] > 0, bands[:3]
    elif rule == "value":  # this value only ever marks space
        valid = bands[0] != value
    else:                  # "border": the value also marks saturated cloud tops,
        labels, _ = ndimage.label(bands[0] == value)  # so only regions touching the border are space
        edge = np.unique(np.r_[labels[0], labels[-1], labels[:, 0], labels[:, -1]])
        valid = ~np.isin(labels, edge[edge > 0])
    time = datetime.strptime(tag, "%Y:%m:%d %H:%M:%S") if tag else None
    return {"name": name, "bands": bands, "valid": valid, "transform": transform, "crs": crs, "time": time}


def box_values(f, box, band=0, n=200):
    """Valid pixel values of file `f` (from read) inside a lon/lat box, sampled on an n x n lattice."""
    west, east, south, north = box
    lon, lat = np.meshgrid(np.linspace(west, east, n), np.linspace(south, north, n))
    x, y = Transformer.from_crs("EPSG:4326", f["crs"], always_xy=True).transform(lon, lat)
    if f["crs"].coordinate_operation.method_name == "Equidistant Cylindrical":
        half_world = np.pi * f["crs"].ellipsoid.semi_major_metre  # file2 runs past 180 E
        x = np.where(x < f["transform"].c, x + 2 * half_world, x)
    ok = np.isfinite(x) & np.isfinite(y)
    rows, cols = rasterio.transform.rowcol(f["transform"], x[ok], y[ok])
    rows, cols = np.asarray(rows), np.asarray(cols)
    h, w = f["valid"].shape
    inside = (rows >= 0) & (rows < h) & (cols >= 0) & (cols < w)
    rows, cols = rows[inside], cols[inside]
    keep = f["valid"][rows, cols]
    return f["bands"][band][rows[keep], cols[keep]].astype(float)


def ir_colormap(n=512):
    """Grey for warm (dark) to cold (white), then a purple ramp for ice-topped cloud below 235 K."""
    lo, hi = BT_RANGE
    bt = np.linspace(lo, hi, n)
    warm = plt.get_cmap("Greys")(0.05 + 0.85 * (bt - GLACIATED) / (hi - GLACIATED))
    cold = plt.get_cmap("Purples")(0.25 + 0.75 * (GLACIATED - bt) / (GLACIATED - lo))
    return ListedColormap(np.where((bt < GLACIATED)[:, None], cold, warm))


def wv_colormap():
    """Water vapour, as forecasters view it: dry (dark) to moist / high cloud (white). High DN = cold = moist."""
    return LinearSegmentedColormap.from_list("wv", ["#1b1f3b", "#35608f", "#8fbbe0", "#ffffff"])


def pretty(text):
    """ASCII metadata text (kept plain for NetCDF tools) to figure text: 'um' -> 'µm', ' deg' -> '°'."""
    return text.replace(" um", " µm").replace(" deg", "°")


def map_gridlines(ax):
    """Light lat/lon gridlines with degree labels on the left and bottom."""
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color="0.5", alpha=0.6, linestyle=":")
    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = gl.ylabel_style = {"size": 7, "color": "0.3"}
    return gl


def out_dir(step):
    d = OUT / step
    d.mkdir(parents=True, exist_ok=True)
    return d
