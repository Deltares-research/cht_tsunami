"""GEM Global Active Faults database download and loading utilities.

Downloads the GEM fault database from S3 on first use and caches it
locally. Provides a GeoDataFrame with fault traces and attributes
(dip, rake, slip rate, etc.) compatible with Okada parameters.
"""

from pathlib import Path
from typing import Any, Dict, Optional, Union

import geopandas as gpd
import numpy as np

# Default S3 location for the GEM faults file
S3_BUCKET = "delftdashboard"
# Endpoint URL for S3-compatible stores; set to None for AWS S3
S3_ENDPOINT = "https://s3.deltares.nl"
S3_KEY = "data/tsunami/gem_active_faults.geojson"
S3_REGION = "eu-west-1"
FILE_NAME = "gem_active_faults.geojson"


def get_faults(
    path: Union[str, Path],
    s3_bucket: str = S3_BUCKET,
    s3_key: str = S3_KEY,
    s3_region: str = S3_REGION,
    check_online: bool = True,
) -> Optional[gpd.GeoDataFrame]:
    """Load the GEM Global Active Faults database, downloading from S3 if needed.

    Parameters
    ----------
    path : str or Path
        Local directory to store the GeoJSON file.
    s3_bucket : str
        S3 bucket name.
    s3_key : str
        S3 object key.
    s3_region : str
        AWS region.
    check_online : bool
        If True, download from S3 when the local file is missing.

    Returns
    -------
    gpd.GeoDataFrame or None
        Fault traces with attributes, or None if unavailable.
    """
    path = Path(path)
    local_file = path / FILE_NAME

    if not local_file.exists() and check_online:
        _download_from_s3(local_file, s3_bucket, s3_key, s3_region)

    if local_file.exists():
        try:
            gdf = gpd.read_file(local_file)
            gdf["index"] = range(len(gdf))
            return gdf
        except Exception as e:
            print(f"Could not read GEM faults file: {e}")
            return None

    print(f"GEM faults file not found: {local_file}")
    return None


def _download_from_s3(
    local_file: Path,
    s3_bucket: str,
    s3_key: str,
    s3_region: str,
) -> None:
    """Download the GEM faults GeoJSON from S3 (unsigned/public access).

    Parameters
    ----------
    local_file : Path
        Target file path.
    s3_bucket : str
        S3 bucket name.
    s3_key : str
        S3 object key.
    s3_region : str
        AWS region.
    """
    try:
        import boto3
        from botocore import UNSIGNED
        from botocore.client import Config

        print(f"Downloading GEM faults from s3://{s3_bucket}/{s3_key} ...")
        local_file.parent.mkdir(parents=True, exist_ok=True)

        s3 = boto3.client(
            "s3",
            endpoint_url=S3_ENDPOINT or None,
            region_name=s3_region,
            config=Config(signature_version=UNSIGNED),
        )
        s3.download_file(s3_bucket, s3_key, str(local_file))
        print(f"Downloaded to {local_file}")

    except Exception as e:
        print(f"Could not download GEM faults from S3: {e}")


def get_okada_params_from_fault(gdf: gpd.GeoDataFrame, index: int) -> Dict[str, Any]:
    """Extract Okada-compatible parameters from a GEM fault feature.

    Computes source location (midpoint), strike (from geometry azimuth),
    fault length (from geometry in km), depth (from seismogenic range),
    dip, and rake from the GEM database attributes.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        GEM faults GeoDataFrame (as returned by :func:`get_faults`).
    index : int
        Row index of the fault to extract.

    Returns
    -------
    dict
        Okada parameters with keys: ``longitude``, ``latitude``,
        ``strike``, ``length``, ``depth``, ``dip``, ``rake``,
        ``name``, ``slip_type``.  Values are ``None`` when not
        available from the database.
    """
    if gdf is None or index >= len(gdf):
        return {}

    row = gdf.iloc[index]
    geom = row.geometry

    # Midpoint of the fault trace
    midpoint = geom.interpolate(0.5, normalized=True)

    # Strike from geometry azimuth (start → end)
    coords = list(geom.coords)
    strike = None
    if len(coords) >= 2:
        dx = coords[-1][0] - coords[0][0]
        dy = coords[-1][1] - coords[0][1]
        strike = round((90.0 - np.degrees(np.arctan2(dy, dx))) % 360, 1)

    # Length from geometry (degrees → km, rough mid-latitude estimate)
    length = None
    if geom.length > 0:
        mid_lat = midpoint.y
        km_per_deg = 111.0 * np.cos(np.radians(mid_lat))
        # Use average of lon and lat scaling
        length = round(geom.length * (111.0 + km_per_deg) / 2, 1)

    # Depth from seismogenic depth range
    upper = parse_gem_tuple(row.get("upper_seis_depth"))
    lower = parse_gem_tuple(row.get("lower_seis_depth"))
    if upper is not None and lower is not None:
        depth = round((upper + lower) / 2, 1)
    elif upper is not None:
        depth = upper
    else:
        depth = None

    return {
        "longitude": round(midpoint.x, 4),
        "latitude": round(midpoint.y, 4),
        "strike": strike,
        "length": length,
        "depth": depth,
        "dip": parse_gem_tuple(row.get("average_dip")),
        "rake": parse_gem_tuple(row.get("average_rake")),
        "net_slip_rate": parse_gem_tuple(row.get("net_slip_rate")),
        "name": row.get("name", "Unknown"),
        "slip_type": row.get("slip_type", ""),
    }


def parse_gem_tuple(value) -> Optional[float]:
    """Parse a GEM attribute tuple string like ``'(38,,)'`` or ``'(1.55,0.8,2.22)'``.

    Returns the first (preferred) value as a float, or None if not parseable.

    Parameters
    ----------
    value : any
        Raw attribute value from the GEM GeoDataFrame.

    Returns
    -------
    float or None
        Parsed preferred value.
    """
    if value is None:
        return None
    s = str(value).strip("() ")
    parts = s.split(",")
    for part in parts:
        part = part.strip()
        if part:
            try:
                return float(part)
            except ValueError:
                continue
    return None
