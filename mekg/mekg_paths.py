"""Locations of data artifacts; builder code stays inside the repository."""
from pathlib import Path
import os

def workspace():
    return Path(os.environ.get('MEKG_WORKSPACE',Path(__file__).resolve().parents[1]/'data/mekg')).resolve()
