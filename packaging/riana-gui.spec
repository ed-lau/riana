# -*- mode: python ; coding: utf-8 -*-

"""PyInstaller spec for the standalone Riana GUI bundle.

Build (from the repo root, in an env with ``.[gui]`` + PyInstaller installed):

    pyinstaller packaging/riana-gui.spec --noconfirm

Output:
    dist/riana-gui/              one-dir bundle (all platforms)
    dist/Riana.app              macOS app bundle (double-clickable; wraps the above)

The bundle is CLI+GUI: no args -> GUI, args -> the ``riana`` CLI (see riana_gui.py).
"""

import os
import sys

from PyInstaller.utils.hooks import (
    collect_all, collect_data_files, collect_submodules,
)

# SPECPATH is the directory containing this spec file (injected by PyInstaller).
here = SPECPATH  # noqa: F821  (PyInstaller global)
launcher = os.path.join(here, "riana_gui.py")

# Package data files, keeping their in-package paths: the coefficient CSVs
# (riana/data/coefficients/*.csv) and the GUI icon (riana/gui/resources/*.png),
# both loaded at runtime via importlib.resources. NB: pass no `includes=` filter —
# its fnmatch is anchored oddly and drops files in subdirectories; the riana package
# ships only these CSV/PNG data files, so an unfiltered collect grabs exactly them.
datas = collect_data_files("riana")

# The ProcessPoolExecutor workers re-import riana.* by name, so every submodule
# must be bundled even though it is not statically imported from the launcher.
hiddenimports = (
    collect_submodules("riana")
    + collect_submodules("qasync")
    + collect_submodules("pyqtgraph")
)

# IsoSpecPy carries a compiled extension + data; pull the whole package in.
binaries = []
for _pkg in ("IsoSpecPy",):
    _b, _d, _h = collect_all(_pkg)
    binaries += _b
    datas += _d
    hiddenimports += _h

a = Analysis(
    [launcher],
    pathex=[],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    # Trim weight: no Tk GUI, no test tooling in the shipped bundle.
    excludes=["tkinter", "pytest", "_pytest"],
    noarchive=False,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="riana-gui",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,          # windowed app (no terminal on double-click)
    disable_windowed_traceback=False,
    argv_emulation=True,    # macOS: deliver Finder file-open args as argv
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    name="riana-gui",
)

# macOS: wrap the one-dir bundle in a double-clickable .app. icon=None uses the
# default; drop a packaging/riana.icns here and set icon= to brand the Dock tile
# (the window/Dock icon is also set at runtime from riana/gui/resources/riana.png).
if sys.platform == "darwin":
    app = BUNDLE(
        coll,
        name="Riana.app",
        icon=None,
        bundle_identifier="net.laulab.riana",
        info_plist={
            "NSHighResolutionCapable": True,
            "CFBundleShortVersionString": os.environ.get("RIANA_VERSION", "0.0.0"),
            "LSMinimumSystemVersion": "11.0",
        },
    )
