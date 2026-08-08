# Packaging — standalone Riana GUI bundle

Build a self-contained Riana GUI that non-Python users can launch without a
Python install (a macOS `.app`, or a one-dir bundle on Linux/Windows), using
[PyInstaller](https://pyinstaller.org).

## Build

From the repo root, in an environment with the GUI + packaging extras:

```bash
pip install -e ".[gui,packaging]"
packaging/build_app.sh
```

or invoke PyInstaller directly:

```bash
pyinstaller packaging/riana-gui.spec --noconfirm --clean
```

Output:

- `dist/riana-gui/` — the one-dir bundle (all platforms). Run
  `dist/riana-gui/riana-gui`.
- `dist/Riana.app` — a double-clickable macOS app (wraps the one-dir bundle).

## What the bundle is

The frozen binary is **CLI + GUI** (`packaging/riana_gui.py`):

- launched with **no arguments** (a Finder/Explorer double-click) → opens the GUI;
- launched **with arguments** → the full `riana` CLI, e.g.
  `dist/riana-gui/riana-gui integrate ...` or `riana-gui --version`.

`multiprocessing.freeze_support()` is called first so the GUI's
`ProcessPoolExecutor` workers (which re-execute this binary) run their payload
instead of relaunching the app.

## What the spec bundles

- **Data:** the coefficient tables (`riana/data/coefficients/*.csv`) and the GUI
  icon (`riana/gui/resources/*.png`), collected with their in-package paths so the
  runtime `importlib.resources` loads still resolve.
- **Hidden imports:** all `riana.*` submodules (the pool workers import them by
  name), plus `qasync` / `pyqtgraph`; `IsoSpecPy` is pulled in whole
  (`collect_all`) for its compiled extension. PySide6 is handled by PyInstaller's
  bundled hook.

## Verification status

The build and a **CLI smoke test** (`riana-gui --version`, no display needed) are
the automated checks. **Launching the GUI window must be verified manually** on a
machine with a display — CI here is headless. When you do:

1. `packaging/build_app.sh`
2. `dist/riana-gui/riana-gui --version` → prints the version (freeze + data files OK).
3. Double-click `dist/Riana.app` (macOS) or run `dist/riana-gui/riana-gui` → the
   window opens; exercise Integrate/Model/Protein once (the ProcessPoolExecutor
   path is the thing a frozen build most often breaks).

## Optional: macOS Dock icon

`icon=None` in the spec uses the default tile. To brand it, add a
`packaging/riana.icns` (convert `riana/gui/resources/riana.png` with
`sips`/`iconutil`) and set `icon="packaging/riana.icns"` in the `BUNDLE(...)` call.
The in-app window/Dock icon is already set at runtime from the packaged PNG.
