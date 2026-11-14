"""Batch refits for the wls-var / D / E evaluation (WP0 / WP2).

Builds relabeled manifests with ``wlsvar_manifests.py`` and runs ``riana fit`` on
each into ``<out>/<test>_<seed>/``. T1 assignments are de-duplicated: an
assignment and its A/B mirror are the same test, so seeds that repeat one are
skipped (boomi has 2⁹ / 2 = 256 distinct assignments). Existing outputs are kept,
so the batch can be resumed. T2 permutations are de-duplicated the same way.

``--stage2 STORE`` replaces ``riana fit`` with ``wlsvar_stage2.py``: only the curve fit and
bootstrap are redone, on the per-point fs of a depth-1 fit of the source manifest
(STORE), under the new labels; the options after ``--`` then go to ``wlsvar_stage2.py``.

python tests/benchmark/wlsvar_refits.py t1 SRC_MANIFEST OUT --n 30 -- <riana fit options>
python tests/benchmark/wlsvar_refits.py t3 SRC_MANIFEST OUT -- <riana fit options>
python tests/benchmark/wlsvar_refits.py t2 SRC_MANIFEST OUT --n 20 --stage2 STORE -- -d 3 -W 14
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

from wlsvar_manifests import _write_manifest, read_manifest, t1, t2, t3  # noqa: E402


def _signature(rows) -> tuple:
    a = sorted((r.identity.labeling_time, r.identity.data_file) for r in rows
               if r.identity.condition == "A")
    b = sorted((r.identity.labeling_time, r.identity.data_file) for r in rows
               if r.identity.condition == "B")
    return tuple(min(a, b))


def _fit(manifest: Path, out: Path, fit_args: list[str], stage2: tuple | None = None) -> float:
    if (out / "riana_fit_peptides.txt").exists():
        print(f"  {out} exists, skipped")
        return 0.0
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if stage2:
        store, src = stage2
        cmd = [sys.executable, str(HERE / "wlsvar_stage2.py"), str(store), str(src),
               str(manifest), str(out), *fit_args]
    else:
        cmd = ["riana", "fit", "--manifest", str(manifest), "-o", str(out), *fit_args]
    with open(out.parent / f"{out.name}.log", "w") as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
    dt = time.time() - t0
    print(f"  {out}: {dt / 60:.1f} min", flush=True)
    return dt


def main() -> None:
    argv = sys.argv[1:]
    fit_args = argv[argv.index("--") + 1:] if "--" in argv else []
    argv = argv[:argv.index("--")] if "--" in argv else argv
    p = argparse.ArgumentParser()
    p.add_argument("test", choices=["t1", "t2", "t3"])
    p.add_argument("src", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--n", type=int, default=20, help="assignments / permutations")
    p.add_argument("--stage2", type=Path, metavar="STORE",
                   help="refit stage 2 only, from this depth-1 fit of SRC_MANIFEST")
    a = p.parse_args(argv)
    stage2 = (a.stage2, a.src) if a.stage2 else None
    rows = read_manifest(a.src, stage="integrate")
    a.out.mkdir(parents=True, exist_ok=True)
    if a.test == "t3":
        for name, rs in zip(["odd", "even"], t3(rows)):
            m = a.out / f"t3_{name}.tsv"
            _write_manifest(m, rs)
            _fit(m, a.out / f"t3_{name}", fit_args, stage2)
        return
    seen, seed, done = set(), 0, 0
    while done < a.n:
        seed += 1
        rs = t1(rows, seed) if a.test == "t1" else t2(rows, seed)
        if a.test == "t1":
            sig = _signature(rs)
        else:  # the individuals in each condition; a relabeling of the groups is the same test
            groups = {}
            for r in rs:
                groups.setdefault(r.identity.condition, set()).add(
                    (r.identity.labeling_time, r.identity.individual))
            sig = min(tuple(sorted(g)) for g in groups.values())
        if sig in seen:
            continue
        seen.add(sig)
        m = a.out / f"{a.test}_{seed:03d}.tsv"
        _write_manifest(m, rs)
        _fit(m, a.out / f"{a.test}_{seed:03d}", fit_args, stage2)
        done += 1


if __name__ == "__main__":
    main()
