"""Locating prepared receptors under 03_receptors/, wherever they sit.

The 12 receptor directories are grouped by target family (ttbk/, mao/). Nothing here knows
which family a structure belongs to, and that is the point: resolution is by glob at either
depth, so regrouping a structure cannot make it silently unfindable, and a half-migrated tree
still resolves. The family rule lives in exactly one place -- receptor_group in
scripts/receptor_paths.sh -- and is used only when creating a receptor directory.

Importable as `import project_paths` from anything run as `python scripts/<x>.py`, since
Python puts the script's own directory on sys.path.
"""
import glob
import os

RECEPTOR_ROOT = "03_receptors"


def receptor_dir(pdb, root=RECEPTOR_ROOT):
    """Directory holding this receptor, or None.

    Raises on an ambiguous hit: two copies of one receptor would make it undefined which
    preparation a score was computed against, and silently picking one is how this project
    has produced false positives before.
    """
    hits = [p for p in (glob.glob(os.path.join(root, pdb)) +
                        glob.glob(os.path.join(root, "*", pdb))) if os.path.isdir(p)]
    hits = sorted(set(hits))
    if not hits:
        return None
    if len(hits) > 1:
        raise RuntimeError(f"receptor {pdb!r} resolves to {len(hits)} directories: {hits}")
    return hits[0]


def receptor_file(pdb, name="receptor.pdbqt", root=RECEPTOR_ROOT):
    """Path to a file inside a receptor's directory, or None if the receptor is absent.

    Returns the path even when the file itself does not exist, so callers keep their own
    'is it there?' check and their own error message.
    """
    d = receptor_dir(pdb, root)
    return None if d is None else os.path.join(d, name)
