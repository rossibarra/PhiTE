"""Shared fixture helper: stamp a test ancestral table as a valid v2 table."""

import json
from pathlib import Path

import numpy as np

from normalize_tes.build_ancestral_states import ARRAY_NAMES, SCHEMA_VERSION
from normalize_tes.vcf_eligibility import _sha256_array


def stamp_ancestral_table(table: Path) -> Path:
    """Record the v2 schema and the current digest of each array.

    Fixtures edit their arrays after writing them, so this is called once the
    arrays are final (and again after any later edit).
    """
    table = Path(table)
    metadata_path = table / "metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["schema_version"] = SCHEMA_VERSION
    metadata["array_sha256"] = {
        name: _sha256_array(np.load(table / f"{name}.npy", allow_pickle=False))
        for name in ARRAY_NAMES
    }
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
    return table
