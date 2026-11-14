import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))

from engine.diff import compare_releases  # noqa: E402
from engine.release import load_release  # noqa: E402


def main():
    rel_a = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "data" / "tcga-brca" / "raw"
    rel_b = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "data" / "tcga-brca" / "releases" / "practice-B"
    report = compare_releases(load_release(rel_a), load_release(rel_b))
    report = {"release_a": str(rel_a), "release_b": str(rel_b), **report}
    print(json.dumps(report, indent=2, default=str))


if __name__ == "__main__":
    main()
