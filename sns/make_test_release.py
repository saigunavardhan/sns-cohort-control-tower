import argparse
import json
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from chaos import DATA, Base, Chaos  # noqa: E402
from engine.impact import write_layer  # noqa: E402
from engine.release import load_release  # noqa: E402


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--max-ops", type=int, default=4)
    args = p.parse_args()

    rel_a = load_release(DATA / "raw")
    chaos = Chaos(Base(rel_a), args.seed)
    chaos.scenario(args.max_ops)

    out = Path(args.out)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    for name, layer in rel_a.items():
        if name in chaos.dropped_files:
            continue
        if name in chaos.tables:
            write_layer(layer, chaos.tables[name], out / name)
        else:
            shutil.copy2(layer.path, out / name)

    key_path = out.parent / f"{out.name}.ANSWER_KEY.json"
    key_path.write_text(json.dumps({"seed": args.seed, "changes": chaos.changes}, indent=2), encoding="utf-8")
    print(json.dumps({"release": str(out), "answer_key": str(key_path), "n_changes": len(chaos.changes)}, indent=2))


if __name__ == "__main__":
    main()