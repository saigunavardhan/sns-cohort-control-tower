import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from chaos import score  # noqa: E402


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--key", required=True)
    p.add_argument("--findings", required=True)
    args = p.parse_args()

    key = json.loads(Path(args.key).read_text(encoding="utf-8"))
    doc = json.loads(Path(args.findings).read_text(encoding="utf-8"))
    rows, false_findings, extra_ids = score(key["changes"], doc)

    print(json.dumps({
        "injected": [f"{c['type']}:{c['file']}" for c in key["changes"]],
        "caught": f"{sum(r['caught'] for r in rows)}/{len(rows)}",
        "missed": [f"{r['type']}:{r['file']}" for r in rows if not r["caught"]],
        "false_findings": false_findings,
        "extra_sample_ids": extra_ids,
        "verdict": doc["verdict"],
        "findings": [f"{f['id']} {f['kind']} {f['severity']}: {f['summary']}" for f in doc["findings"]],
    }, indent=2))


if __name__ == "__main__":
    main()