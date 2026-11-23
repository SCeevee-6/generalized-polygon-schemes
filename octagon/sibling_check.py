"""Run the bare-graph sibling obstruction on every cover in a stage-three file.

python octagon/sibling_check.py octagon/stage3_fp1_switched.json

Each cover is tested with the witness recheck switched on, whatever stage three
said about it. Writes sibling_<input name> next to the input.
"""
import json
import sys
from pathlib import Path

import sibling_obstruction as SO

src = Path(sys.argv[1])
rows = json.loads(src.read_text())
out = []
for i, rec in enumerate(rows):
    count, line, missing = SO.test(rec["vol"], True)
    out.append(dict(row=i, stage3_rejected=rec["rejected"], blocks=count,
                    anchor=line, missing=[int(m) for m in missing]))
src.with_name("sibling_" + src.name).write_text(json.dumps(out))
print("covers", len(out), "left after stage 3", sum(not r["rejected"] for r in rows),
      "left after sibling obstruction", [d["row"] for d in out if d["anchor"] is None])
