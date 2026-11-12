#DIR="/home/venkatesh/Documents/IISc/data/savedata_parametric_fast_2026-10-28_20-49-39/100Hz_withoutLineNoise_keeptrials_no"
DIR = "/home/venkatesh/Documents/IISc/data/savedata_parametric_fast_2026-10-28_20-53-55"

python3 - "$DIR" <<'PY'
from pathlib import Path
from datetime import datetime
from collections import Counter
import sys

root = Path(sys.argv[1])
files = sorted(
    (p for p in root.rglob("*") if p.is_file()),
    key=lambda p: p.stat().st_mtime
)

if not files:
    print("No files found.")
    sys.exit()

def timestamp(p):
    return datetime.fromtimestamp(p.stat().st_mtime)

start = timestamp(files[0])
end = timestamp(files[-1])
duration = (end - start).total_seconds()

print("Total files:", len(files))
print("First saved file:", start)
print("Last saved file: ", end)
print("Elapsed between first and last saved files:",
      round(duration / 60, 2), "minutes")

print("\nFile counts by extension:")
print(Counter(p.suffix for p in files))

print("\nFile counts by parent folder:")
for folder, count in Counter(str(p.parent.relative_to(root))
                             for p in files).items():
    print(count, folder)

print("\nFirst 15 files:")
for p in files[:15]:
    print(timestamp(p), p.relative_to(root))

print("\nLast 15 files:")
for p in files[-15:]:
    print(timestamp(p), p.relative_to(root))
PY