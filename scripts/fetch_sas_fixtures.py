"""Download only public pyreadstat test fixtures and provenance; never patient files."""
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

root = Path(__file__).resolve().parents[1]
folder = root / "tests/fixtures/sas"
folder.mkdir(parents=True, exist_ok=True)
# Paths verified against the upstream tests/test_basic.py. No API token is required.
selected = ["test_data/basic/sample.sas7bdat", "test_data/basic/dates.sas7bdat",
            "test_data/basic/sample_bincompressed.sas7bdat", "test_data/missing_data/missing_test.sas7bdat"]
manifest = []
for p in selected:
    url = "https://raw.githubusercontent.com/Roche/pyreadstat/master/" + p
    body = urlopen(url, timeout=30).read()
    if len(body) > 2_000_000:
        continue
    target = folder / (p.replace("/", "__"))
    target.write_bytes(body)
    manifest.append({"file": target.name, "url": url, "upstream_ref": "master; content hash frozen below",
                     "bytes": len(body), "sha256": hashlib.sha256(body).hexdigest()})
for filename in ["LICENSE"]:
    body = urlopen("https://raw.githubusercontent.com/Roche/pyreadstat/master/" + filename, timeout=30).read()
    (folder / "UPSTREAM_LICENSE").write_bytes(body)
(folder / "manifest.json").write_text(json.dumps(manifest, indent=2))
print(f"Saved {len(manifest)} public fixtures")
