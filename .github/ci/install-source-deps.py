"""Install the four source-only dependencies using fixed versions and SHA-256."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import urllib.request

ci = Path(__file__).resolve().parent
cache = Path(os.environ.get('CI_DOWNLOAD_CACHE', 'ci-output/dependency-cache')).resolve()
cache.mkdir(parents=True, exist_ok=True)
for package in json.loads((ci/'source-dependencies.json').read_text()):
    target = cache/package['filename']
    def valid():
        return target.exists() and hashlib.sha256(target.read_bytes()).hexdigest() == package['sha256']
    if not valid():
        failures = []
        for url in package['urls']:
            try:
                with urllib.request.urlopen(url, timeout=120) as response:
                    data = response.read()
                if hashlib.sha256(data).hexdigest() != package['sha256']:
                    raise ValueError('SHA-256 mismatch')
                target.write_bytes(data)
                break
            except Exception as error:
                failures.append(f'{url}: {error}')
        if not valid():
            raise SystemExit('\n'.join(failures))
    subprocess.run(['R','CMD','INSTALL','--preclean',str(target)],check=True)
print('All source dependency archives verified and installed.')
