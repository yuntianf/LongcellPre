"""The same CI entry point runs on GitHub Actions and locally under Linux."""
from pathlib import Path
import hashlib
import json
import os
import re
import subprocess
import sys

source = Path(__file__).resolve().parents[2]
ci = source/'.github/ci'
output = Path(os.environ.get('LONGCELL_CI_OUTPUT', source/'ci-output')).resolve()
output.mkdir(parents=True, exist_ok=True)
library = output/'library'
library.mkdir(exist_ok=True)
env = os.environ.copy()
env.update(R_LIBS_USER=str(library), _R_CHECK_FORCE_SUGGESTS_='false',
           OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', TZ='UTC', LC_ALL='C.UTF-8')
# Disable package-index checks during R CMD check; dependencies are pinned above.
profile = output/'profile.R'
profile.write_text('options(repos=NULL)\n')
env['R_PROFILE_USER'] = str(profile)
status = {}
def run(command, label):
    print(label, flush=True)
    with (output/(label+'.log')).open('w') as log:
        result = subprocess.run(command, cwd=output, env=env, stdout=log, stderr=subprocess.STDOUT)
    status[label] = result.returncode
    (output/'status.json').write_text(json.dumps(status, indent=2)+'\n')
    if result.returncode:
        print((output/(label+'.log')).read_text()[-12000:])
        raise SystemExit(f'{label} failed with exit code {result.returncode}')

run(['Rscript','-e','stopifnot(!length(find.package("Longcellsrc",quiet=TRUE))); '
                  'stopifnot(all(nzchar(Sys.which(c("minimap2","samtools","bedtools","bash","python")))))'], 'preflight')
run(['R','CMD','build','--no-build-vignettes',str(source)], 'build')
version = re.search(r'^Version:\s*(\S+)',(source/'DESCRIPTION').read_text(),re.M).group(1)
tarball = output/f'LongcellPre_{version}.tar.gz'
run(['R','CMD','INSTALL','--preclean','--install-tests','--library='+str(library),str(tarball)], 'install')
run(['Rscript',str(ci/'test-installed.R'),str(output/'test-results.rds'),str(source/'DESCRIPTION'),
     str(output/'session.txt'),str(output/'r-packages.csv')], 'regression-tests')
run(['Rscript',str(ci/'native-equivalence.R'),'updated',str(ci/'fixtures')], 'native-equivalence')
run(['Rscript',str(ci/'iso-equivalence.R'),'updated',str(ci/'fixtures/iso-baseline.rds')], 'iso-equivalence')
run(['R','CMD','check','--no-manual','--ignore-vignettes',str(tarball)], 'package-check')
check_log = (output/'LongcellPre.Rcheck/00check.log').read_text()
run([sys.executable,str(ci/'check-diagnostics.py'),str(output/'LongcellPre.Rcheck/00check.log')], 'check-diagnostics')
warnings = re.findall(r'^\* (.*?) \.\.\. WARNING\s*$',check_log,re.M)
run(['Rscript',str(library/'LongcellPre/exec/RunLongcellPre.R'),'--help','--full'], 'cli-help')
digest = hashlib.sha256(tarball.read_bytes()).hexdigest()
(output/'SHA256SUMS').write_text(f'{digest}  {tarball.name}\n')
(output/'validation.json').write_text(json.dumps({'version':version,'all_checks_passed':True,
    'known_warning_sections':warnings,'source_sha256':digest,
    'platform':'Linux R 4.4 / Bioconductor 3.20'},indent=2)+'\n')
print('PASS: package checks, unskipped regression tests, exact native fixtures and CLI help.',flush=True)
