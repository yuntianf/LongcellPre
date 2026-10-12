"""Reject new/changed check warnings; retain one exact, documented legacy warning."""
from pathlib import Path
import json
import re
import sys

def check(log, allowed):
    if re.search(r'^\* [^\n]+ \.\.\. ERROR\s*$',log,re.M):
        raise ValueError('R CMD check reports an error')
    warnings=dict(re.findall(r'^\* ([^\n]+) \.\.\. WARNING\n(.*?)(?=^\* |\Z)',log,re.M|re.S))
    unexpected=[name for name,body in warnings.items()
                if re.sub(r'\s+',' ',body).strip()!=allowed.get(name)]
    if unexpected:
        raise ValueError('New or changed R CMD check warnings: '+', '.join(unexpected))
    return list(warnings)

if __name__=='__main__':
    allowed=json.loads((Path(__file__).with_name('known-warnings.json')).read_text())
    try:
        warnings=check(Path(sys.argv[1]).read_text(),allowed)
    except ValueError as error:
        raise SystemExit(str(error))
    print('Check diagnostic gate passed. Known warning sections:',warnings)
