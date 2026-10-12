"""Write a genuine gzip FASTQ copy, detecting input compression by content."""
import argparse
import gzip
import os
from pathlib import Path
import shutil
import tempfile

def prepare(source, destination):
    source, destination = Path(source), Path(destination)
    if destination.exists():
        raise ValueError(f'Refusing to overwrite {destination}')
    with source.open('rb') as stream:
        compressed = stream.read(2) == b'\x1f\x8b'
    opener = gzip.open if compressed else open
    with opener(source, 'rb') as stream:
        if stream.read(1) != b'@':
            raise ValueError('Input does not start with a FASTQ header')
    descriptor, staged = tempfile.mkstemp(prefix='.longcell-',dir=destination.parent)
    os.close(descriptor)
    try:
        with opener(source,'rb') as inp, open(staged,'wb') as raw:
            with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as out:
                shutil.copyfileobj(inp,out,length=1024*1024)
        # Exclusive publication avoids overwriting an output created meanwhile.
        os.link(staged,destination)
    finally:
        os.unlink(staged)

if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source')
    parser.add_argument('destination')
    args=parser.parse_args()
    prepare(args.source,args.destination)
    print(f'Wrote gzip FASTQ: {args.destination}')
