#!/usr/bin/env python3
"""Extract the two explicit witnesses from the supplied November 5 PDF.

This one-off recovery utility accepts a pdftotext -layout transcription.
The recovered integers are subsequently verified without the PDF.
"""
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def recover(text):
    section = text.split('Theorem J.1.')[1].split('Proof.')[0]
    rows = []
    for line in section.splitlines():
        nums = re.findall(r'(?<!\d)-?\d+(?!\d)', line.split('=')[-1].replace('−', '-'))
        if len(nums) == 25 and set(nums) <= {'0', '1', '-1'}:
            rows.append(list(map(int, nums)))
    assert len(rows) == 25, len(rows)
    digits = re.findall(r'^\s*([0-5]{15})\s*$', text, re.M)
    assert len(digits) == 15, len(digits)
    return rows, [list(map(int, row)) for row in digits]


if __name__ == '__main__':
    w, h = recover(Path(sys.argv[1]).read_text())
    target = ROOT / 'data/revision'
    target.mkdir(parents=True, exist_ok=True)
    for name, value in [('bw25_16_witness', w), ('sixth_root_order15', h)]:
        (target / (name + '.json')).write_text(json.dumps(value, indent=2) + '\n')
    print('Recovered 625 signed entries and 225 sixth-root exponents.')
