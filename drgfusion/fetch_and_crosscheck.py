"""Fetch Brouwer's DRG and SRG tables and cross-check every strongly regular fusion against the SRG table.

1. Downloads the machine-readable DRG tables (https://aeb.win.tue.nl/drg/) into cache/drg/ (file names get a
   suffix _ext for the extended tables IA4a, IA5a, because IA4A and IA4a clash on case-insensitive disks).
2. Runs drg_table_fusions.py on them (closed fusions, diameter <= 7).
3. Downloads the SRG table pages (v <= 1300, https://aeb.win.tue.nl/graphs/srg/) into cache/srg/, reads the
   row colour (green: exists, red: does not exist, yellow: open) and the comment, and looks up each fused SRG.
Writes srg_crosscheck.csv (every fused SRG with v <= 1300) and nonexistent_fusions.csv (rows whose fused SRG
is marked as not existing). cache/ is not committed.

Usage: python fetch_and_crosscheck.py   (run from drgfusion/)
"""
import os, re, csv, html, subprocess, sys
import requests

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache')
DRG = ['IA3P', 'IA4A', 'IA4B', 'IA4C', 'IA4P', 'IA5A', 'IA5B', 'IA5C', 'IA5P', 'IA6A', 'IA6B', 'IA6C', 'IA6P',
       'IA7A', 'IA7C', 'IA7P', 'IA8', 'IA9', 'IA10', 'IA11', 'IA12', 'IA13', 'IA4a', 'IA5a']
COLOUR = {'#b0ffb0': 'exists', '#ffb0b0': 'does not exist', '#ffffb0': 'open'}


def fetch(url, path):
    if not os.path.exists(path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        r.encoding = 'utf-8'
        with open(path, 'w', encoding='utf8') as f:
            f.write(r.text)
    return open(path, encoding='utf8').read()


def srg_table():
    idx = fetch('https://aeb.win.tue.nl/graphs/srg/srgtab.html', os.path.join(CACHE, 'srg', 'srgtab.html'))
    pages = sorted(set(re.findall(r'href="(srgtab\d+-\d+\.html)"', idx)))
    rec = {}
    for pg in pages:
        h = fetch('https://aeb.win.tue.nl/graphs/srg/' + pg, os.path.join(CACHE, 'srg', pg))
        v = None
        for m in re.finditer(r'<tr([^>]*)>(.*?)</tr>', h, re.S):
            attrs, body = m.groups()
            cells = [html.unescape(re.sub('<[^>]+>', '', c)).strip()
                     for c in re.findall(r'<td[^>]*>(.*?)</td>', body, re.S)]
            if len(cells) < 5:
                continue
            col = re.search(r'bgcolor="([^"]+)"', attrs)
            try:
                if cells[1]:
                    v = int(cells[1])
                k, lam, mu = int(cells[2]), int(cells[3]), int(cells[4])
            except ValueError:
                continue
            refs = re.findall(r'href="srgtabrefs.html#([^"]+)"', body)
            rec[(v, k, lam, mu)] = (COLOUR.get(col.group(1) if col else '', 'unknown'), cells[-1], ' '.join(refs))
    return rec


def main():
    tdir = os.path.join(CACHE, 'drg')
    for t in DRG:
        fetch('https://aeb.win.tue.nl/drg/' + t, os.path.join(tdir, t + ('_ext' if t[-1] == 'a' else '') + '.txt'))
    out = os.path.join(HERE, 'table_fusions.csv')
    subprocess.run([sys.executable, os.path.join(HERE, 'drg_table_fusions.py'), tdir, out], check=True)
    rec = srg_table()
    rows, bad = [], []
    for r in csv.DictReader(open(out)):
        if not r['srg_parameters']:
            continue
        for sp in r['srg_parameters'].split(' / '):
            key = tuple(int(x) for x in re.findall(r'\d+', sp))
            if key[0] > 1300:
                continue
            st = rec.get(key, ('not listed (trivial or infeasible)', '', ''))
            row = [r['table'], r['array'], r['fusion'], sp, st[0], st[1], st[2]]
            rows.append(row)
            if st[0] == 'does not exist':
                bad.append(row)
    hdr = ['table', 'array', 'fusion', 'fused_srg', 'srg_table_status', 'srg_table_comment', 'srg_table_refs']
    for fn, data in (('srg_crosscheck.csv', rows), ('nonexistent_fusions.csv', bad)):
        with open(os.path.join(HERE, fn), 'w', newline='', encoding='utf8') as f:
            w = csv.writer(f)
            w.writerow(hdr)
            w.writerows(data)
    print(f'{len(rows)} fused SRGs with v <= 1300; {len(bad)} of them marked as not existing; '
          f'{len({b[1] for b in bad})} distinct arrays')


if __name__ == '__main__':
    main()
