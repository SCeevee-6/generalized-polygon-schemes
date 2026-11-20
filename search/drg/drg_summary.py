"""Combine the census outputs into one line per array (drg_summary.csv).

Inputs: survivors.csv, drg_results.jsonl (drg_search.py), drg_affine_results.jsonl (drg_affine.py),
array_status.csv (status of unrealised arrays from the literature, entered by hand).
Columns: v, array, groups (primitive groups that passed the filter), realising (those that carry a graph with this
array), status, source. Hamming arrays with q != 4 are left out (see drg_filter.py).

Usage: python drg_summary.py
"""
import os, csv, json
from collections import defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    groups = defaultdict(set); real = defaultdict(set); checked = defaultdict(set)
    for r in csv.DictReader(open(os.path.join(HERE, 'survivors.csv'))):
        if not r['hamming']:
            groups[(int(r['v']), r['array'])].add(int(r['i']))
    for fn, key in (('drg_results.jsonl', 'n_hits'), ('drg_affine_results.jsonl', 'status')):
        for l in open(os.path.join(HERE, fn)):
            r = json.loads(l)
            for e in r['results']:
                k = (r['v'], e['array'])
                if e.get('status') in ('LIMIT',):
                    continue
                checked[k].add(r['i'])
                if (key == 'n_hits' and e['n_hits']) or (key == 'status' and e['status'] == 'REALISED'):
                    real[k].add(r['i'])
    st = {r['array']: (r['status'], r['source']) for r in csv.DictReader(open(os.path.join(HERE, 'array_status.csv')))}
    rows = []
    for k in sorted(groups):
        v, a = k
        unchecked = groups[k] - checked[k]
        if real[k]:
            s, src = 'exists (vertex-primitive)', ''
        elif unchecked:
            s, src = 'not finished', ''
        else:
            s, src = st.get(a, ('no vertex-primitive graph', 'not in the sage-drg nonexistence lists or the DistanceRegular.org catalogue'))
        rows.append({'v': v, 'array': a, 'groups': len(groups[k]), 'unchecked': len(unchecked),
                     'realising': ' '.join(map(str, sorted(real[k]))), 'status': s, 'source': src})
    with open(os.path.join(HERE, 'drg_summary.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    from collections import Counter
    print(len(rows), 'arrays;', dict(Counter(r['status'] for r in rows)))


if __name__ == '__main__':
    main()
