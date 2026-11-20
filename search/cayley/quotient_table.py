"""Quotient obstructions for the arrays that pass the Krein test, each checked with CP-SAT (quotient_lift.py) and
SCIP (quotient_mip.py). Writes quotient_table.csv.

Every abelian group of order 625 has a quotient of order 125, and every abelian group of order 750 or 1000
has its Sylow 5-subgroup (order 125) as a quotient. So if the quotient model is infeasible for Z5^3, Z25 x Z5
and Z125 with fibre f = v/125, no abelian group of that order carries the graph.
For the 729 array, Z9 x Z9 x Z3 (fibre 3) is a quotient of Z9^3 and of Z9^2 x Z3^2.
"""
import csv, json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
CHECKS = [("52,42,16;1,6,28", q, 5) for q in ("5,5,5", "25,5", "125")] + \
         [("49,40,22;1,5,28", q, 6) for q in ("5,5,5", "25,5", "125")] + \
         [("54,42,11;1,3,44", q, 8) for q in ("5,5,5", "25,5", "125")] + \
         [("63,48,24;1,6,28", q, 8) for q in ("5,5,5", "25,5", "125")] + \
         [("56,42,20;1,6,28", "9,9,3", 3)]


def run(script, args, out):
    p = os.path.join(HERE, 'out_table', out)
    if not os.path.exists(p):
        with open(p + '.err', 'w') as fe:
            subprocess.run([sys.executable, script] + args + ['--out', p], cwd=HERE, stdout=fe, stderr=fe)
    return json.load(open(p)) if os.path.exists(p) and os.path.getsize(p) else None


def main():
    os.makedirs(os.path.join(HERE, 'out_table'), exist_ok=True)
    rows = []
    for arr, q, f in CHECKS:
        tag = arr.split(';')[0].replace(',', '-') + '_' + q.replace(',', 'x')
        # quotient_lift.py builds the quotient model from G = (quotient) x Z_f; the model itself only uses
        # the quotient group and f, so any G with that quotient gives the same model.
        ns = q + ',' + str(f)
        cp = run('quotient_lift.py', [arr, ns, '--mode', 'count', '--maxq', '1', '--qtl', '3600',
                                      '--quot', q + ',1'], tag + '_cpsat.json')
        mi = run('quotient_mip.py', [arr, q, str(f), '--tl', '3600'], tag + '_scip.json')
        st_cp = (cp['quotient_status'] if cp['n_quotient'] == 0 else 'FEASIBLE') if cp else 'ERR'
        rows.append(['{' + arr + '}', q, f, st_cp, mi['status'] if mi else 'ERR'])
        print(rows[-1], flush=True)
    with open(os.path.join(HERE, 'quotient_table.csv'), 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['array', 'quotient', 'fibre', 'cpsat', 'scip']); w.writerows(rows)


if __name__ == '__main__':
    main()
