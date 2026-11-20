"""Run abelian_cayley.py over jobs.json; one JSON per job in out/, skipped if present."""
import json, os, subprocess, sys, time
here = os.path.dirname(os.path.abspath(__file__))
tl = sys.argv[1] if len(sys.argv) > 1 else '600'
maxat = int(sys.argv[2]) if len(sys.argv) > 2 else 10**9
for arr, ns, na, elem in json.load(open(os.path.join(here, 'jobs.json'))):
    if na > maxat: continue
    tag = arr.split(';')[0].replace(',', '-') + '_' + ns.replace(',', 'x')
    out = os.path.join('out', tag + '.json')
    if os.path.exists(os.path.join(here, out)): continue
    t = time.time()
    with open(os.path.join(here, 'out', tag + '.err'), 'w') as fe:
        subprocess.run([sys.executable, 'abelian_cayley.py', arr, ns, '--limit', '200', '--tl', tl,
                        '--workers', '8', '--out', out], cwd=here, stdout=fe, stderr=fe)
    with open(os.path.join(here, 'run_all.log'), 'a') as f:
        f.write(f'{tag} {na} {round(time.time() - t)}s\n')
