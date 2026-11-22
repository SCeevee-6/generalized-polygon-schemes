# The 60 extra triple relations of GO(2,4): exact form, and which ones hold in every GO(2,4)

## Background

`triple_note.md` found linear relations among the triple intersection numbers [ijk] of the Ree-Tits octagon O(2) that do not follow from the intersection numbers and the Krein equation q^4_44 = 0:

| type | relations |
|---|---:|
| (2,3,3) | 15 |
| (2,3,4) | 12 |
| (2,4,4) | 9 |
| (3,3,3) | 14 |
| (3,3,4) | 7 |
| (3,4,4) | 3 |
| total | 60 |

Here the type (U,V,W) is (d(u,v), d(u,w), d(v,w)), and [ijk] counts the points x with d(u,x) = i, d(v,x) = j, d(w,x) = k. This note writes the relations down exactly, then sorts them into those I can prove for every GO(2,4) and those I can't.

## 1. Exact form (`triple_relations.py`)

For each type, the method is:
- The LP system is rebuilt over Q.
- Columns with an index 4 are put first, so that the RREF pivots on them. The free variables are then all counts [ijk] with i, j, k ≤ 3, i.e. points within distance 3 of all of u, v, w, which have unique geodesics to them. The counts with an index 4 are affine in these and the intersection numbers.
- The relations are an integer basis (RREF, primitive rows) of the linear forms on the free variables that are constant on all T-orbit vectors.

Every relation is checked exactly on every T-orbit vector.

| type | free dim | T-orbits | distinct vectors | affine span | relations | max abs. coefficient |
|---|---:|---:|---:|---:|---:|---:|
| (2,3,3) | 17 | 6 | 4 | 2 | 15 | 4 |
| (2,3,4) | 13 | 3 | 2 | 1 | 12 | 3 |
| (2,4,4) | 11 | 8 | 4 | 2 | 9 | 6 |
| (3,3,3) | 19 | 18 | 6 | 5 | 14 | 50 |
| (3,3,4) | 15 | 28 | 14 | 8 | 7 | 8 |
| (3,4,4) | 12 | 40 | 18 | 9 | 3 | 7 |

Files:
- `out_orbits/triple_relations.csv`: one equation per row, e.g. `+[212] = 0`.
- `out_orbits/triple_relations_long.csv`: the integer coefficient matrices in long format (type, relation, unknown, coefficient, plus a `rhs` row).

## 2. Relations that hold in every GO(2,4)

Fourteen of the 60 are proved below.

**Zeros: a mechanical test (`triple_zero_proofs.py`).** Let x be within distance 3 of u, v and w.
- The geodesics from x to u, v, w are unique, and their union T_x is a tree: two of them that separate and meet again would close a circuit of length ≤ 12 in the incidence graph, which has girth 16.
- T_x is fixed by the three Gromov products, i.e. the lengths of the common initial segments. For a tree, the two smallest of these are equal.
- The path in T_x between two of u, v, w is non-backtracking. If its length is at most 8, it is a geodesic and gives the distance exactly. If it is longer than 8, it closes a circuit with a geodesic, so the two lengths add up to at least 16.
- At each branch node, the incoming edge plus the outgoing ones cannot exceed the degree (s+1 = 3 for a line, t+1 = 5 for a point).

If no choice of Gromov products passes all of this, then [ijk] = 0 for that type in every GO(2,4). The test never contradicts the Ree-Tits data, and in the six types it proves every zero that the Ree-Tits octagon shows:

| type | zeros proved for every GO(2,4) |
|---|---|
| (2,3,3) | [113], [122], [212], [231], [232], [321], [322] |
| (2,3,4) | [123], [222] |
| (3,3,3) | [222] |
| (3,3,4) | [222] |

**Two counts equal to t, by hand.**

* *Type (2,3,4): [313] = t = 4.* Here d(u,v) = 2 and d(v,w) = 4.
  - Let p be the middle point of u–v and L = pv. On L, p is collinear with u and the other points are at distance 2 from u.
  - On each of the other t lines through v, all s points are at distance 3 from u, by (G1).
  - w is opposite v, so each line M on v is at distance 7 from w in the incidence graph. M therefore has exactly one point at distance 3 from w, and that point is not v.

  This gives exactly t points x with d(u,x) = 3, d(v,x) = 1, d(w,x) = 3. The same argument gives [313] = t = 4 in type (2,4,4).
* *Type (2,4,4): [133] = t = 4.* Here d(u,v) = 2 and d(u,w) = 4. This is the same argument with u in place of v: the points collinear with u off the line toward v are at distance 3 from v, and each of the t lines on u off that line has exactly one point at distance 3 from w.

Adding these constants to the LP lowers its free dimension by

| type | drop | relations |
|---|---:|---:|
| (2,3,3) | 7 | 15 |
| (2,3,4) | 3 | 12 |
| (2,4,4) | 2 | 9 |
| (3,3,3) | 1 | 14 |
| (3,3,4) | 1 | 7 |
| (3,4,4) | 0 | 3 |

That makes **14 of the 60 relations proved for every GO(2,4)**. One further constant holds on every Ree-Tits orbit, [332] = 16 in type (2,3,4), and is not implied by the proved ones. I have not proved it.

## 3. The other 46: where the octagon's own structure comes in (`triple_shapes.py`)

**Labels.** The *base configuration* of a triple is the union of the geodesics (in the incidence graph) between its non-opposite pairs. Each x within distance 3 of all three points is labelled by what its geodesics add to the base configuration:
- *tree*: no new circuit;
- *c16*: new circuits of length 16 only (apartments);
- *long*: a new circuit of length more than 16 is needed.

Distinct vectors per base shape (orbits, distinct vectors):

| type | base shapes |
|---|---|
| (2,3,3) | apartment (16-circuit): 2 orbits, 1 vector; tripod: 4 orbits, **3** vectors |
| (2,3,4) | two shapes: 2 orbits, 1 vector and 1 orbit, 1 vector |
| (2,4,4) | one shape (the base is the u–v geodesic plus the isolated point w): 8 orbits, 4 vectors |
| (3,3,3) | tripod 6/1; 18-circuit 6/**2**; three different 16-circuit shapes 2/1 each |
| (3,3,4) | 2/1, 20/**7**, 6/**6** |
| (3,4,4) | one shape, 40 orbits, 18 vectors |

So in (2,3,3), (3,3,3) and (3,3,4) there are Ree-Tits triples with the same base configuration and different vectors. In the tripod case of (2,3,3) the tree part is the same for all orbits. The difference is entirely in the c16 and long parts, i.e. in how many x close one or two further apartments or longer circuits. The counts therefore depend on more than the local configuration. Any argument that only looks at the geodesics spanned by u, v, w and x (like the tree test) cannot give these relations.

For types with an opposite pair the base shape doesn't see the opposite pair at all, so a finer invariant would be needed there. I haven't tried one.

Status of the remaining 46:
- They hold in the Ree-Tits octagon.
- They are not consequences of the LP, the Krein equation or the local arguments above.
- I can't say whether they are specific to the Ree-Tits octagon. It is the only known GO(2,4), so there is nothing to compare against.

**The (2,2,4) vs (2,3,3) explanation.** The earlier guess was that (2,3,3) has several vectors because its geodesics can close a circuit, unlike (2,2,4). That doesn't hold up. The (2,3,3) triples whose geodesics form an apartment all share one vector. The extra vectors come from the tripod triples, which have no circuit at all.

## Files

- `triple_relations.py` → `out_orbits/triple_relations.csv`, `triple_relations_long.csv`, `triple_relations_output.txt`
- `triple_zero_proofs.py` → `out_orbits/triple_zero_proofs_output.txt` (all 17 types; usage `python triple/triple_zero_proofs.py [s t]`)
- `triple_shapes.py` → `out_orbits/triple_shapes_output.txt`, `triple_shapes_<type>.csv` (per orbit, class and [ijk])
