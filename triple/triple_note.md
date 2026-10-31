# Triple intersection numbers for GO(2,4) and the vanishing Krein parameter (angle 5)

Labels: **known**, **check** (computational confirmation of something known).

## Literature

- Uniqueness of GO(2,4) is open. Partial results: De Bruyn (the Tits octagon is the only GO(2,4) containing a
  sub-GO(2,1)); Cohen, O'Brien, Shpectorov, J. Algebra (2014) (unique under a condition on a point stabiliser); and
  uniqueness among those with a vertex-transitive automorphism group (cited in the same paper). **known**
- q^4_44 = 0 for extremal generalized octagons, and the triple-intersection-number method of Coolsaet and Jurišić for
  vanishing Krein parameters, are both in the literature (De Bruyn, EJC 29(4) P4.46 (2022), item (4) of the
  introduction; van Dam–Koolen–Tanaka survey §6.3.1). Vidali's `sage-drg` automates the method. **known**

## Computation

`go24_triple_lp.py` sets up, for each of the 17 triple types (d(u,v), d(u,w), d(v,w)), the linear system for the
numbers [i j k] (marginals equal to intersection numbers, plus the Krein equation
sum Q_{i4} Q_{j4} Q_{k4} [i j k] = 0), and computes the range of every [i j k] by LP, with and without the Krein
equation. It then counts the same numbers for one triple of each type in the actual octagon (ATLAS 1755-point action).

## Result

- The Krein equation cuts the solution space by exactly one dimension for every type, leaving 2 to 19 free
  parameters.
- It fixes values that are otherwise free only for the type (1,1,2) (a path v - u - w with v, w not collinear):
  [3 3 3] = 0, [3 3 4] = [3 4 3] = 64, [3 4 4] = 384, and [4 3 3] = [4 3 4] = [4 4 3] = [4 4 4] = 256. These are
  exactly what the near-polygon property of the octagon gives (a point at distance 3 from u reaches u through one
  line only), and the octagon itself has these values. **check**
- For all other types the Krein equation only narrows the LP ranges; the octagon's values lie inside them.

## Assessment

At the level of triple intersection numbers the vanishing Krein parameter reproduces the octagon geometry for one
triple type and gives no further restriction. This does not move the uniqueness question. Going further would need
quadruple counts or a different idea; I stopped here.
