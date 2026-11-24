# Incidence directions and five-line star partitions

The no-extension system in [theta_line_margins.md](theta_line_margins.md)
uses point and line distance tables around the 32 named incidence vertices.
There are further restrictions on its line types, on the individual
five-line star partitions, and on the integrality of orbit totals. All
can be computed without knowing a completion or choosing automorphisms
of the ambient octagon.

## Directions on a prospective line

Use incidence distance throughout this section. For a prospective line
$L$ and a named vertex $v$, let $k=d(L,v)$. If $k=0$, its three points have
distance one. If $k=8$, all three have distance seven. If $0<k<8$, exactly
one has distance $k-1$ and the other two have distance $k+1$. This follows
from the unique geodesic below eight and the girth and diameter. The rule
must hold also when $v$ is one of the eleven named lines. It eliminates
2,390 of the 8,982 line types from the previous system.

Here is another necessary rule, for a root vertex $x$ and two different
named vertices $u,v$ at distances $h,k$ from $x$, both in $1,\ldots,7$.
Their paths from $x$ are unique. Since their combined length is at most
14, they share an initial segment and cannot separate and meet again.
If their first edges are different, the path via $x$ has length $h+k$.
If their first edge is the same but their second edges are different,
it has length $h+k-2$. In either case the resulting nonbacktracking path
of length $a$ implies

$$d(u,v)=a\quad(a\le8),\qquad d(u,v)\ge16-a\quad(a>8).$$

If the first two edges agree, only the upper bound $d(u,v)\le h+k-4$
is used. If one of the two paths ends at their common first edge, then
$d(u,v)=h+k-2$ instead. When only agreement of the first edge is known,
the upper bound is $h+k-2$.

Apply these tests both at $L$ and at each of its points. At $L$ the first
point of every short geodesic is known from its profile. At a point $x$
on $L$ we know whether the first line towards a named vertex is $L$;
when it is and the distance is at least two, the profile of the other
two points of $L$ identifies the next point. For a pair of named targets
whose geodesics both start on another line, this one line gives no
restriction. All comparisons use the given 32-vertex distance matrix.
Another 102 line types are eliminated, leaving 6,490 in 596 orbits.

## Star partitions

For a point profile $p$, include all 32 coordinates at incidence distance
one through seven in a bit mask $M_p$. Every incident line must account
for the coordinates for which it is the unique line one step closer.
The five inward masks must partition $M_p$, including possible repeated
empty masks. Each allowed line type gives its mask at a place carrying
profile $p$. The [filter](theta_incidence_stars.py) checks that a five-mask
partition is possible; none of the 753 point profiles left by the earlier
tests are removed.

For each of the 78 point profile orbits choose its first representative
in the profile enumeration. Enumerate every five-mask partition using
the allowed masks at that profile. A new nonnegative multiplicity for
each of these partitions counts the points of the representative profile
having that partition. The sum of these multiplicities is the point
profile multiplicity. For each possible mask, count the line places with
that mask and require that their number equals the total number of slots
with that mask in these partitions. This is stronger than demanding only
the aggregate inward count for each coordinate separately.

The equations can again be averaged over the twelve relabellings of the
seed. At a representative the star partition multiplicities may be
averaged over its stabilizer; relabelling then supplies the other profiles
of its orbit. No integral symmetry assumption is made. There are 3,139
five-mask partitions at the chosen representatives, producing 1,243
additional rows. The complete rational system has 2,625 rows in 3,813
columns, including the previous point and line marginal equations on
the reduced support.

## Necessary integrality congruences

Let $T_i$ and $S_j$ be **total** multiplicities in point profile orbit $i$
and line type orbit $j$, respectively. They must be integers in an actual
completion even if the geometry has no symmetries. The common averaged
multiplicity of a named profile is its total divided by its orbit size.
Every orbit size divides twelve. Replacing the averaged multiplicities
by these quotients and multiplying by twelve therefore gives an integer
coefficient system in $T,S$.

In a mixed degree or inward-count row there is just one point coefficient,
say $aT_i+\sum_j b_jS_j=0$. If $g=\gcd_j b_j$, the total $T_i$ must be
divisible by $g/\gcd(a,g)$. Taking the least common multiple over the
relevant rows gives the following nontrivial divisors in zero-based orbit
order (the same order as the certificate).

| Point orbit | Divisor |
|---|---:|
| 0 | 3 |
| 1, 3, 6, 9, 10, 21, 22, 23 | 2 |

For stronger congruences write $T_i=v_i U_i$, with these divisors $v_i$
and $v_i=1$ elsewhere. Divide each row with line coefficients by their
gcd (which now divides all coefficients and the right-hand side), and
divide each point-only row by the gcd of its point coefficients. Reduce
modulo two and eliminate the line variables first. There remain 52
independent parity equations on the $U_i$. Doing the same reduction on
point-only rows gives 43, so nine more independent parity equations are
obtained from the requirement of integral line orbit totals. The
[verifier](theta_divisibility.py) does this elimination by XOR on integer
bit sets, not by floating-point rank calculation.

The star-partition variables themselves do not enter these congruences.
An integral completion would still have to supply actual line orbit
totals, not merely satisfy these necessary congruences.

## Exact witness for the strengthened relaxation

The [rational certificate](theta_star_margins.json) satisfies the full
star-partition system, has **integral point orbit totals** adding to 1,755,
and satisfies all individual divisibility conditions and all 52 parity
equations. It was found by a mixed integer search for the point totals,
followed by exact rational reconstruction with those totals fixed. It
lists 867 positive numerators and a common denominator of 17 digits.
Verification uses only the enumerations above, modular elimination and
integer multiplication of the equation rows by the numerators; it does
not use the numerical solver.

```text
$ python thin/theta_star_margins.py
incidence line rule 6592 of 8982
incidence compatible lines 6490 of 6592
remaining points 753 78 lines 6490 596
profiles 753 78 line types 6490 596
matrix (1382, 674) 108474
star patterns 3139 new rows 1243
full matrix (2625, 3813) 128936
exact rational star witness 867 occupied columns
integer point orbit totals, all individual divisibility conditions satisfied
failed parity equations 0 of 52
```

This is a fractional assignment of line types and star partitions, not
an integral assignment of all profiles and lines. It also does not impose
the remaining distances and cycles of an incidence graph. Thus neither
the strengthened star-partition equations nor these congruences force a
seventh point. Whether a no-extension seed can satisfy full integrality,
let alone occur in an arbitrary octagon, is still open.

[Verifier](theta_star_margins.py)
