# Smith factor 2 for opposition matrices

`smith_two.py` computes, for the point and line opposition matrices of O(2), the binary rank r0 and the number e1 of
invariant factors exactly divisible by 2. It uses the map x ↦ A·lift(x)/2 mod 2 from ker(A mod 2) to
F_2^n / im(A mod 2); this map is well defined (changing the lift by 2z changes the image by A z) and commutes with
the group, and its rank is e1. Results: points r0 = 26, e1 = 272; lines r0 = 246, e1 = 350, matching the earlier
full Smith forms.

So the image of the map is a 272-dimensional F_2 T-module for points and a 350-dimensional one for lines. For O(8)
the matrix has order about 1.2·10^9 and nothing here reaches it. What would be needed is the composition structure of
the 2-modular permutation module of 2F4(8) on points (and on lines) restricted to ker and coker of A mod 2. I have
not checked the literature for it. One data point (q = 2) does not determine a formula: 272 = 26 + 246 and
350 = 246 + 104 are consistent with several guesses in terms of the 2-modular irreducibles.
