"""Class representatives of T = 2F4(2)' as words in the ATLAS standard generators a, b.

Words copied from the ATLAS v3 page for 2F4(2)' (section "Conjugacy classes"), first
alternative of each class. A digit after a letter, ')' or ']' is an exponent.
[x, y] = x^-1 y^-1 x y.
"""
ATLAS_CLASS_WORDS = [
    ("1A", "a2"), ("2A", "a"), ("2B", "[a,bab]2"), ("3A", "b"),
    ("4A", "(ababab2)3"), ("4B", "ababababab2abab2ab2"), ("4C", "ababab2ab2"),
    ("5A", "abab2"), ("6A", "ababababab2"), ("8A", "(ab)5ab2(ababab2)2ab2"),
    ("8B", "(ab)5ab2ab2(ababab2)2"), ("8C", "abababab2"),
    ("8D", "abab(abababab2)2"), ("10A", "abababab2abab2ab2"), ("12A", "ababab2"),
    ("12B", "abababab2ab2"), ("13A", "ab"), ("13B", "(ab)2"),
    ("16A", "(ab)5ab2ab(ab2)3"), ("16B", "(ab)5(ab2)3abab2"),
    ("16C", "(ab)3(ab2)3ababab2"), ("16D", "abababab2abab(ab2)3"),
]


def _parse(s, i=0):
    """Return (list of atoms, next index); atom = ('g', name) | ('seq', atoms) | ('comm', x, y), with exponent."""
    out = []
    while i < len(s) and s[i] not in ")],":
        c = s[i]
        if c in "ab":
            atom = ("g", c); i += 1
        elif c == "(":
            inner, i = _parse(s, i + 1); assert s[i] == ")"; i += 1
            atom = ("seq", inner)
        elif c == "[":
            x, i = _parse(s, i + 1); assert s[i] == ","; y, i = _parse(s, i + 1); assert s[i] == "]"; i += 1
            atom = ("comm", x, y)
        else:
            raise ValueError("bad char %r in %r" % (c, s))
        j = i
        while j < len(s) and s[j].isdigit():
            j += 1
        e = int(s[i:j]) if j > i else 1
        i = j
        out.append((atom, e))
    return out, i


def parse_word(s):
    w, i = _parse(s.replace(" ", ""))
    assert i == len(s.replace(" ", "")), s
    return w


def evaluate(word, a, b, mul, inv, one):
    """Evaluate a parsed word with user-supplied group operations."""
    def ev(seq):
        r = one
        for atom, e in seq:
            if atom[0] == "g":
                x = a if atom[1] == "a" else b
            elif atom[0] == "seq":
                x = ev(atom[1])
            else:
                x_, y_ = ev(atom[1]), ev(atom[2])
                x = mul(mul(mul(inv(x_), inv(y_)), x_), y_)
            for _ in range(e):
                r = mul(r, x)
        return r
    return ev(word)


def to_gap(s):
    """GAP expression for a word, e.g. 'abab2' -> 'a*b*a*b^2'."""
    def g(seq):
        parts = []
        for atom, e in seq:
            if atom[0] == "g":
                t = atom[1]
            elif atom[0] == "seq":
                t = "(" + g(atom[1]) + ")"
            else:
                t = "Comm(" + g(atom[1]) + "," + g(atom[2]) + ")"
            parts.append(t if e == 1 else t + "^" + str(e))
        return "*".join(parts) if parts else "One(a)"
    return g(parse_word(s))
