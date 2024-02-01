load('../utils.sage')
f = x*y*z*(x^2-y^2)*(y^2-z^2)*(x^2-z^2) ; C = Curve(f)

"""
How to generate generic lines?
L is generic to C if for every intersection point P we have that the intersection multiplicity is minimal?
"""
limit = 5
triples = [(a, b, c) for a in range(1, limit + 1)
                    for b in range(1, limit + 1)
                    for c in range(1, limit + 1)
                    if gcd_of_three(a, b, c) == 1]

Lines = triples
kbar = QQbar
Rbar.<x,y,z> = PolynomialRing(kbar,3)
P2bar=ProjectiveSpace(Rbar)
GenericLines = []
from tqdm import tqdm
for a,b,c in tqdm(Lines, desc="Generating Lines"):
    l = a*x+b*y+c*z ; Lbar=Curve(Rbar(l),P2bar) ; Cbar = Curve(Rbar(f),P2bar)
    CcupL = Curve(Rbar(f)*Rbar(l),P2bar)
    CcapL = Cbar.intersection(Lbar)
    cap_pts = CcapL.rational_points()
    culprits = 0
    for pt in cap_pts:
        # Check if pt is ordinary in C\cupL
        is_ordinary = CcupL.is_ordinary_singularity(pt)
         # Check if the point was a singularity before
        tangents = Cbar.tangents(pt) ; number_tangents = len(tangents)
        if (not is_ordinary):
            print(f"{(a,b,c)}: Point {pt} is not ordinary singularity in CcupL")
            culprits += 1
        if (number_tangents>1):
            print(f"{(a,b,c)}: Point {pt} has {number_tangents} tangents in C")
            culprits += 1
    if culprits==0:
        GenericLines.append((a,b,c))

with open("generic_lines_ex1.txt", "w") as file:
    for tup in GenericLines:
        file.write(','.join(map(str, tup)) + '\n')


