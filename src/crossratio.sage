load('utils.sage')
# Bag of free curves (#TODO make a function to fill the bag, for now we use handpicked ones)
C1 = x*z ; C2 = z^2-x*y ; f = (C2)*(C2+C1)*(C2-C1) ; C = Curve(f)


#Lines = [(1,2,1),(1,2,3),(1,2,4),(1,2,5),(1,3,1),(1,3,5),(2,1,3),(2,1,5),(2,3,5)]

#for a,b,c in Lines:
#    print(a,b,c)
#    l = a*x+b*y+c*z
#    Zpts = example2(f,l)

#a,b,c = (1,2,3) # test
#l = a*x+b*y+c*z
#Zpts = example2(f,l)

# CHECK CROSS RATIO OF PARTICULAR EXPERIMENT
#kbar = QQbar
#Rbar.<x,y,z> = PolynomialRing(kbar,3)
#P2bar=ProjectiveSpace(Rbar)
#fourpoints = [P2bar(Zpts[0]),P2bar(Zpts[1]),P2bar([-2,1,0]),P2bar([0,1,0])]
#CR = crossratio(fourpoints[0],fourpoints[1],fourpoints[2],fourpoints[3])
#print(CR)

"""
How to generate generic lines?
L is generic to C if for every intersection point P we have that the intersection multiplicity is minimal?
"""
limit = 3
triples = [(a, b, c) for a in range(1, limit + 1)
                    for b in range(1, limit + 1)
                    for c in range(1, limit + 1)
                    if gcd_of_three(a, b, c) == 1]

#Lines = [triple[0]*x+triple[1]*y+triple[2]*z for triple in triples]
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
    for pt in cap_pts:
        # Check if pt is ordinary in C
        is_ordinary = CcupL.is_ordinary_singularity(pt)
        #print(f"Point {pt} is ordinary? {is_ordinary}")
        # Check if L is a tangent line to pt
        tangents = Cbar.tangents(pt)
        is_tangent = l in tangents
        #print(f"C has tangents {tangents}, is L tangent? {is_tangent}")
        if (is_ordinary) and (is_tangent):
            print(f"{tuple(a,b,c)} is generic")
            GenericLines.append(tuple(a,b,c))

with open("generic_lines.txt", "w") as file:
    for tup in GenericLines:
        file.write(','.join(map(str, tup)) + '\n')


