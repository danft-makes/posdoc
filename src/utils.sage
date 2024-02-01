k = QQ
R.<x,y,z> = PolynomialRing(k,3)
M = R.derivation_module()
e = x*M.gen(0) + y*M.gen(1) + z*M.gen(2)
P2 = ProjectiveSpace(R)

def format_for_log(data):
    if isinstance(data, dict):
        return '\n'.join(f"{k}: {format_for_log(v)}" for k, v in data.items())
    elif isinstance(data, list):
        return ', '.join(format_for_log(v) for v in data)
    else:
        return str(data)

def gcd_of_three(a, b, c):
    """Calculate the GCD of three numbers."""
    return gcd(gcd(a, b), c)

def NewDer(dtuple):
    d = dtuple[0]*M.gen(0) + dtuple[1]*M.gen(1) + dtuple[2]*M.gen(2)
    return d

def ToDer0(f,d):
    # M is the derivation module of f and d is in M
    if d(f) in Ideal(f):
        k = d(f)/f
        k = R(k) # coerce from fraction field to polynomial ring, else k*e doesnt work
        newd = d - (k*e)/(f.degree())
        return newd
    else:
        print("d is not a derivation of f")
        return 0

def UnionDer(d,f,g):
    # d is a derivation in Der_0(f)
    gd = g*d - (d(g)*e)/(g.degree()+f.degree())
    return gd

def ZeroesDer(d,g):
    coef = d.list()
    dx = coef[0] ; dy = coef[1] ; dz = coef[2]
    kbar=QQbar
    Rbar.<x,y,z>=PolynomialRing(kbar,3)
    P2bar = ProjectiveSpace(Rbar)
    Igdbar = Ideal(Rbar(dx),Rbar(dy),Rbar(dz))
    Xgdbar = P2bar.subscheme(Igdbar)
    Z = Xgdbar.rational_points()
    gbar = Rbar(g)
    #print(f"The total zeroscheme has length {len(Z)} and is: {Z}. Intersecting with {gbar}...")
    ZcapL = [pt for pt in Z if gbar(list(pt))==0] # I think we dont need to verify this, because Z == ZcapL always #TODO
    #print(f"The zeroscheme has length {len(ZcapL)} and is: {ZcapL}.\n")
    return ZcapL

def ComputeZcapJ(Z,f):
    kbar = QQbar
    Rbar.<x,y,z> = PolynomialRing(kbar,3)
    fbar = Rbar(f)
    P2bar = ProjectiveSpace(Rbar)
    Z = [P2bar(pt) for pt in Z] # it seems we need to enforce Z in P2bar again...
    IJbar = JacobianIdeal(fbar)
    XJbar = P2bar.subscheme(IJbar)
    XJbar_pts = XJbar.rational_points()
    #print(f"Checking if {Z} intersects {XJbar_pts}")
    ZcapJ = [pt for pt in Z if pt in XJbar_pts]
    return ZcapJ


# Canonical derivation of f and g:
def GetPartials(f):
    partials = f.jacobian_ideal().gens()
    return partials
def GenCanonicalDerivationOfPencil(f,g):
    fx, fy, fz = GetPartials(f)
    gx, gy, gz = GetPartials(g)
    d = [fy*gz-gy*fz,fz*gx-fx*gz,fx*gy-fy*gx]
    d = NewDer(d)
    return d
def JacobianIdeal(f):
    partials = GetPartials(f)
    J = Ideal(partials)
    return J
def JacobianSyzygy(f):
    J = JacobianIdeal(f)
    Syz = J.syzygy_module()
    return Syz
def MinimumDegreeColumn(Syz):
    gens = Syz.rows()
    mindegree = min([max(generator[0].degree(),generator[1].degree(),generator[2].degree()) for generator in gens])
    mindeg_generators = [generator for generator in gens if max(generator[0].degree(),generator[1].degree(),generator[2].degree())==mindegree]
    return {"mindegree_generators":mindeg_generators, "degree":mindegree, "len":len(mindeg_generators)}
def MinimumDegreeSyzygy(f):
    Syz = JacobianSyzygy(f)
    mindeg = MinimumDegreeColumn(Syz)
    return mindeg
def process_lines(poly, d, Lines):
    Z_dict = {}
    ZcapJ_dict = {}
    ZcapJCL_dict = {}

    for line in Lines:
        ld = UnionDer(d, poly, line)
        Z = ZeroesDer(ld, line)
        if Z:
            Z_dict[str(line)] = Z
        ZcapJ = ComputeZcapJ(Z, poly)
        if ZcapJ:
            ZcapJ_dict[str(line)] = ZcapJ
        ZcapJCL = ComputeZcapJ(Z, poly * line)
        if ZcapJCL:
            ZcapJCL_dict[str(line)] = ZcapJCL
    return Z_dict, ZcapJ_dict, ZcapJCL_dict

from tqdm import tqdm
import json
def batch_experiment(FC, Lines):
    with open('experiment_data.log', 'w') as log_file:
        for poly in tqdm(FC, desc="Processing Polynomials"):
            experiment_data = {}
            DerDict = MinimumDegreeSyzygy(poly)
            a = DerDict["degree"]
            d = NewDer(DerDict["mindegree_generators"][0])

            Z_dict, ZcapJ_dict, ZcapJCL_dict = process_lines(poly, d, Lines)
            experiment_data[str(poly)] = {
                "a": a,
                "d": d,
                "Z": Z_dict,
                "ZcapJ": ZcapJ_dict,
                "ZcapJCL": ZcapJCL_dict
            }

            # log
            log_file.write(f"Poly: {poly}\n")
            log_file.write(format_for_log(experiment_data))
            log_file.write("\n\n")
    return experiment_data

def single_experiment(poly,line):
    # yanked this function to another one to work on pencil of cubics, TODO refactor this into a single experiment
    d = MinimumDegreeSyzygy(poly)
    a = d["degree"]
    d = NewDer(d["mindegree_generators"][0])
    ld = UnionDer(d, poly, line)
    Zpts = ZeroesDer(ld, line)
    kbar = QQbar
    Rbar.<x,y,z> = PolynomialRing(kbar,3)
    polybar = Rbar(poly)
    for pt in Zpts:
        # check if the point satisfies poly(pt)==0
        insideC = check_pt_inside_C(polybar,pt) 
        print(f"Point {pt} is inside C? {insideC}")
    #ZcapJ = ComputeZcapJ(Z, poly)
    #ZcapJCL = ComputeZcapJ(Z, poly * line)
    ## Check cross-ratio
    return Zpts

def example1(poly,line):
    """
    The purpose of this function is to investigate the crossratio given by distinguished points coming from the free arrangement given by 9 lines and a generic line L.
    We compute the Bourbaki scheme given by the derivation 'ld', which should be always 3 points. We then pick other 1 point defined by the intersection of L and the z = 0 (could be x=0 or y=0...).
    We hope that the crossratio is constant for any generic L (we also assume that L is generic, careful to not feed example1 with nongeneric lines).
    """
    # Example of free 9 lines and a generic line
    C = Curve(poly) ; L = Curve(line)
    # Choose one of the distinguished lines (z=0)
    SpecialLine = P2.line_through(P2(1,0,0),P2(0,1,0))
    P = L.intersection(SpecialLine).rational_points()[0]
    d = MinimumDegreeSyzygy(poly)
    d = NewDer(d["mindegree_generators"][0])
    ld = UnionDer(d, poly, line)
    Zpts = ZeroesDer(ld, line)
    kbar = QQbar
    Rbar.<x,y,z> = PolynomialRing(kbar,3)
    polybar = Rbar(poly)
    for pt in Zpts:
        # check if the point satisfies poly(pt)==0
        insideC = check_pt_inside_C(polybar,pt) 
        print(f"Point {pt} is inside C? {insideC}")
    #ZcapJ = ComputeZcapJ(Z, poly)
    #ZcapJCL = ComputeZcapJ(Z, poly * line)
    ## Check cross-ratio
    CRP = crossratio(P,Zpts[0],Zpts[1],Zpts[2])
    print(f"Cross-Ratio with respect to {P} is {CRP}")
    return [Zpts,CRP]

def example2(poly,line):
    """
    The purpose of this function is to investigate the crossratio given by distinguished points coming from the union of 3 osculating conics and a generic line L.
    We compute the Bourbaki scheme given by the derivation 'ld', which should be always 2 points. We then pick other 2 points defined by the intersection of L and the two special lines x=0, z=0.
    We hope that the crossratio is constant for any generic L (we also assume that L is generic, careful to not feed example2 with nongeneric lines).
    """
    # Example of 3 osculating conics and a generic line
    C = Curve(poly) ; L = Curve(line)
    # Define singular points
    P = P2(0,1,0) ; Q = P2(1,0,0)
    Ptangents = C.tangents(P) ; Qtangents = C.tangents(Q)
    # The 'distinguished' lines are given by x=0, the unique tangent of P, and z=0, the line containing P and Q. We find the intersections with 'line' now
    uniqueTangent = Ptangents[0] ; uniqueTangent = Curve(uniqueTangent)
    lineThroughSingularities = P2.line_through(P,Q)
    L_x = L.intersection(uniqueTangent).rational_points()[0]
    L_z = L.intersection(lineThroughSingularities).rational_points()[0]
    print(L_x,L_z)
    d = MinimumDegreeSyzygy(poly)
    d = NewDer(d["mindegree_generators"][0])
    ld = UnionDer(d, poly, line)
    Zpts = ZeroesDer(ld, line)
    kbar = QQbar
    Rbar.<x,y,z> = PolynomialRing(kbar,3)
    polybar = Rbar(poly)
    for pt in Zpts:
        # check if the point satisfies poly(pt)==0
        insideC = check_pt_inside_C(polybar,pt) 
        print(f"Point {pt} is inside C? {insideC}")
    #ZcapJ = ComputeZcapJ(Z, poly)
    #ZcapJCL = ComputeZcapJ(Z, poly * line)
    ## Check cross-ratio
    CR = crossratio(L_z,L_x,Zpts[1],Zpts[0])
    print(f"Cross-Ratio is {CR}")
    return [Zpts,CR]

def crossratio(A,B,C,D):
    # Points P are of the form (a:b:c), we use the trick of projecting to z = 0 from (0:0:1) to justify that the projection becomes (a:b:0), so we just ignore the last coordinate and we think we are in P1
    AC = detcr(A,C) ; BC = detcr(B,C) ; BD = detcr(B,D) ; AD = detcr(A,D)
    CR = (AC/BC)*(BD/AD)
    return CR

def detcr(P,Q):
    P = list(P) ; Q = list(Q)
    det = P[0]*Q[1]-Q[0]*P[1]
    return det

def check_pt_inside_C(poly,pt):
    insideC = poly(list(pt))==0
    return insideC

def eigenscheme(d):
    """
    Computes the eigenscheme defined by the derivation d
    """
    kbar = QQbar
    Rbar.<x,y,z> = PolynomialRing(kbar,3)
    P2bar = ProjectiveSpace(Rbar)
    dlist = d.list() ; d1 = dlist[0] ; d2 = dlist[1] ; d3 = dlist[2]
    R1 = x*d2 - y*d1
    R2 = x*d3 - z*d1
    R3 = y*d3 - z*d2
    I = ideal(Rbar(R1),Rbar(R2),Rbar(R3))
    X = P2bar.subscheme(I)
    return X.rational_points()


def splitting(poly,line):
    # For now this is only verifying if a the points are ordinary
    kbar = QQbar
    Rbar.<x,y,z> = PolynomialRing(kbar,3)
    polybar = Rbar(poly) ; linebar = Rbar(line)
    C = Curve(polybar)
    L = Curve(linebar)
    CcapL = C.intersection(L) ; pts = CcapL.rational_points()
    CcupL = Curve(polybar*linebar)
    for pt in pts:
        ordinary = CcupL.is_ordinary_singularity(pt)
        print(f"Point {pt} is ordinary? {ordinary}")
    return pts
