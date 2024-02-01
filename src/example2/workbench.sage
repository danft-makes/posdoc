load('../utils.sage')
C1 = x*z ; C2 = z^2-x*y ; f = (C2)*(C2+C1) ; C = Curve(f)

# tests
d = MinimumDegreeSyzygy(f)['mindegree_generators'][0]
print(d)

#kbar = QQbar ; Rbar.<x,y,z> = PolynomialRing(kbar,3) ; P2bar = ProjectiveSpace(Rbar)
#IJ = JacobianIdeal(Rbar(f)) ; J = P2.subscheme(IJ) ; Zsection = J.rational_points()
e1tangent = C.tangents(P2(1,0,0))
e2tangent = C.tangents(P2(0,1,0))

print(e1tangent,e2tangent)
