load('../utils.sage')
C1 = x*y*z ; C2 = x^2+y^2+z^2 ; f = C1^2+C2^3 ; C = Curve(f)

d = MinimumDegreeSyzygy(f)['mindegree_generators'][0]
d = NewDer(d)
epts = eigenscheme(d)
print(epts)
a,b,c = (1,2,3) ; l = a*x+b*y+c*z
ld = UnionDer(d,f,l)
print(ld)
lepts = eigenscheme(ld)
print(lpets)
