load('utils.sage')
# Bag of free curves (#TODO make a function to fill the bag, for now we use handpicked ones)
C1 = x*z ; C2 = z^2-x*y
FC = [(C1)*(C1+C2)*(C1-C2),(C2)*(C2+C1)*(C2-C1),x*y*z*(x^2-y^2)*(y^2-z^2)*(x^2-z^2)]

h = FC[1]
a,b,c = (1,2,3) # test
l = a*x+b*y+c*z
#pts = splitting(h,l)
Zpts = single_experiment(h,l)

# CHECK CROSS RATIO OF PARTICULAR EXPERIMENT
kbar = QQbar
Rbar.<x,y,z> = PolynomialRing(kbar,3)
P2bar=ProjectiveSpace(Rbar)
fourpoints = [P2bar(Zpts[0]),P2bar(Zpts[1]),P2bar([-2,1,0]),P2bar([0,1,0])]
CR = crossratio(fourpoints[0],fourpoints[1],fourpoints[2],fourpoints[3])
print(CR)
