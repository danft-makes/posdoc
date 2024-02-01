load('../utils.sage')
f = x*y*z*(x^2-y^2)*(y^2-z^2)*(x^2-z^2) ; C = Curve(f)

# Load the variable back
Lines = []
with open("generic_lines_ex1.txt", "r") as file:
    for line in file:
        # Convert each line back into a tuple
        Lines.append(tuple(map(int, line.strip().split(','))))

with open('example1.log', 'w') as log_file:
    for a,b,c in tqdm(Lines):
        #print(a,b,c)
        l = a*x+b*y+c*z
        Zpts, CRP = example1(f,l)
        log_file.write(f"{(a,b,c)}\n")
        log_file.write(f"Cross-Ratio: CRP = {CRP}\n")


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
