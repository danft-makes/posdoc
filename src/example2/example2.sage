load('../utils.sage')
C1 = x*z ; C2 = z^2-x*y ; f = (C2)*(C2+C1)*(C2-C1) ; C = Curve(f)

# Load the variable back
Lines = []
with open("generic_lines_ex2.txt", "r") as file:
    for line in file:
        # Convert each line back into a tuple
        Lines.append(tuple(map(int, line.strip().split(','))))

with open('example2.log', 'w') as log_file:
    for a,b,c in tqdm(Lines):
        #print(a,b,c)
        l = a*x+b*y+c*z
        Zpts, CR = example2(f,l)
        log_file.write(f"{(a,b,c)}\n")
        log_file.write(f"Cross-Ratio: {CR}\n")

