load('utils.sage')
# Bag of free curves (#TODO make a function to fill the bag, for now we use handpicked ones)
C1 = x*z ; C2 = z^2-x*y
FC = [(C1)*(C1+C2)*(C1-C2),(C2)*(C2+C1)*(C2-C1),x*y*z*(x^2-y^2)*(y^2-z^2)*(x^2-z^2)]

limit = 5
#triples = [(a, b, c) for a in range(-limit, limit + 1)
#                    for b in range(-limit, limit + 1)
#                    for c in range(-limit, limit + 1)
#                    if gcd_of_three(abs(a), abs(b), abs(c)) == 1 and (a, b, c) != (0, 0, 0)]
triples = [(a, b, c) for a in range(1, limit + 1)
                    for b in range(1, limit + 1)
                    for c in range(1, limit + 1)
                    if gcd_of_three(a, b, c) == 1]

Lines = [triple[0]*x+triple[1]*y+triple[2]*z for triple in triples]

experiment_data = compute_experiment_data(FC, Lines)
