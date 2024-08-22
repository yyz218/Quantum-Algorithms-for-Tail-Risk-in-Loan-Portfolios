# Use to generate random loan portfolio and save the loans in input_file.txt

import random

def random_loans(_):
    # set range for each parameters
    L_i = _
    v_i = Z_i = 1
    P_i = random.randint(10, 150) * 100
    r_i = round(random.uniform(0.26, 0.5), 3)
    gamma_i = round(random.uniform(0.05, 0.2), 3)
    λ_i = round(random.uniform(0.2, 0.3), 2)
    mu_i = round(random.uniform(0.00001, 0.1), 3) 
    sigma_i = round(random.uniform(0.1, 0.5), 3)
    tau_i = round(random.uniform(0.81, 0.88), 2)
    return [L_i, v_i, P_i, r_i, gamma_i, λ_i, mu_i, sigma_i, tau_i, Z_i]

def save_file(filename, num_sets):
    with open(filename, 'w') as file:
        # set the initial market parameters, Z, mu, and sigma
        inp = [1, 0.1, 0.3162]
        line = '{} {} {}\n'.format(*inp)
        file.write(line)
        for _ in range(num_sets):
            params = random_loans(_)
            # Format each parameter for alignment to review in txt file
            line = '{:5} {:1} {:5} {:5} {:5} {:4} {:5} {:5} {:4} {:2}\n'.format(*params)
            file.write(line)

num = input("Input the number of loans needed: ")
save_file('input_file.txt',int(num))
