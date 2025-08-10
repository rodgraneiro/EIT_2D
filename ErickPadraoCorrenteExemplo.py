

n_eletrodos = 8
pula =2

for i in range(n_eletrodos):
    elet_in = i + 1
    elet_out = (i+(pula+1)) % n_eletrodos +1
    print(f'{elet_in}; {elet_out}')