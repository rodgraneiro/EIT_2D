##########################################################
# Esse codigo cria um arquivo .py com os dados iniciais
# de um teste
##########################################################

import os

nome_arquivo = input("digite o nome do arquivo: ")
print(nome_arquivo)

raio = float(input("digite o raio do objeto (nr type float ex.: 0.1): "))
n_eletrodos = int(input("digite o Nr de eletrodos (nr type int ex.: 32): "))
ponto_medio_x = str(0)
ponto_medio_y = str(0)
lc1 = float(input("digite o lc (characteristic length) da cuba (nr type float ex.: 1.0e-1): "))
lc2 = float(input("digite o lc (characteristic length) da anomalia (nr type float ex.: 1.0e-1): "))

#arquivo = open(nome_arquivo  + '.py', 'w+')
arquivo = open(os.path.join(".\\Configs\\", nome_arquivo + ".py"), "w+")
arquivo.writelines(u'nome_arquivo = ')
arquivo.writelines(u'"')
arquivo.writelines(nome_arquivo)
arquivo.writelines(u'" \n')
arquivo.writelines(u'raio = ')
arquivo.writelines(str(raio))
arquivo.writelines(u'\n')

arquivo.writelines(u'n_eletrodos = ')
arquivo.writelines(str(n_eletrodos))
arquivo.writelines(u'\n')

arquivo.writelines(u'ponto_medio_x = ')
arquivo.writelines(ponto_medio_x)
arquivo.writelines(u'\n')


arquivo.writelines(u'ponto_medio_y = ')
arquivo.writelines(ponto_medio_y)
arquivo.writelines(u'\n')

arquivo.writelines(u'lc1 = ')
arquivo.writelines(str(lc1))
arquivo.writelines(u'\n')

arquivo.writelines(u'lc2 = ')
arquivo.writelines(str(lc2))

arquivo.close()   