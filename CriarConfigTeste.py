##########################################################
# Esse codigo cria um arquivo .py com os dados iniciais
# de um teste
##########################################################

nome_arquivo = input("Digite o nome do arquivo de configuração: ")
print(nome_arquivo)

raio = float(input("Digite o raio do objeto (nr type float ex.: 0.1): "))
n_eletrodos = int(input("Digite o Nr de eletrodos (nr type int ex.: 32): "))
ponto_medio_x = str(0)
ponto_medio_y = str(0)
lc1 = float(input("Digite o lc (characteristic length) da cuba (nr type float ex.: 1.0e-1): "))
lc2 = float(input("Digite o lc (characteristic length) da anomalia (nr type float ex.: 1.0e-1): "))

anomalia_raio = raio/4
max_ptoCentral = raio/2
anomalia_lados = int(input("Digite o Nr de lados da anomalia (nr type int ex.: 4): "))
anomalia_rotacao = float(input("digite o ângulo de rotação da anomalia em graus (float ex.:22.5): "))
print("Digite a coordenada x do ponto central da anomalia entre zero e ",  max_ptoCentral)
x_ptoCentral = float(input("Digite o coordenada x (float ex.:0.01): "))
print("Digite a coordenada y do ponto central da anomalia entre zero e ",  max_ptoCentral)
y_ptoCentral =  float(input("digite o coordenada x (float ex.:0.01): "))



arquivo = open('D:/GIT_EIT_2D/EIT_2D/Configs/'+ nome_arquivo  + '.py', 'w+')

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
arquivo.writelines(str(ponto_medio_x))
arquivo.writelines(u'\n')


arquivo.writelines(u'ponto_medio_y = ')
arquivo.writelines(str(ponto_medio_y))
arquivo.writelines(u'\n')


arquivo.writelines(u'lc1 = ')
arquivo.writelines(str(lc1))
arquivo.writelines(u'\n')


arquivo.writelines(u'lc2 = ')
arquivo.writelines(str(lc2))
arquivo.writelines(u'\n')

arquivo.writelines(u'anomalia_lados = ')
arquivo.writelines(str(anomalia_lados))
arquivo.writelines(u'\n')

arquivo.writelines(u'anomalia_rotacao = ')
arquivo.writelines(str(anomalia_rotacao))
arquivo.writelines(u'\n')

arquivo.writelines(u'x_ptoCentral = ')
arquivo.writelines(str(x_ptoCentral))
arquivo.writelines(u'\n')

arquivo.writelines(u'y_ptoCentral = ')
arquivo.writelines(str(y_ptoCentral))
arquivo.writelines(u'\n')

arquivo.close()   