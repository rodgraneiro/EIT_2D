# EIT_2D


Este código, ainda em fase de desenvolvimento,  resolve o 'Problema Direto' e o 'Problema Inverso' para determinar a condutividade  $\sigma$ discretizada numa lalha de Elementos Finitos 2D, por meio da minimização pelo método 'Gauss-Newton' da equação: 



$$
\hat{\sigma}_{k+1} = \hat{\sigma}_k + \alpha_k
\left( J_k^T W_1 J_k + \lambda^2 L_2^T L_2 \right)^{-1}
\cdot \left( J_k^T W_1 (z - h(\hat{\sigma}_k)) - \lambda^2 L_2^T L_2 (\hat{\sigma}_k - \sigma^*) \right)
$$


onde,
- $\sigma$ é condutividade $S/m$;
  - $\hat{\sigma}_{k+1}$ condutividade estimada da próxima iteração:
  - $\hat{\sigma}_k$ condutividade estimada atual;
  - $\sigma^*$ condutividade 'a priori'
- $\alpha$ é uma constante de regularização;
- $\lambda$ é uma constante de regularização;
- $h$ é a matriz de condutividades (rigidez);
- $W_1$ é a matriz Identidade;
- $J$ é a matriz jacobiana;
- $L_2$ é a matriz de um Filtro Passa-Alta.

---

Para rodar esse código em python são necessários as seguintes bibliotecas;
- numpy;
- matplotlib;
- meshio;
- time;
- PyQtgraph;
- QtCore;
- PyQt5;
- PySide6.

---

Para rodar esse código também são necessários os sequintes arquivos:
- EIT_2D_main.py programa principal;
- EIT_functions_2D.py contém as funções necessárias;
- EIT_mesh_2D.py com as seguintes malhas 1D:
  - DoisTriangulos_4elem.msh 4 elementos 5 nós;
  - octogono_2208_mod.msh 146 elementos 90 nós ou.

    
Nessa etapa, o programa resolve somente o "Problema Direto". 
O programa plota uma MEFinitos como as condutividade dos elementos e
imprime as tensões calculadas nos eletrodos para cada padrão de corrente

---

Para escolher uma malha, descomente a seguinte linha de comando no arquivo EIT_1D.py:

` ` `
#opcao = input('Escolha a malha(1, 2): ')
` ` ` 

ou altere a seguinte linha de comando como desejado.

` ` `
opcao = '1' # as opcoes são 1 e 2
` ` ` 



