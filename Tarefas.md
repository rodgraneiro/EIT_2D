Trarefas
#04set25
Sugestão de ata para reunião PhD Edson 04/09/25.

- Discutir dúvidas sobre Classes em python;

- Discutir dúvida: Como implementar uma malha Gmsh com Hua?

- Apresentar código para criar malhas adaptado para rodar com classes e eletrodo Hua;

- Temas gerais. 


#21ago25
(21/08) Sugestão de tarefas para próximas semanas:
 - deixar todas as malhas (e os respectivos .geo) em uma pasta 'malhas' no repositório, e adaptar os códigos para usar essas malhas;
 - implementar physical groups nas malhas 2D, de modo que os valores de resistividade ou condutividade sejam passados por physical group, e não por número do elemento, assim os valores do problema direto podem ser passados independente do número de elementos da malha;
 - verificar qual o maior número possível (aproximado) de elementos que você pode rodar no seu computador;
 - verificar como implementar injeção de corrente em uma superfície, e não em um ponto (modelo de eletrodo HUA), e começar a implementar no 2D;
 - verificar como usar a API do gmsh para python, para criar as malhas automaticamente e sem abrir a interface gráfica do gmsh:
     - tutoriais python em "https://gitlab.onelab.info/gmsh/gmsh/-/tree/master/tutorials/python";
     - para gravar o msh no formato 2.2, antes de gravar incluir o comando "gmsh.option.setNumber("Mesh.MshFileVersion",2.2)";

# 24jul25
## Tarefa1: assisttir aulas github: 
Curso de Git, Github usando o GitDesktop: https://www.youtube.com/watch?v=xEKo29OWILE&list=PLHz_AreHm4dm7ZULPAmadvNhH6vk9oNZA
## Tarefa2: ler artigos:

# 17ago25
## Tarefa: ler artigo:
Oi Edson, preparando a aula de hj me deparei com este artigo de revisão:
https://www.sciencedirect.com/science/article/pii/S1388245720305678
Ele fala alguma coisa sobre a medição da anisotropia do músculo dentro do contexto da Miografia de Impedância Elétrica. Acho bom dar uma olhada.

# 14ago25
## ler artigos:
2007_Felipe_Abascal_Validation_of_a_finite-element_solution_for_EIT_in_an_anisotropic_medium
2017_Gonzalez_Isotropic and anisotropic total variation regularization in electrical impedance tomography
2015_Zhou_Comparison of total variation algorithms for electrical impedance tomography
EIT_in_anisotropic_media_with_known_eigenvector


# 10ago25
## Tarefa1: estudar tutoriais do PyQtgraph:
link: https://www.pyqtgraph.org/
## Tarefa2: estudar tutoriais pyEIT
link: https://github.com/eitcom/pyEIT
## Tarefa3: ler artigo:
2012_Evaluating_EIT_image_resolution_using_anatomical_atlas_TMSi2012_Camargo