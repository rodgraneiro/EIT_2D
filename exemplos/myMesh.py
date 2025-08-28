import numpy as np
import meshio
import matplotlib.pyplot as plt
import matplotlib.tri as tri


class myMesh():
    def __init__(self, msh_file_name, show_info=False, physical_tags=None):
        self.__mshdata = meshio.read(msh_file_name)
        self.Ygeo = None
        self.Y = None

        if physical_tags is None:
            self.tags = {}
            self.tags['gnd'] = 10000          # nó GND
            self.tags['electrodes'] = 10001   # primeiro eletrodo
            self.tags['elements'] = 1000      # corpo... objetos são 1001, 1002 etc
        else:
            self.tags = physical_tags
            # ToDo: verificar se tags necessárias existem


        # Check msh dimension
        if 'gmsh:physical' in self.__mshdata.cell_data_dict.keys():
            if 'tetra' in self.__mshdata.cell_data_dict["gmsh:physical"].keys():
                self.dim = 3 # 3D mesh
                self.element_type = 'tetra'
            elif 'triangle' in self.__mshdata.cell_data_dict["gmsh:physical"].keys():
                self.dim = 2 # 2D mesh
                self.element_type = 'triangle'
            else:
                raise Exception("myMesh(): dimension identification error.")
        else:
            raise Exception("myMesh(): invalid mesh (no physical entities found).")
        
        # Getting coordinates and topology
        self.coordinates = self.__mshdata.points       # Monta matriz coordenadas, sempre 3 dimensões (mesmo 2D)
        self.topology = self.__mshdata.cells_dict[self.element_type]      # Monta matriz topologia
        self.physical_groups = self.__mshdata.cell_data_dict["gmsh:physical"][self.element_type]
        self.physical_tags = np.unique(self.physical_groups)

        self.number_of_nodes = self.coordinates.shape[0]
        self.number_of_elements = self.topology.shape[0]

        self.sigma_vec = np.ones(self.number_of_elements)*(-1)   # deixa inicialmente como -1 por segurança
        self.sigma_defined = False


        # Getting electrodes data
        self.hua_electrode_available = False
        if (self.dim == 3) and ('triangle' in self.__mshdata.cell_data_dict["gmsh:physical"].keys()):
            self.hua_electrode_available = True   # ToDo: implementar leitura das faces
        elif (self.dim == 2) and ('line' in self.__mshdata.cell_data_dict["gmsh:physical"].keys()):
            self.hua_electrode_available = True   # ToDo: implementar leitura das linhas

        node_physical_groups = self.__mshdata.cell_data_dict["gmsh:physical"]["vertex"]
        self.gnd_node = -1 # '-1' para identificar quando não encontra no msh
        if self.tags['gnd'] in node_physical_groups:
            self.gnd_node = np.squeeze(self.__mshdata.cells_dict["vertex"][node_physical_groups==self.tags['gnd']])
        self.electrodes_center_node = np.squeeze(self.__mshdata.cells_dict["vertex"][node_physical_groups>=self.tags['electrodes']])
        self.number_of_electrodes = len(self.electrodes_center_node)

        if show_info:
            print(f'Mesh {msh_file_name} ({self.dim}D) with {self.number_of_nodes} nodes and {self.number_of_elements} elements.')
            print(f'{self.number_of_electrodes} Electrodes found: nodes {self.electrodes_center_node}')
            if self.gnd_node != -1:
                print(f'GND node found: {self.gnd_node}')
            else:
                print('GND node not found.')
    

    # sigma pode ser um vetor de tamanho [self.number_of_elements]
    # ou pode ser um dicionário com itens: sigma['tag'] = valor
    def apply_sigma(self, sigma):
        if len(sigma) == self.number_of_elements:
            self.sigma_defined = True
            self.sigma_vec = sigma
        elif len(sigma) == len(self.physical_tags):
            # verifica se todas as tags foram fornecidas
            for tag in self.physical_tags:
                if tag not in sigma.keys():
                    raise Exception(f"myMesh.apply_sigma(): Tag {tag} not found!")
            self.sigma_defined = True
            for tag in sigma.keys():
                self.sigma_vec[self.physical_groups==tag] = sigma[tag]

    
    # Plota Grafico de uma condutividade na malha
    def plot(self, figsize=(6, 5), title=''):
        if self.dim==2:
            if not(self.sigma_defined):
                sigma_vec = np.zeros(self.number_of_elements)
            else:
                sigma_vec = self.sigma_vec

            if title == '':
                title = r"Condutividade ($\sigma$)"
            
            x, y = self.coordinates[:, 0], self.coordinates[:, 1]                                                 
            triang = tri.Triangulation(x, y, self.topology)
            fig, ax = plt.subplots(figsize=figsize)
            tpc = ax.tripcolor(triang, facecolors=sigma_vec, edgecolors='k',cmap='Blues' )
            fig.colorbar(tpc, ax=ax, label=r'$\sigma$ [S/m]')
            ax.set_title(title)
            plt.xlabel("[m]")
            plt.ylabel("[m]")
            plt.tight_layout()
            plt.axis('equal')
            plt.show()
        else: 
            print("Plot for 3D meshes not available.") # ToDo: implementar visualização para malhas 3D


    ###############################################################################
    # Esta função calcula a matriz local de condutividade da malha de elementos
    # finitos
    #
    #            sigma_i                                      
    # Y_local = ------------- * [matriz 2d ***completar***]             
    #             4*A_i                                             
    ###############################################################################
    # OBS: só parte geométrica (que nunca muda), sem multiplicar por sigma (que muda)
    def calc_Ygeo_local_triangle(self, elem):
        noh1 = int(self.topology[elem][0])
        noh2 = int(self.topology[elem][1])
        noh3 = int(self.topology[elem][2])
        
        x = [self.coordinates[noh1][0], self.coordinates[noh2][0], self.coordinates[noh3][0]]
        y = [self.coordinates[noh1][1], self.coordinates[noh2][1], self.coordinates[noh3][1]]

        triangulo = np.array([ [1, x[0],  y[0]], [1, x[1],  y[1]], [1, x[2],  y[2]]], dtype=np.float64)
        area_triangulo = abs((np.linalg.det(triangulo))/2)
    
        C_11 =  ( y[1]-y[2])**2 + (x[2]-x[1])**2
        C_12 =  ( y[1]-y[2])*(y[2]-y[0])+(x[2]-x[1])*(x[0]-x[2])          #A_21 = A_12
        C_13 =  ( y[1]-y[2])*(y[0]-y[1])+(x[2]-x[1])*(x[1]-x[0])          #A_31 = A_13
        C_22 =  (y[2]-y[0])**2 + (x[0]-x[2])**2                           #A_32 = A_23
        C_23 =  (y[2]-y[0])*(y[1]- y[1])+(x[0]-x[2])*(x[1]-x[0])          # Na tese está: (Bt_m*Bt_o)+(Gm_m*Gm_o)+(Dt_m*Dt_o)
        C_33 =  (y[0]- y[1])**2 + (x[1]-x[0])**2  

        
        matriz_local = (1.0/(4.0*area_triangulo))*np.array([[C_11, C_12, C_13], 
                                        [C_12, C_22, C_23],      
                                        [C_13, C_23, C_33]       
                                        ])

        return matriz_local


    ###############################################################################
    # Calcula Ygeo local para todos os elementos
    ###############################################################################
    def calc_Ygeo(self):
        if self.Ygeo is None: # Se já calculou antes, não precisa recalcular...
            if self.dim == 2:
                self.Ygeo = np.zeros((self.number_of_elements,3,3))
                for elem in range(self.number_of_elements):
                    self.Ygeo[elem] = self.calc_Ygeo_local_triangle(elem)
            else:
                raise Exception(f"myMesh.calc_Ygeo(): Ainda não implementado para 3D!") # ToDo: implementar 3D


    ###############################################################################
    # Esta função calcula a matriz global de condutividade da malha de elementos
    # finitos
    ###############################################################################
    def calc_Y_global(self):
        if self.sigma_defined:
            self.calc_Ygeo()

            self.Y = np.zeros((self.number_of_nodes, self.number_of_nodes))

            for e, conn in enumerate(self.topology):
                for i in range(3):
                    for j in range(3):
                        self.Y[conn[i], conn[j]] += self.sigma_vec[e]*self.Ygeo[e][i, j]
        else:
            raise Exception(f"myMesh.calc_Y_global(): Sigma not defined!")
    