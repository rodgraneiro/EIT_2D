clear all; close all; clc;

run('D:\EIDORS\eidors-v3.10\eidors\startup.m');
addpath('D:\EIDORS\eidors-v3.10\eidors\examples');

%-----------------------------------------------------------
% 0) Inicializacao do EIDORS
%-----------------------------------------------------------
try
    graphics_toolkit fltk
catch
end

#imdl = mk_common_model('c2d2c', 4);
imdl = mk_common_model('a2c0', 4);
fmdl = imdl.fwd_model;

nodes    = fmdl.nodes;
elems    = fmdl.elems;
boundary = find_boundary(elems);

n_nodes = size(nodes,1);
n_elems = size(elems,1);
n_bnd   = size(boundary,1);
n_elec  = length(fmdl.electrode);

% ----------------------------------------------------------
% Encontrar o nó mais próximo do centro geométrico
% ----------------------------------------------------------
centro = mean(nodes, 1);
dist2 = (nodes(:,1) - centro(1)).^2 + (nodes(:,2) - centro(2)).^2;
[~, idx_gnd] = min(dist2);

disp(['Nó mais próximo do centro: ', num2str(idx_gnd)]);
disp('Coordenadas do nó:');
disp(nodes(idx_gnd,:));

% -----------------------------------
% classificar arestas de contorno
% edge_tag = 0 -> contorno sem eletrodo
% edge_tag = k -> eletrodo k
% -----------------------------------
edge_tag = zeros(n_bnd,1);

for e = 1:n_bnd
    n1 = boundary(e,1);
    n2 = boundary(e,2);

    for k = 1:n_elec
        enodes = fmdl.electrode(k).nodes(:);
        if ismember(n1,enodes) && ismember(n2,enodes)
            edge_tag(e) = k;
            break;
        end
    end
end

fid = fopen('a2c0_EIDORS.msh','w');

% =========================
% MeshFormat
% =========================
fprintf(fid, '$MeshFormat\n');
fprintf(fid, '2.2 0 8\n');
fprintf(fid, '$EndMeshFormat\n');

% =========================
% PhysicalNames
% agora inclui:
% 1 ponto gnd
% 1 contorno
% n_elec eletrodos
% 1 domínio
% total = n_elec + 3
% =========================
fprintf(fid, '$PhysicalNames\n');
fprintf(fid, '%d\n', n_elec + 3);

% dimensão 0 = ponto
fprintf(fid, '0 %d "gnd"\n', 10000);

% dimensão 1 = linhas
fprintf(fid, '1 %d "contorno"\n', 1);

for k = 1:n_elec
    fprintf(fid, '1 %d "eletrodo_%d"\n', 5000+k, k);
end

% dimensão 2 = superfície
fprintf(fid, '2 %d "dominio"\n', 1000);
fprintf(fid, '$EndPhysicalNames\n');

% =========================
% Nodes
% =========================
fprintf(fid, '$Nodes\n');
fprintf(fid, '%d\n', n_nodes);

for i = 1:n_nodes
    fprintf(fid, '%d %.16g %.16g %.16g\n', i, nodes(i,1), nodes(i,2), 0.0);
end

fprintf(fid, '$EndNodes\n');

% =========================
% Elements
% agora inclui:
% n_bnd segmentos
% n_elems triângulos
% 1 ponto gnd
% total = n_bnd + n_elems + 1
% =========================
fprintf(fid, '$Elements\n');
fprintf(fid, '%d\n', n_bnd + n_elems + 1);

elem_id = 1;

% ---- segmentos de contorno (type 1)
for e = 1:n_bnd
    n1 = boundary(e,1);
    n2 = boundary(e,2);

    if edge_tag(e) == 0
        continue;   %  pula contorno comum
    end

    phys_tag = 5000 + edge_tag(e);   % eletrodo_k
    geom_tag = phys_tag;

    fprintf(fid, '%d 1 2 %d %d %d %d\n', elem_id, phys_tag, geom_tag, n1, n2);
    elem_id = elem_id + 1;
end

% ---- triângulos (type 2)
for i = 1:n_elems
    n1 = elems(i,1);
    n2 = elems(i,2);
    n3 = elems(i,3);

    phys_tag = 1000;
    geom_tag = 1000;

    fprintf(fid, '%d 2 2 %d %d %d %d %d\n', elem_id, phys_tag, geom_tag, n1, n2, n3);
    elem_id = elem_id + 1;
end

% ---- ponto gnd (type 15)
% associa o physical_tag = 10000 ao nó idx_gnd
fprintf(fid, '%d 15 2 %d %d %d\n', elem_id, 10000, 10000, idx_gnd);
elem_id = elem_id + 1;

fprintf(fid, '$EndElements\n');

fclose(fid);

disp('Arquivo .msh exportado com sucesso.');
disp(['Nó GND exportado com physical_tag = 10000: nó ', num2str(idx_gnd)]);
