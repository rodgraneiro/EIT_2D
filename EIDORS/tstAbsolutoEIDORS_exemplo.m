%===========================================================
% Exemplo mínimo (EIDORS) – reconstrução ABSOLUTA de σ
% (funciona em Matlab e normalmente em Octave)
%===========================================================
% Objetivo:
%   - gerar dados absolutos (vi) a partir de um "fantoma" com σ absoluto
%   - reconstruir σ absoluto com Gauss-Newton (inv_solve), usando CEM e prior
%
% Observação importante:
%   - isso é ABSOLUTO (usa vi), NÃO diferencial (vh-vi).
%===========================================================

clear all; close all; clc;
% Ajuste o caminho do seu EIDORS se necessário
run('D:\EIDORS\eidors-v3.10\eidors\startup.m');
addpath('D:\EIDORS\eidors-v3.10\eidors\examples');
eidors_cache clear;
graphics_toolkit fltk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARTE 1: GERAÇÃO DE MEDIDAS SIMULADAS (solução do problema direto)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------------
% 1) Modelo FEM + eletrodos
%-----------------------------
% c2d0c: circular 2D (distmesh), 16 eletrodos
imdl = mk_common_model('c2d4c', 16);
fmdl = imdl.fwd_model;

% Para reconstrução ABSOLUTA, é MUITO recomendável usar CEM:
% (impedância de contato > 0 em todos os eletrodos)
zc = 0.01; % ohm*m^2 (valor "chute" para exemplo; calibração importa muito!)
for i = 1:numel(fmdl.electrode)
    fmdl.electrode(i).z_contact = zc;
end

%-----------------------------
% 2) Definir estímulos / medições
%-----------------------------
% Padrão típico: adjacente (corrente entre eletrodos vizinhos)
% e medir tensões em todos os pares possíveis (exceto eletrodos de corrente).
% (Se sua versão do EIDORS tiver padrões diferentes, este é o mais comum)
n_elec = numel(fmdl.electrode);
stim = mk_stim_patterns(n_elec, 1, [0 1], [0 1], {}, 1);
fmdl.stimulation = stim;

%-----------------------------
% 3) Gerar "dados medidos" absolutos (vi)
%-----------------------------
% Fundo absoluto (atribuição de distribuição uniforme de condutividade sigma0 à malha fmdl)
sigma0 = 3.0;  % S/m (exemplo)

img0 = mk_image(fmdl, sigma0);

% Inserir uma anomalia (σ mais alta) – simples e "mínimo"
% A função mk_c2f_circ_mapping existe em algumas configs; pra ser robusto,
% vamos fazer uma anomalia por "região" usando o centro dos elementos:
%ctr = interp_mesh(fmdl);              % centroids/estruturas auxiliares
%cc  = ctr.elem_centre;                % [ne x 2] centros dos elementos
%r   = sqrt(cc(:,1).^2 + cc(:,2).^2);
%theta = atan2(cc(:,2), cc(:,1));

% --- Centros dos elementos (centróides) sem interp_mesh ---
E = fmdl.elems;   % [ne x 3] nós de cada triângulo
P = fmdl.nodes;   % [nn x 2] coordenadas (x,y) dos nós

cc = (P(E(:,1),:) + P(E(:,2),:) + P(E(:,3),:)) / 3;  % [ne x 2] centróides





sigma_true = sigma0 * ones(size(img0.elem_data));
% anomalia: um "blob" em (x≈0.3, y≈0.0) com raio ~0.15
x0 = 0.30; y0 = 0.00; rad = 0.15;
d = sqrt( (cc(:,1)-x0).^2 + (cc(:,2)-y0).^2 );
sigma_true(d < rad) = 2.0;  % S/m (anomalia mais condutiva)

img_true = img0;
img_true.elem_data = sigma_true;

% Simular tensões ABSOLUTAS
vi = fwd_solve(img_true);

% (Opcional) adicionar ruído pequeno para parecer mais real
noise_level = 0.001; % 0.1% (exemplo)
#vi.meas = vi.meas .* (1 + noise_level*randn(size(vi.meas)));

% Simular tensões ABSOLUTAS
#vi = fwd_solve(img_true);

% (Opcional) adicionar ruído pequeno para parecer mais real
#noise_level = 0.001; % 0.1% (exemplo)

if isstruct(vi) && isfield(vi,'meas')
    vi.meas = vi.meas .* (1 + noise_level*randn(size(vi.meas)));
else
    % vi veio como vetor numérico: adiciona ruído e empacota como "data"
    vi_num = vi;
    vi_num = vi_num .* (1 + noise_level*randn(size(vi_num)));

    vi = eidors_obj('data','vi');
    vi.meas = vi_num;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARTE 2: GERAÇÃO DA IMAGEM DE EIT (solução do problema inverso)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------------
% 4) Configurar inversão ABSOLUTA (Gauss-Newton)
%-----------------------------
imdl_base = mk_common_model('c2d0c', 16);
fmdl_base = imdl_base.fwd_model;

#imdl_base = eidors_obj('inv_model', 'ABS_GN_minimo');
imdl_base.fwd_model = fmdl_base;

% Solver Gauss-Newton (absoluto)
imdl_base.solve = @inv_solve_gn;            % iterativo GN
imdl_base.reconst_type = 'absolute';

% Jacobiano padrão
imdl_base.jacobian = @jacobian_adjoint;
imdl_base.RtR_prior = @prior_laplace;   %@prior_tikhonov;      % prior mínimo (0ª ordem)
imdl_base.hyperparameter.value = 1e-5;      % ajuste conforme necessário

% Chute inicial (bkgnd) – SUPER importante no absoluto
imdl_base.jacobian_bkgnd.value = sigma0;

% Controle de iterações (para ficar “mínimo” e estável)
imdl_base.inv_solve_gn.max_iterations = 20;
imdl_base.inv_solve_gn.tol = 1e-10;
imdl_base.inv_solve_gn.verbose = 2;

%-----------------------------
% 5) Reconstrução ABSOLUTA
%-----------------------------
img_rec = inv_solve(imdl_base, vi);

%-----------------------------
% 6) Visualização
%-----------------------------
figure;
subplot(1,2,1);
show_fem(img_true, 1);
title('Verdadeiro \sigma (absoluto)');

subplot(1,2,2);
show_fem(img_rec, 1);
title('Reconstruído \sigma (absoluto)');

% Comparação numérica simples (média do fundo e máximo)
fprintf('sigma0 = %.3f S/m\n', sigma0);
fprintf('Verdadeiro:  mean=%.3f  max=%.3f\n', mean(img_true.elem_data), max(img_true.elem_data));
fprintf('Reconstr.:  mean=%.3f  max=%.3f\n', mean(img_rec.elem_data),  max(img_rec.elem_data));

%===========================================================
% DICAS (curtas) se ficar “ruim”:
%  - ajuste z_contact (zc) e hyperparameter
%  - aumente max_iterations
%  - troque prior_tikhonov por prior_laplace (suaviza melhor)
%  - garanta que seu padrão de estimulação é o mesmo do seu hardware/dados
%===========================================================

