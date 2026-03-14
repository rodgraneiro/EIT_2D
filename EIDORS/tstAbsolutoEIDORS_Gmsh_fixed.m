%===========================================================
% tstAbsolutoEIDORS_Gmsh_fixed.m
% EIDORS v3.10 (Octave/Matlab) – Reconstrucao ABSOLUTA de sigma
%===========================================================
% Objetivo:
%   - Gerar dados absolutos (vh) com um modelo DIRETO (malha Gmsh mais densa)
%   - Reconstruir sigma absoluto com um modelo INVERSO (malha Gmsh mais grossa)
%   - Garantir que o fwd_model do INVERSO tenha os campos exigidos pelo EIDORS:
%       * solve, system_mat, jacobian, stimulation, electrodes, etc.
%
% Principais correcoes vs tstAbsolutoEIDORS_Gmsh.m:
%   (1) O fwd_model do INVERSO (inv_fwd_model) agora define:
%         inv_fwd_model.solve      = @fwd_solve_1st_order;
%         inv_fwd_model.system_mat = @system_mat_1st_order;
%         inv_fwd_model.jacobian   = @jacobian_adjoint;
%       Isso evita: "error: structure has no member 'solve'" vindo de fwd_solve().
%   (2) O inv_model (imdl_inv) e criado como eidors_obj('inv_model',...)
%       e NAO como um struct "cru".
%   (3) Variaveis img_true e sigma0 foram normalizadas (sigma_base).
%===========================================================

clear all; close all; clc;

run('D:\EIDORS\eidors-v3.10\eidors\startup.m');
addpath('D:\EIDORS\eidors-v3.10\eidors\examples');
%-----------------------------------------------------------
% 0) Inicializacao do EIDORS
%-----------------------------------------------------------
% (Opcional) Escolha de toolkit no Octave (se existir)
try
    graphics_toolkit fltk
catch
end

% Garanta que o EIDORS esteja no path:
% run('D:\EIDORS\eidors-v3.10\eidors\startup.m'); % ajuste se necessario

%-----------------------------------------------------------
% 1) Carregar malhas Gmsh (direta e inversa)
%-----------------------------------------------------------
% Assumindo que voce ja tem rotinas/arquivos que carregam a malha
% e retornam uma struct S com campos:
%   S.pts : [N x 2] coordenadas
%   S.tri : [Ne x 3] conectividade (1-based)
%   S.phys: tags (opcional)
%
% >>> Ajuste estes loads para o seu caso real <<<
%
% Exemplo (se voce salvou em .mat):
% load('gmsh_direct.mat','S');       % malha densa (direto)
% load('gmsh_inverse.mat','S_base'); % malha grossa (inverso)

% ---- NO SEU ARQUIVO ORIGINAL, parece que voce ja tem S e S_base prontos.
% Se for o caso, comente os loads acima e mantenha a forma como voce cria S/S_base.

% Se o seu script anterior ja define S e S_base, deixe como esta.
% Caso contrario, descomente e ajuste os loads.


%S = load('dezesseis_triangulos_22jan25.mat');
S = load('circ16_3_anomalia6_v1_2208.mat');
%S = load('circ16_3_anomalia6_v3.mat');
%S = load('circ16_3_anomalia6_v3.mat');


%S_base = load('quatro_base_22jan25.mat');
S_base = load('circ16_base_v2208.mat');
%S_base = load('circ16_base_coarse.mat');

phys = S.phys;                     % Ne x 1

%-----------------------------------------------------------
% 2) Construir fwd_model do DIRETO (malha densa)
%-----------------------------------------------------------
fwd_model = eidors_obj('fwd_model','gmsh_direct');
fwd_model.nodes = S.pts;      % Nx2
fwd_model.elems = S.tri;      % Ne x 3
fwd_model.gnd_node = 17;      % AJUSTE para o seu modelo (ou detecte automaticamente)

% Eletrodos pontuais (1 no por eletrodo)
n_elec = 16;

for i = 1:n_elec
    fwd_model.electrode(i).nodes = i;      % <-- AJUSTE: mapeie para os nos reais dos eletrodos
    fwd_model.electrode(i).z_contact = 0;  % pontual / sem CEM
end

% Stimulation / medidas
stim = mk_stim_patterns(n_elec, 1, [0 1], [0 1], {}, 1);
fwd_model.stimulation = stim;

% Campos essenciais para o EIDORS (importante!)
fwd_model.solve      = @fwd_solve_1st_order;
fwd_model.system_mat = @system_mat_1st_order;
fwd_model.jacobian   = @jacobian_adjoint;

%-----------------------------------------------------------
% 3) Gerar dados ABSOLUTOS (vh) do modelo direto
%-----------------------------------------------------------
sigma_true = 3.0;                 % "verdadeiro" (exemplo)
img_true   = mk_image(fwd_model, sigma_true);
img_true.elem_data(phys == 1000) = 3.0;  % body
img_true.elem_data(phys == 1001) = 2.0;  % anomalia1 (exemplo)
img_true.elem_data(phys == 1002) = 2.0;  % anomalia2
img_true.elem_data(phys == 1003) = 2.0;  % anomalia3

vh = fwd_solve(img_true);         % dados absolutos (tensoes)

%-----------------------------------------------------------
% 4) Construir fwd_model do INVERSO (malha grossa)
%-----------------------------------------------------------
inv_fwd_model = eidors_obj('fwd_model','gmsh_inverse');
inv_fwd_model.nodes = S_base.pts;
inv_fwd_model.elems = S_base.tri;
inv_fwd_model.gnd_node = 17;      % AJUSTE coerente com o inverso

for i = 1:n_elec
    inv_fwd_model.electrode(i).nodes = i;     % <-- AJUSTE: mapeie para os nos reais dos eletrodos
    inv_fwd_model.electrode(i).z_contact = 0;
end

inv_fwd_model.stimulation = stim;






% >>>>> ESTA ERA A CAUSA DO ERRO "no member 'solve'" <<<<<
inv_fwd_model.solve      = @fwd_solve_1st_order;
inv_fwd_model.system_mat = @system_mat_1st_order;
%inv_fwd_model.jacobian   = @jacobian_adjoint;

%-----------------------------------------------------------
% 5) Montar inv_model (Gauss-Newton) de forma "EIDORS-like"
%-----------------------------------------------------------
imdl_inv = eidors_obj('inv_model','inv_gmsh_abs');
imdl_inv.fwd_model = inv_fwd_model;

imdl_inv.solve        = @inv_solve_gn;
imdl_inv.reconst_type = 'absolute';

% Jacobiano (escolha 1)
%imdl_inv.jacobian     = @jacobian_adjoint;

% Regularizacao
imdl_inv.RtR_prior    = @prior_laplace;
%imdl_inv.RtR_prior    =@prior_tikhonov;
imdl_inv.hyperparameter.value = 1e-5;

% Fundo do jacobiano (importantissimo em absoluto)
sigma_base = 1.0;
imdl_inv.jacobian_bkgnd.value = sigma_base;









% Opcoes do GN (dependem do seu build do EIDORS; se nao existirem, serao ignoradas)
imdl_inv.inv_solve_gn.jacobian        = @jacobian_adjoint; % garante compatibilidade com parse_options
imdl_inv.inv_solve_gn.max_iterations  = 1;
imdl_inv.inv_solve_gn.min_value       = 0.5;   % evita sigma negativa
imdl_inv.inv_solve_gn.verbose         = 2;


%-----------------------------------------------------------
% 6) Reconstrucao ABSOLUTA
%-----------------------------------------------------------

img_rec = inv_solve(imdl_inv, vh);

%-----------------------------------------------------------
% 7) Visualizacao
%-----------------------------------------------------------
figure;
subplot(1,2,1);
show_fem(img_true, 1);
title('Verdadeiro \sigma (absoluto)');

subplot(1,2,2);
show_fem(img_rec, 1);
title('Reconstruido \sigma (absoluto)');

fprintf('sigma_true = %.3f S/m\n', sigma_true);
fprintf('sigma_base = %.3f S/m\n', sigma_base);


img0 = mk_image(imdl_inv.fwd_model, sigma_base);
img0.fwd_model.jacobian = @jacobian_adjoint;
J = calc_jacobian(img0);


disp('Jacobiano J =');
disp(J);

%fprintf('size(J) = %d x %d\n', size(J,1), size(J,2));
fprintf('rank(J) = %d\n', rank(J));

s = svd(J);

disp('Valores singulares de J =');
disp(s);

figure;
semilogy(s,'o-');
grid on;
xlabel('indice');
ylabel('valor singular');
title('Espectro singular do Jacobiano (adjoint)');

