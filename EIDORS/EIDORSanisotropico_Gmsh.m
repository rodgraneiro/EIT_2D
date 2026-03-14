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


%S = load('circ16_3_anomalia6_v1_2208.mat');
S = load('fantomaAnisotropia01.mat');
%S = load('quatro_4e_em_XY.mat');
%
%S = load('circ4_anomZero_Hua.mat');
%S = load('circ16_base_v2208.mat');

phys = S.phys;                     % Ne x 1

%-----------------------------------------------------------
% 2) Construir fwd_model do DIRETO (malha densa)
%-----------------------------------------------------------
fwd_model = eidors_obj('fwd_model','gmsh_direct');
fwd_model.nodes = S.pts;      % Nx2
fwd_model.elems = S.tri;      % Ne x 3
fwd_model.gnd_node = 1;      % AJUSTE para o seu modelo (ou detecte automaticamente)

% Eletrodos pontuais (1 no por eletrodo)
n_elec = 4;

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
% (opcional) para remover o warning do normalize_measurements
%fwd_model.normalize_measurements = 0;


%-----------------------------------------------------------
% 3) Gerar dados ABSOLUTOS (vh) do modelo direto
%-----------------------------------------------------------
%sigma_true = 0.1;                 % "verdadeiro" (exemplo)
%img_true   = mk_image(fwd_model, sigma_true);

%idxb = (phys == 6);

%img_true.elem_data(idxb, 1) = 3.0;  % sigma_xx
%img_true.elem_data(idxb, 2) = 0.0;  % sigma_xy
%img_true.elem_data(idxb, 3) = 3.0;  % sigma_yy

%idx = (phys == 7);

%img_true.elem_data(idx, 1) = 2.0;  % sigma_xx
%img_true.elem_data(idx, 2) = 0.0;  % sigma_xy
%img_true.elem_data(idx, 3) = 2.0;  % sigma_yy

n_elem = size(fwd_model.elems,1);

img_true = mk_image(fwd_model, 1);        % valor placeholder
img_true.elem_data = zeros(n_elem,3);     % força tensor Ne×3

% fundo (por exemplo sigma=3 isotrópico)
img_true.elem_data(:,1) = 1.0;   % sigma_xx
img_true.elem_data(:,2) = 0.0;   % sigma_xy
img_true.elem_data(:,3) = 0.1;   % sigma_yy

% quadrado: sigma=2 isotrópico
%idx = (phys == 7);
%img_true.elem_data(idx,1) = 1.0;
%img_true.elem_data(idx,2) = 0.0;
%img_true.elem_data(idx,3) = 1.0;



vh = fwd_solve(img_true);
%fprintf('tensões  = %.3f V\n', vh.meas);
%figure;
%img_xx = img_true;
%img_xx.elem_data = img_true.elem_data(:,1);

imgA = mk_image(img_true, 1.0);

S = calc_system_mat(imgA);

K = full(S.E);

#disp(K(1:10,1:10))  % imprime parte da matriz

#disp(K);


elem = 1;  % número do elemento

nodes = imgA.fwd_model.nodes;
elems = imgA.fwd_model.elems;

idx = elems(elem,:);
coords = nodes(idx,:);

area = polyarea(coords(:,1), coords(:,2));

%disp(area)

figure;
img_xx = img_true; img_xx.elem_data = img_true.elem_data(:,1);
show_fem(img_xx); caxis([0 3]); colorbar; title('\sigma_{xx}');

%v = img_xx.elem_data;  % ou sigma_xx
%fprintf('min=%g max=%g  nNaN=%d  nElem=%d\n', min(v), max(v), sum(isnan(v)), numel(v));
%show_fem(img_xx);
%ax = gca;
%set(gca, 'clim', [0 3]);   % em Octave geralmente funciona
%caxis([0 3]);   % evita enganar o olho
%colorbar;
%title('\sigma_{xx}');

figure;
img_yy = img_true;
img_yy.elem_data = img_true.elem_data(:,3);
show_fem(img_yy);
%ax = gca;
%ax.CLim = [0 3];
set(gca, 'clim', [0 3]);   % em Octave geralmente funciona
%caxis([0 3]);   % evita enganar o olho
colorbar;
title('\sigma_{yy}');

figure;
img_A = img_true;
img_A.elem_data = abs(img_true.elem_data(:,1) - img_true.elem_data(:,3));
show_fem(img_A);
%ax = gca;
%ax.CLim = [0 3];
set(gca, 'clim', [0 3]);   % em Octave geralmente funciona
%caxis([0 3]);   % evita enganar o olho
colorbar;
title('|σ_{xx} - σ_{yy}|');


fprintf('tensões  = %.8f V\n', vh.meas);


% imprimir todos os estímulos de corrente e medição de tensões
for i = 1:length(fwd_model.stimulation)
    fprintf('\nEstimulo %d\n',i);
    disp(fwd_model.stimulation(i))
end







