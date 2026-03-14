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


imdl = mk_common_model('c2d2c', 4);
fmdl = imdl.fwd_model;


%n_elec = length(fmdl.electrode);
%I = 1e-3;   % 1 mA

%for k = 1:n_elec
%    stim = zeros(n_elec,1);

%    stim(k) =  I;
%    stim(mod(k-1 + n_elec/2, n_elec) + 1) = -I;

%    fmdl.stimulation(k).stim_pattern = stim;
%end
%for k = 1:length(fmdl.stimulation)
%    fprintf('\nEstimulo %d\n',k);
%    full(fmdl.stimulation(k).stim_pattern)
%end


I = 1e-3;   % 1 mA

% estabelece corrente de 1mA
for k = 1:length(fmdl.stimulation)
    fmdl.stimulation(k).stim_pattern = ...
        I * fmdl.stimulation(k).stim_pattern / ...
        max(abs(fmdl.stimulation(k).stim_pattern));
end

nodes = fmdl.nodes;      % [N x 2] ou [N x 3]
elems = fmdl.elems;

% cria uma imagem FEM com condutividade uniforme igual a 1 em todos os elementos da malha.
img = mk_image(fmdl, 1);

n_elem = size(fmdl.elems,1);

% condutividades anisotrópicas
sigma_xx = 1.0 * ones(n_elem,1);
sigma_yy = 0.10 * ones(n_elem,1);
sigma_xy = 0.0 * ones(n_elem,1);

%define o tensor de condutividade anisotrópica em cada elemento da malha.
img.elem_data = [sigma_xx, sigma_xy, sigma_yy];

% cálculo do problema direto
vh = fwd_solve(img);

% cria figura anisotrópica sxx
figure;
img_xx = img;
img_xx.elem_data = sigma_xx;
show_fem(img_xx);
colorbar;
title('\sigma_{xx}');

fprintf('tensões  = %.8f V\n', vh.meas);

% cria figura anisotrópica syy
figure;
img_yy = img;
img_yy.elem_data = sigma_yy;
show_fem(img_yy);
colorbar;
title('\sigma_{yy}');


% cria figura de anisotropia
figure;
anisotropia = img;
anisotropia.elem_data = abs(sigma_xx - sigma_yy);
show_fem(anisotropia);
colorbar;
title('anisotropia');

#fprintf('anisotropia  = %.3f V\n', anisotropia.elem_data);

%salva a matriz nodes/elems em um arquivo de texto chamado nodes.txt
#dlmwrite('nodes.txt', nodes, 'delimiter',' ');
#dlmwrite('elems.txt', elems, 'delimiter',' ');


vh.meas

% imprimir todos os estímulos de corrente e medição de tensões
for i = 1:length(fmdl.stimulation)
    fprintf('\nEstimulo %d\n',i);
    disp(fmdl.stimulation(i))
end

