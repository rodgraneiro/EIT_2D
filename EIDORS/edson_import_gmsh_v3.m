%graphics_toolkit fltk
clc
S = load('circ16_3_anomalia6_v3.mat');
whos('-file','circ16_3_anomalia6_v3.mat')

% conferir campos principais
size(S.pts)
size(S.tri)
size(S.phys)


fwd_model = eidors_obj('fwd_model','gmsh_import');
fwd_model.nodes = S.pts;      % Nx2
fwd_model.elems = S.tri;      % Ne x 3 (já em indexação MATLAB 1-based)
fwd_model.gnd_node = 1;       % provisório (depois a gente ajusta)

% --- gnd (nó central)
fwd_model.gnd_node = 17;

% --- eletrodos pontuais (1..16)
for i = 1:16
    fwd_model.electrode(i).nodes = i;      % nó i é o eletrodo i
    fwd_model.electrode(i).z_contact = 0;  % comece com 0 (ajuste depois se quiser)
end



fwd_model.gnd_node = 17; % provisório

[stim, msel] = mk_stim_patterns(16, 1, '{ad}','{ad}', {}, 1);
fwd_model.stimulation = stim;

fwd_model.solve = @fwd_solve_1st_order;
fwd_model.system_mat = @system_mat_1st_order;

%figure; show_fem(fwd_model); axis equal;


fwd_model.jacobian = @jacobian_adjoint;

#img = mk_image(fwd_model, 1.0);  % fundo 1.0 S/m (exemplo)
img = mk_image(fwd_model, 3.0);
J = calc_jacobian(img);   % normalmente retorna M x N (M=medições, N=elementos)
size(J)

s = svd(J);
semilogy(s,'o-'); grid on;
xlabel('índice'); ylabel('valor singular (log)');
title('Espectro singular do Jacobiano');

% rank efetivo com limiar (exemplo)
tol = 1e-6 * s(1);
rank_eff = sum(s > tol);
fprintf('rank efetivo ~ %d (tol=%g)\n', rank_eff, tol);

%figure; show_fem(img); axis equal;

phys = S.phys;                     % Ne x 1
img.elem_data(phys == 1000) = 3.0;  % body
img.elem_data(phys == 1001) = 3.0;  % anomalia1 (exemplo)
img.elem_data(phys == 1002) = 3.0;  % anomalia2
img.elem_data(phys == 1003) = 3.0;  % anomalia3


vh  = fwd_solve(img);
disp('vh.meas (tensões/medições do forward):');
disp(vh.meas);
plot(vh.meas), grid on
title('Medições simuladas')

figure; show_fem(img, 3.0); axis equal;
hold off

imgi = mk_image(fwd_model, 3.0);
#imgi.elem_data(100) = 2.0;  % só pra testar (depois fazemos uma região)


imgi.elem_data(phys == 1000) = 3.0;  % body
imgi.elem_data(phys == 1001) = 2.0;  % anomalia1 (exemplo)
imgi.elem_data(phys == 1002) = 2.0;  % anomalia2
imgi.elem_data(phys == 1003) = 2.0;  % anomalia3
#show_fem(imgi);
vi = fwd_solve(imgi)

% --- INV MODEL (criar aqui, imediatamente antes de usar)
inv_mdl = struct();
inv_mdl.type = 'inv_model';
inv_mdl.name = 'inv_model_circ16_importado';
inv_mdl.fwd_model = fwd_model;

disp("DEBUG inv_mdl fields:");
disp(fieldnames(inv_mdl));

fwd_model.jacobian = @jacobian_adjoint;
inv_mdl.fwd_model  = fwd_model;
inv_mdl.reconst_type = 'difference';


inv_mdl.solve = @inv_solve_diff_GN_one_step;

% jacobiano (OBRIGATÓRIO)
inv_mdl.jacobian = @jacobian_adjoint;

#inv_mdl.RtR_prior = @prior_laplace;
inv_mdl.RtR_prior = @prior_gaussian_HPF

inv_mdl.hyperparameter.value = 1e-3;
inv_mdl.jacobian_bkgnd.value = 1.0;

% --- checagens para não dar o mesmo erro de novo
disp(fieldnames(inv_mdl));
if ~isfield(inv_mdl,'jacobian')
    error("inv_mdl ainda está sem 'jacobian' (erro no script).");
end


% --- reconstrução


img_rec = inv_solve(inv_mdl, vh, vi);
show_fem(img_rec, 3.0);
