% Programm zur Lösung der Poissongleichung auf dem Einheitskreis mit
% Dirichlet-Nullrandwerten und Fehlerabschätzung

tic
%% Geometrie des Kreises definieren
line = nrbline([0,0,0],[1,0,0]);
kreis = nrbrevolve(line,[0,0,0],[0,0,1]);
kreis = nrbdegelev(kreis, [0,1]);

%% Rechte Seite f der Poisson-Gleichung und exakte Lösung u_ex
f = @(x,y) - y.^2./((x.^2+y.^2).^(3/2)) - x.^2./((x.^2+y.^2).^(3/2));
u_ex = @(x,y) sqrt(x.^2+y.^2)-1;

%% Aufstellen des LGS
geometry    = geo_load(kreis);
rule        = msh_gauss_nodes(geometry.nurbs.order);
[qn,qw]     = msh_set_quad_nodes(geometry.nurbs.knots, rule); 
msh         = msh_cartesian (geometry.nurbs.knots, qn, qw, geometry);
space       = sp_nurbs (geometry.nurbs, msh); 
K           = op_gradu_gradv_tp (space, space, msh, @(x,y) ones(size(x))); 
rhs         = op_f_v_tp (space, msh, @(x,y) f(x,y));

%% Dirichlet-Nullrandbedingungen auf Kreisrand
drchlt_dofs = [];
for i = 4
    drchlt_dofs = union(drchlt_dofs, space.boundary(i).dofs);
end
int_dofs = setdiff([1:space.ndof],drchlt_dofs);

%% Lösen des LGS
u = zeros(space.ndof,1);
u(drchlt_dofs) = 0;
rhs = rhs - K*u;
u(int_dofs) = K(int_dofs, int_dofs)\rhs(int_dofs);

%% Plot der Lösung
output_file = 'kreis';
vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u');

%% Berechnung des Fehlers in der L^2-Norm
l2_err = sp_l2_error (space, msh, u, @(x,y) u_ex(x,y))
toc