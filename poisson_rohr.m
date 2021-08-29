% Programm zur Lösung der Poissongleichung auf einem Rohr mit
% Dirichlet-Nullrandwerten und ausführlichen Erklärungen

%% Geometrie des Rohrs definieren
curve1 = nrbcirc (1);
curve2 = nrbcirc (2);
ring  = nrbruled (curve1, curve2);
pipe  = nrbextrude (ring, [0 0 5]);
pipe = nrbdegelev(pipe, [2,2,2]);

%% Rechte Seite f der Poisson-Gleichung 
f = @(x,y,z) (8-9*sqrt(x.^2+y.^2)).*sin(2*atan(y./x))./(x.^2+y.^2);

%% Aufstellen des LGS

% Einlesen der NURBS-Geometrie
geometry = geo_load(pipe);

% Definition der Quadraturregel (Gauß-Legendre-Quadratur)
rule    = msh_gauss_nodes(geometry.nurbs.order);
[qn,qw] = msh_set_quad_nodes(geometry.nurbs.knots, rule); 

% Mesh-Struktur, die Informationen zur Geometrie und Quadratur enthält
msh = msh_cartesian (geometry.nurbs.knots, qn, qw, geometry);

% Definition des diskreten NURBS Raumes V_h 
space = sp_nurbs (geometry.nurbs, msh); 

% Steifigkeitsmatrix und Lastvektor
K = op_gradu_gradv_tp (space, space, msh, @(x,y,z) ones(size(x))); 
rhs = op_f_v_tp (space, msh, @(x,y,z) f(x,y,z));

%% Dirichlet-Nullrandbedingungen fast überall
drchlt_dofs = [];
% Dirichlet-Randbedigungen auf fast dem gesamten Rand
for i = [1,2,3,5,6]
    drchlt_dofs = union(drchlt_dofs, space.boundary(i).dofs);
end
int_dofs = setdiff([1:space.ndof],drchlt_dofs);

%% Zusätzlich: inhomogene Neumann-Randbedigungen auf Mantelfläche
G = zeros (space.ndof, 1);
g = @(x,y,z) sin(y).*cos(x)-3*x.^2
G(space.boundary(4).dofs) = op_f_v_tp (space.boundary(4), msh.boundary(4), @(x,y,z) g(x,y,z));

% Kommentiert man diesen Teil aus, also verzichtet auf die Neumann-Randbedingungen, müsste man in Zeile 26 noch den Index 4
% mit in den Dirichlet-Nullrandteil einbauen und statt Zeile 42 Zeile 43
% als rechte Seite verwenden.

%% Lösen des LGS
u = zeros(space.ndof,1);
u(drchlt_dofs) = 0;
rhs = rhs + G - K*u;
% rhs = rhs - K*u;
u(int_dofs) = K(int_dofs, int_dofs)\rhs(int_dofs);

%% Exportieren der Lösung
output_file = 'pipe';
vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20), linspace(0, 1, 20)};
fprintf ('The result is saved in the file %s \n \n', output_file);

% Export als VTS Datei
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')