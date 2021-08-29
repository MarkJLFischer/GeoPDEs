% Programm zur Erstellung des Doppel-T-Trägers und der Einheitskugel

%% Doppel-T-Träger
pnts_r =  [7 7 4 4 7 7; 0 -1 -1 -4 -4 -5];
right = nrbmak(pnts_r, [0 0 1/5 2/5 3/5 4/5 1 1]);

pnts_l = [0 0 3 3 0 0; 0 -1 -1 -4 -4 -5];
left = nrbmak(pnts_l, [0 0 1/5 2/5 3/5 4/5 1 1]);

srf = nrbruled(left, right);
vol = nrbextrude(srf, [0 0 7]);
nrbctrlplot(vol)

%% Einheitskugel
% crv1 = nrbcirc(1, [], 0, pi);
% crv2 = nrbline([-1 0 0], [1, 0 ,0]);
% srf = nrbruled(crv1, crv2);
% vol = nrbrevolve(srf, [0 0 0], [1 0 0]);
% nrbctrlplot(vol)