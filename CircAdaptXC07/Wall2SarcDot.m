function Wall2SarcDot
% function Wall2SarcDot
% After multiple patch geometry has been solved, for walls and patches
% stresses and areas are calculated for output. Needed to calculate the
% derivatives of sarcomere state variables.
% Theo Arts, Maastricht University, March 9, 2024

global P;

T         = P.Wall.T*P.Wall.Wall2Patch;
Ap        = P.Patch.Ap0 + T .* P.Patch.DADT; % Patch area
P.Patch.T = T ; % Patch tension
P.Patch.Ap= Ap; % Midwall patch area
P.Patch.Ef= 0.5*log(Ap./P.Patch.ApRef); % Sarcomere strain

SarcDot; % Stress Sf in patches, time derivatives XbDot and LsiDot
end
