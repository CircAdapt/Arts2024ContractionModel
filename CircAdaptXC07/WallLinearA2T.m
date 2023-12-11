function WallLinearA2T
% function WallLinearA2T
% Determines linear function patch area Ap = fu Tension T
% Ap(T)=Ap0+T*DADT
% Constants Ap0 and DADT are determined from Patch state variables
% and reference geometry
% Theo Arts, Maastricht University, April 15, 2018

global P;

% Patch Ap= Ap0+DADT*T, provides Ap0 and DADT
ApRef     = P.Patch.ApRef; % midwall area for ref sarcomere length 2mu
VWall     = P.Patch.VWall; % wall volume
hi        = 0.5*P.Patch.Lsi; %extension factor relative to Ls=2um
Wall2Patch= P.Wall.Wall2Patch; %linking matrix Wall-Patch
ApDead    = P.Wall.ApDead; % non-contractile, stiff wall area
hSe       = P.Patch.hSe;

Ap        = hi.^2 .* ApRef; % midwall patch area
Ef        = log(hi); % natural fiber strain
P.Patch.Ef= Ef; % fiber strain Ef with zero length SE-element
SarcSf; % sarcomere strain->stress function
Sf        = P.Patch.Sf; % sarcomere stress
DSfDEf    = P.Patch.DSfDEf; % stiffness

DEf0 = -Sf./DSfDEf; % zero tension strain relative to hi
DEfSe= hSe./hi; % series elastic element strain
Corr = 1+(2*DEfSe-DEf0)*1.16; % non-linearity correction for stiffness
Ap0  = Ap.*exp(2*DEf0); % zero tension area
DADT = 4*Corr.*Ap.^2./(DSfDEf.*VWall);
% Correction DADT linearizes tension-area relation. A as fu(T) is
% non-linear. Corr shifts best fit to the mid of the working range

P.Patch.DADT= DADT; % area compliance patches
P.Patch.Ap0 = Ap0; % zero tension area patches

% Wall is composed of patches: Also for wall: Aw(T)=Aw0+DADT*T
P.Wall.VWall= VWall*Wall2Patch'; % wall volume
P.Wall.Aw0  = ApDead + Ap0*Wall2Patch'; %zero tension area
P.Wall.DADT = DADT*Wall2Patch'; %wall compliance

end

