function WallLinearA2T
% function WallLinearA2T
% Input: sarcomere length at zero Xb tension (Lsi= 2 hi)
% Determines linear function for patch area Ap(T)=Ap0+T*DADT around
% isometric contraction state
% Output: Ap0 and DADT
% Theo Arts, Maastricht University, March 9, 2024

global P;

% Patch Ap= Ap0+DADT*T, provides Ap0 and DADT
ApRef     = P.Patch.ApRef; % midwall area for ref sarcomere length 2mu
VWall     = P.Patch.VWall; % wall volume
hi        = 0.5*P.Patch.Lsi; % extension factor relative to Ls=2um
Wall2Patch= P.Wall.Wall2Patch; %link matrix Wall<->Patch
ApDead    = P.Wall.ApDead; % non-contractile, stiff wall area
hSe       = P.Patch.hSe; % series elasticity

Ap        = hi.^2 .* ApRef; % midwall patch area
Ef        = log(hi); % natural fiber strain
P.Patch.Ef= Ef; % fiber strain Ef with zero length SE-element
SarcSf; % sarcomere strain->stress function
Sf        = P.Patch.Sf; % sarcomere stress
DSfDEf    = P.Patch.DSfDEf; % stiffness

DEf0 = -Sf./DSfDEf; % zero tension strain relative to hi
DEfSe= hSe./hi; % series elastic element strain
Ap0  = Ap.*exp(2*DEf0); % zero tension area
SfIso=Sf+DEfSe.*DSfDEf; %Isometric tension, selected center of work region
DADT = 4.*Ap.^2./((DSfDEf-2*SfIso).*VWall); %derivative at isomitric Ls

% Correction DADT linearizes tension-area relation. A as fu(T) is
% non-linear. Corr shifts best fit to the mid of the working range

P.Patch.DADT= DADT; % area compliance patches
P.Patch.Ap0 = Ap0; % zero tension area patches

% Wall is composed of patches: Also for wall: Aw(T)=Aw0+DADT*T
P.Wall.VWall= VWall*Wall2Patch'; % wall volume
P.Wall.Aw0  = ApDead + Ap0*Wall2Patch'; %zero tension area
P.Wall.DADT = DADT*Wall2Patch'; %wall compliance

end
