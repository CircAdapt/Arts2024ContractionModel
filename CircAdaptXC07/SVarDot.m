function OutDotT=SVarDot(SVarT)
% function OutDotT=SVarDot(SVarT)
% SVarT = column vector= transpose of state variables SVar vector
% OutDotT= column vector of derivatives
% State Variables -> Time derivatives of State Variables
% Ready to be used in MatLab function odeCA() to solve set of Diff Eq
% Theo Arts, Maastricht University, Jun 8, 2018
%====

global P

P.SVar= SVarT'; % store state variables SVar
SVar2P; % state variables SVar -> physiologic representation P.xx
P.tDot=ones(size(P.t)); % time derivative of time = 1

MemAlloc     ; % necessary memory allocations
ArtVenV2p    ; % arteries to veins compliant network
WallLinearA2T; % patch and wall: Am= Am0+T*DADT
ChamberV2p   ; % Chamber, pTrans and Wall.Am,T
TriSegV2p    ; % TriSeg, pTrans and Wall.Am,T
Wall2SarcDot ; % filling Sarc with LsiDot,XbDot
TubeV2p      ; % Delayed source transmural pressures and resistances of tube
pBag         ; % Bag pressures render external pressures on elements
pNodeVDot    ; % Node pressures and volume derivatives VDot's
ValveqDot    ; % flow derivatives qDot
TubeDelays   ; % renders wave amplitude just at moment of reflection, delays

P2SVarDot   ; % transfer of derivatives to P-structure

OutDotT= real(P.SVarDot)'; % odeXX requires OutDotT to be a column vector
end
