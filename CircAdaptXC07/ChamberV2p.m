function ChamberV2p
% function ChamberV2p
% A chamber is a cavity encapsulated in a single myocardial wall.
% Input: Chamber volume V
% Output: myofiber stress Sf, wall tension T and cavity pressure p
% using the linearized T(Am) relation
% Theo Arts, Maastricht University, March 7, 2024

global P;
if P.Chamber.n==0 %if there is no chamber
    return
end

Chamber= P.Chamber; % Chamber structure
Wall   = P.Wall   ; % Wall structure

RhoB = P.General.RhoB;
iWall= Chamber.iWall          ; % index of related wall
V    = max(0, Chamber.V)      ; % cavity volumes
VWall= Wall.VWall(iWall)      ; % wall volumes
Aw0  = Wall.Aw0(:,iWall)      ; % zero tension midwall area
DADT = Wall.DADT(:,iWall)     ; % wall compliance
p0   = 0.1*P.General.p0       ; % anti-collapse stiffness parameter

Vm   = V+0.5*VWall            ; % midwall enclosed volume
Cm   = (4/3*pi./Vm).^(1/3)    ; % mid wall curvature
Aw   = (4*pi)./Cm.^2          ; % midwall area

a    = 100; % buckling parameter, decrease of stiffness with buckling
DADT = DADT.*exp(a*max(0,1-Aw./Aw0));%Buckling reduces stiffness
T    = (Aw-Aw0)./DADT; % wall tension

% wall properties
Wall.T(:,iWall) = T    ; % wall tension
Wall.Cm(:,iWall)= Cm   ; % curvature=1/radius
Wall.Aw(:,iWall)= Aw   ; % wall area
pTrans = 2*Cm.*T       ; % transmural pressure with effect of buckling
pTrans(:,iWall)= pTrans; % transmural pressure

% Anti-collapse pressure
eps   = 0.1; % lowest value of VNLo=VCavity/VWall for pressure calculation
VNLo  = max(V./VWall,eps);
dpLo  = p0.*max(0,0.5./VNLo-1).^2;% anti-collapse safety
pTrans= pTrans - dpLo;

% Impedance properties, needed to make node connection
Len= 2*Vm.^(1/3); % cavity length
A  = Vm ./ Len  ; % cross-sectional area
Z0 = sqrt(RhoB./(abs(A).^1.5 .* DADT)); %Cavity wave impedance term
Chamber.A     = A     ; % cross-sectional area for valve
Chamber.Z     = 0.5*Z0; % made quite small
Chamber.pTrans= pTrans;

P.Wall   = Wall;
P.Chamber= Chamber;

end
