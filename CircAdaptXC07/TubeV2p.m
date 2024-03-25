function TubeV2p
% function TubeV2p
% Tube volume V -> transmural pressure pTrans, wave propagation velocity
% c0, proximal and distal zero-flow pressure pSProx, pSDist and source
% impedances ZR (prox), ZL (dist), using delayed pressures.
% Theo Arts, Maastricht University, March 11, 2024

global P

dt    = P.General.dt   ; % integration timestep
it    = round(P.t/dt)+1; % time sample counter
q     = P.Tube.q(it,:) ; % lenght averaged tube flow
Len   = P.Tube.Len     ; % representative length of blood vessels
p0    = P.Tube.p0      ; % working pressure
A0    = P.Tube.A0      ; % tube cross-section at p0
Aw    = P.Tube.AWall   ; % wall cross-section
k     = P.Tube.k       ; % tube stiffness coefficient
RhoB  = P.General.RhoB ; % blood density
fWom  = 2.5/P.General.tCycleRef;% applied major Womersley frequency
EtaB  = P.General.EtaB ; %blood viscosity

% Transmural pressure and wave velocity/impedance= fu(cross-section)
A     = P.Tube.V./Len; % vessel cross-section
hRPois= 4*pi*EtaB*Len./A.^2; % half Poiseuille resistance
m     = k/3-1; % cross-sectional stiffnes parameter
a     = (A0./Aw).^0.3; % a thick wall (small a) counteracts vessel collapse
am    = a.*m;
ARef  = A0./am; % reference A, A=ARef for p=0;
fa    = am.^m;
pRef  = p0./(fa-1./fa); % reference p
AN    = A./ARef; % normalized A
f     = AN.^m;
pN    = f-1./f;
KAN   = m.*(f+1./f);
pTrans= pN.*pRef; % transmural pressure
c0    = sqrt(KAN.*pRef/RhoB); % zero flow wave velocity

%Approximation of Womersley attenuation with aWom= r Sqrt(rho w/eta)
aWom  = sqrt(A*(2*fWom*RhoB/EtaB)); % Womersley number
c0    = c0./(1+0.71./aWom+0.71./aWom.^4); % Womersley correction
Att   = fWom*6.3*sqrt(1+0.03125*aWom.^2)./(1+0.25*aWom.^2); % attenuation
Z0    = RhoB * c0 ./ A; % corrected wave impedance with tube flow=0

% Flow dependency of wave velocity and impedance
% getting delayed signals for zero-flow pressure pSProx and pSDist
vb    = q./A; % mean blood velocity
hvdc0 = 0.5*vb./c0; % blood velocity, normalized to wave velocity
% Flow velocity correction on wave velocity and wave impedance
b     = sqrt(1+hvdc0.^2);
bR    = b+hvdc0;
bL    = b-hvdc0;
cR    = sqrt(0.25+(c0.*bR).^2);% wave velocity clipped to >0.5 m/s
cL    = sqrt(0.25+(c0.*bL).^2);% wave velocity clipped to >0.5 m/s
ZR    = Z0./bR;
ZL    = Z0./bL;
DpPois= q.*hRPois; % half pressure drop due to DC-flow 

P.Tube.pSProx = P.Tube.uP(it,:).*sqrt(ZR)+DpPois; % zero flow pressure Prox
P.Tube.pSDist = P.Tube.uD(it,:).*sqrt(ZL)-DpPois; % zero flow pressure Dist
P.Tube.ZR     = ZR     ; %R-wave impedance
P.Tube.ZL     = ZL     ; %L-wave impedance
P.Tube.cR     = cR     ; %R-wave velocity
P.Tube.cL     = cL     ; %L-wave velocity
P.Tube.Att    = Att    ; %Womersley attenuation factor
P.Tube.A      = A      ; %Cross-sectional area
P.Tube.pTrans = pTrans ; %Transmural pressure
end
