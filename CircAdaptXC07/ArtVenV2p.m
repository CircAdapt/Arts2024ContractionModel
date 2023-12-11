function ArtVenV2p
%function ArtVenV2p
% Art/Ven hemodynamics of an organ or group of organs with peripheral
% resistance in between
% volume V-> transmural pressure pTrans(V)
% wave impedance Z(V) and estimate of cross-sectional area A(V)
% Theo Arts, Maastricht University, April 20, 2023

global P;

RhoB  = P.General.RhoB; % blood density
ArtVen= P.ArtVen    ; % ArtVen structure
nAv   = ArtVen.n    ; % number of ArtVen's
VAr   = ArtVen.VAr  ; % arterial volume, state variable
VVe   = ArtVen.VVe  ; % venous volume, state variable
Len   = ArtVen.Len  ; % representative length/size of ArtVen vessel system
p0Av  = ArtVen.p0Av ; % reference pressure drop of microvasculature
q0Av  = ArtVen.q0Av ; % target flow for adaptation of vessel diameter
p0    = ArtVen.p0(:)'   ; % p0=pTrans(A0) reference pressure Ar and Ve 
A0    = ArtVen.A0(:)'   ; % cross-section Ar and Ve, corresponding with p0
Aw    = ArtVen.AWall(:)'; % vessel wall cross-section
k     = ArtVen.k(:)'    ; % stiffness parameter of fibers in vessel wall

% Ar and Ve are merged as a row vector
Ar   = 1:2:2*nAv;
Ve   = Ar+1;
ArVe = [Ar,Ve]; % convert A0 and p0 to row vectors
A0Row= A0(ArVe) ; % working cross-section
AwRow= Aw(ArVe) ; % wall cross-section
p0Row= p0(ArVe) ; % pressure at crossection A0
kRow = k(ArVe)  ; % wall stiffness
AAr  = VAr./Len ; % actual arterial cross-section
AVe  = VVe./Len ; % actual venous cross-section
A    = [AAr,AVe]; % row vector of cross-sectional areas

% Pressure pTrans and wave impedance Z0
m     = kRow/3-1; % cross-sectional stiffnes parameter
a     = (A0Row./AwRow).^0.3; % low a-value counteracts vessel collapse
% Thick wall (=low a) counteracts vessel collapse
am    = a.*m;
ARef  = A0Row./am; % reference A, A=ARef for p=0;
fa    = am.^m;
pRef  = p0Row./(fa-1./fa); % reference p
AN    = A./ARef; % normalized A
f     = AN.^m;
pN    = f-1./f;
KAN   = m.*(f+1./f);
Z0    = sqrt(RhoB*pRef.*KAN)./A; % wave impedance
pTrans= pN.*pRef; % transmural pressure

% Unraveling vessels to Ar and Ve records
RgA     = 1:nAv; % Ar indices
RgV     = RgA+nAv; % Ve indices
pTransAr= pTrans(:,RgA); % transmural pressures prox/Ar
pTransVe= pTrans(:,RgV); % transmural pressures dist/Ve
ZAr     = Z0(:,RgA); % Ar wave impedance
ZVe     = Z0(:,RgV); % Ve wave impedance
Dp      = pTransAr-pTransVe; % pressure drop microcirculation

% Microcirculatory flow, Rpoiseuille corrected for vessel volume change
FAr= VAr./(Len.*A0(Ar)); % Ar Volume factor
FVe= VVe./(Len.*A0(Ve)); % Ve Volume factor
qAv= Dp.*FAr.*FVe.*(q0Av./p0Av); % Flow ~ Volume^2 (Poiseuille)

ArtVen.q       = qAv     ; % microculatory flow
ArtVen.pTransAr= pTransAr; % transmural pressure prox/Ar compartment
ArtVen.pTransVe= pTransVe; % transmural pressure dist/Ve compartment
ArtVen.ZAr     = ZAr     ; % prox/Ar input wave impedance
ArtVen.ZVe     = ZVe     ; % dist/Ve output wave impedance
ArtVen.AAr     = A(:,RgA); % prox/Ar cross-sectional area
ArtVen.AVe     = A(:,RgV); % dist/Ve cross-sectional area

P.ArtVen=ArtVen;

end

