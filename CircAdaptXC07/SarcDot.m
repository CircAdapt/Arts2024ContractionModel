function SarcDot
% function SarcDot
% Calculates myofiber stress Sf with the linearized stress-strain relation.
% Using the just calculated sarcomere length, normalized series-elastic
% element length L is calculated.
% Time derivatives of the patch-related state variables Lsi and Xb are 
% calculated, using the value of L
% The sarcomere is embedded in Patch
% The contractile mechanism is based on the myosin-MechChem-C/D model
% Theo Arts, Maastricht University, March 9, 2024

global P

t   = P.t; %global time
dt  = P.General.dt;
Sarc= P.Patch;

% Timing reference:  heart rate with exercise
TC  = 2*P.General.tCycle/P.General.tCycleRef; % relative cycle time

%==== Input variables
t1 = Sarc.tDep(1,:); % last depolarization time
t2 = Sarc.tDep(2,:); % depolarization time during current cycle

% Ca-Xb properties
vMx = Sarc.vMx .* TC.^Sarc.avMx ; % 0.5*dLs/dt unloaded sarcomere
hSe = Sarc.hSe ; % [um] half sarc length series elasticity

% Ca Pulse
TauCa = Sarc.TauCa.* TC.^Sarc.aTauCa; % duration systolic Ca-pulse
CaS   = Sarc.CaS  .* TC.^(Sarc.aCaS)  ; % systolic [Ca]-source
CaD   = Sarc.CaD  .* TC.^(Sarc.aCaD)  ; % diastolic [Ca]-source
YCaS  = Sarc.YCaS .* TC.^Sarc.aYCaS ; % systolic pulse Ca conductivity
YCaD  = Sarc.YCaD .* TC.^Sarc.aYCaD ; % diastolic Ca conductivity
Ef    = Sarc.Ef    ; % sarcomere strain
Lsi   = Sarc.Lsi   ; % SVar: unloaded sarcomere length
Xb    = Sarc.Xb    ; % [-] ~ number of Xb's
Sf0   = Sarc.Sf    ; % stress with zero SE-length
DSfDEf= Sarc.DSfDEf; % total stiffness

% tc: time after moment of depolarization of the sarcomere
tc   = t-t1 + (t > t2).*(t1-t2);

% Sarcomere length Ls and stress Sf
hi   = Lsi/2; % half unloaded sarc length
Efi  = log(hi); % sarcomere strain
EfMin= Efi-Sf0./DSfDEf; % Ef with Sf=0, Ef<EfMin -> buckling occurs
EfS  = max(Ef,EfMin); % sarcomere strain, also valid with buckling
hS   = exp(EfS); % half sarc length
Sf   = Sf0+(EfS-Efi).*DSfDEf; %sarcomere stress (linearized representation)
% Xb-force is proportional with L
L    = (hS-hi)/hSe; % norm LenSe, L=1->isometric
a    = 0.2; L=a*log(1+exp(L/a)); % avoid L<0 (=buckling)

Ca   = S2Ca(Xb,hS,L); % sarcomeric equilibrium [Ca]

% Systolic Ca conductivity as Fu(t)
T    = max(0,tc./TauCa); % norm. pulse duration
a    = 10; % Steepness of conductivity pulse with unit amplitude
fu   = (t>0)./(1+exp(-a*(1-T)));

% Ca current
YCaSt= YCaS.*fu; % Systolic Ca conductivity pulse
qD   = (CaD-Ca).*YCaD;
qS   = (CaS-Ca).*YCaSt;
qCa  = qD+qS; % 'raw' XbDot

% Extra Ca binding to Tn at low [Ca], while Xb formation is low
qCa=qCa./(exp(-0.5*Xb)+20*exp(-20*Xb)+1);% correction&safety for low Xb

% Stiff ODE correction by soft limitation of XbDot
y1=  0.5*qCa;
y2= -0.3*Xb/dt; % lower limit of derivative, avoid Xb<0
s = y1+y2;
v = y1-y2;
XbDot= s+sqrt(v.^2+0.01);
hDot = (L-1).*vMx;

%=== Stress/Ls and stiffness, collected for output
P.Patch.LsiDot= 2*hDot; % h=half sarc length
P.Patch.XbDot = XbDot;
P.Patch.Ls    = 2*hS;
P.Patch.Sf    = Sf;
P.Patch.Ca    = Ca;

end

% Auxiliary functions ====================
function Ca=S2Ca(Xb,hS,L)
global P

Sarc= P.Patch;

% Sarcomere structure, h.. refers to half sarc length
hAct = Sarc.hAct ; % [um] thin filament length
hMyo = Sarc.hMyo ; % [um] Myo-length including bare zone
hBare= Sarc.hBare; % [um] bare zone length, may be extended with inactive
fD   = Sarc.fD; % fraction of D-zone of total Myo-length

%contraction parameters
hRef = Sarc.hRef ; % MechChem length reference [um]
gTit = Sarc.gTit; % Titin-Xb unlocking slope
LnbC = Sarc.LnbC; % logarithm of Ca sensitivity ratio C/D-zone

LXb  = L.*Xb;
ST   = 0.5*(LXb+sqrt(LXb.^2+0.25));% numerical stability for low Xb
hBEff= max(hBare,hS-hAct); % non-Xb zone around M-line
hD   = (hMyo-hBare)*fD; % length D-zone
hOv  = hMyo-max(hBEff,hS-hAct);
hC   = max(0.1,hOv-hD); % length C-zone

GTit= L.*hS.^gTit./hRef; %Titin unshielding factor
HD  = hD.*GTit; % Normalized D-zone length
HC  = hC.*GTit; % Normalized C-zone length

ST= min(ST,HC+HD-1); %tension MUST be less than 100% activation
s1= ST;
s4= X2S(S2X(X2S(S2X(s1)-HC)-LnbC)-HD);
Err0= ST-s1+s4+LnbC; % LnbC causes overlap of C- and D-zone conditions
s10 = s1;
s1  = ST+s4;
for i=1:2 % 2 iterations are sufficient
    s4=X2S(S2X(X2S(S2X(s1)-HC)-LnbC)-HD);
    s2=X2S(S2X(s1)-HC);
    s3=s2-LnbC;
    s4=X2S(S2X(s3)-HD);
    Err=ST-s1+s4+LnbC;% improved/corrected Y324
    dE=Err0-Err;
    eps=1e-5; % avoids zero division near iteration target+++++
    N=dE.^2+2*eps;
    w0=(eps+Err0.*dE)./N;
    w1=(eps-Err .*dE)./N;
    s1A=w0.*s1+w1.*s10;
    s10=s1;
    Err0=Err;
    s1=s1A;
end
LnCa= -log(2*exp(max(1e-6,-s4))-1)/1.65; % log must be real.
LnCa= min(1,max(-3,LnCa));
Ca  = exp(LnCa);
end

function z=X2S(x)
% Basic function backward
a=7;
zLo= log(4/a)+asinh(x/a+1.0);
zHi= x+1+exp(-x);
N  = sqrt((zHi-zLo).^2+2.0);
V  = (tanh(0.83*(x+1))+1)/2 ;
z=N.*V+zLo;
for i=1:5 % 2 iteration for high accuracy
    [xz,dxdz]=S2X(z);
    dx=x-xz;
    dz=-dx./dxdz;
    z=z-dz;
end
end

function [x,dxdz]=S2X(z)
% Basic function forward
P=exp(-z);
Q=1+P+sqrt((P.*(P+2)));
x=-Q-log(P.*Q);
dxdz=Q;
end
