function TriSegV2p
% function TriSegV2p
% TriSeg is a 3-wall structure (Left,Septal,Right) with 2 cavities (R,L)
% Given the state variables cavity volumes VL and VR -> 
% Calculates: dimensions of the 'double bubble'
% myofiber stress Sf, wall tension Tm and cavity pressures pTrans
% Quantities V and Y are solved from given VL and VR
% V = septal cap volume
% Y = junction radius.
% State variables P.TriSeg V and Y represent initial estimates of
% VS and YS. Only 1 or 2 iterations have been applied.
% Theo Arts, Maastricht University, nov 17, 2022
% (last correction: factor pi in dV, line 163)

global P;
if P.TriSeg.n==0 %if there is no TriSeg
    return
end

TriSeg = P.TriSeg;
Wall   = P.Wall;
Tau    = P.General.dt  ; % lowpass for VS,YS 1st estimate
RhoB   = P.General.RhoB; % blood density
VLT    = TriSeg.VL     ; % left cavity volume
VRT    = TriSeg.VR     ; % right cavity volume
n      = TriSeg.n      ; % number of TriSeg's
iWall  = TriSeg.iWall  ; %related walls
VT     = TriSeg.V      ; % init septal cap volume, to be solved
YT     = TriSeg.Y      ; % init radius junction circle, to be solved
Aw0W   = Wall.Aw0      ; %zero stress wall area
DADTW  = Wall.DADT     ; % wall compliance
VWallW = Wall.VWall    ; % 3 wall volumes

for iTr=1:n %for all TriSeg's
    iW   = iWall(:,iTr)  +(0:2); % 3 walls
    Aw0  = Aw0W(:,iW) ;      % zero stress wall area
    DADT = DADTW(:,iW);      % wall compliance
    VWall= VWallW(iW) ;      % 3 wall volumes
    VWL  = mean(VWall(1:2)); % enclosed wall volume LV cavity
    VWR  = mean(VWall(2:3)); % enclosed wall volume RV cavity
    VL   = VLT(:,iTr)+VWL; % midwall enclosed LV-volume
    VR   = VRT(:,iTr)+VWR; % midwall enclosed RV-volume
    VS   = VT(:,iTr); % septal rightward volume shift
    YS   = YT(:,iTr); % radius LV-RV junction
    
    VS0=VS; % storage of initial value septal cap volume
    YS0=YS; % storage of initial value junction radius
    
    % Calculation iteration VS and YS increment for better equilibrium 
    for j=1:1 % 1-2 iterations appear to be sufficient
        [dVS,dYS]= DvDp(VS,YS,VL,VR,Aw0,DADT);
        VS= VS+dVS; % Improved estimate of TriSeg geometry
        YS= YS+dYS; % only one iteration has been applied
    end
    
    % Limitation of time-derivatives to min and max value
    % for better numerical safety
    AL = Aw0(:,1); %zero stress wall area L
    AS = Aw0(:,2); %zero stress wall area S
    AR = Aw0(:,3); %zero stress wall area R
    F1 = 1.08; % maximum allowed stretch from zero load state
    F3 = F1^3;
    
    YLo= sqrt(AS./(1+AS./min(AL,AR))/pi)/F1; % assume sphere AL+AS
    YHi= F1*sqrt(AS./(1+AS./(max(AL,AR)*F1^2))/pi);
    VLo= -F3*(2*pi/3)*(AL/(2*pi)).^1.5+VL; % given wall area -> max volume
    VHi=  F3*(2*pi/3)*(AR/(2*pi)).^1.5-VR; % given wall area -> max volume
    YS = max(YLo,min(YHi,YS)); % avoid too large steps, numerical safety
    VS = max(VLo,min(VHi,VS));
       
    % final TriSeg geometry requires some recalculations
    VM   = [VS-VL,VS,VS+VR]; %1st estimate cap-volumes
    ARef = pi*YS.^2        ; % normalization reference area
    V    = VM./(ARef.*YS)  ; % normalized capvolumes
    A0   = Aw0./ARef       ; % normalized zero stress area
    dTdA = ARef./DADT      ; % area normalized wall stiffness
    % Auxilary normalized variables
    X    = VdPi2X(V); % normalized cap-height
    XX   = X.^2;
    RR   = XX+1; % normalized wall area
    C    = 2*X./RR; % normalized wall curvature
    
    % Revert normalization
    T     = max(0,(RR-A0).*dTdA);% [N/m] wall tension
    Aw    = ARef.*RR; % wall area
    Cm    = C./YS; % wall curvature
    pTrans= 2*Cm.*T; % transmural pressure with effect of buckling
    
    % Cavity wave impedance properties, needed to make node connection
    Vm   = [VL+VWL,VR+VWR];
    Len  = 2*Vm.^(1/3); % cavity length
    A    = Vm ./ Len; % cross-sectional area
    Z0   = sqrt(RhoB./(4*DADT(:,[1,3]).*abs(A).^1.5)); %Wave impedance/2

    % anti-collapse counter pressure for numerical safety
    p0  = 0.2*P.General.p0;
    eps = 0.1;
    VNLo= max(eps,[VL,VR]./[VWL,VWR]-1);
    dpLo= p0.*max(0,0.3./VNLo-1).^2;% anti-collapse safety

    % writing wall data
    Wall.Aw(:,iW)    = Aw ; % wall area
    Wall.Cm(:,iW)    = Cm ; % wall curvature
    Wall.T(:,iW)     = T  ; % wall tension
    Wall.pTrans(:,iW)= pTrans-mean(pTrans,2); % pTrans
    
    % State variable derivatives
    TriSeg.VDot(:,iTr)= (VS-VS0)/Tau; % V serves as initial estimate of VS
    TriSeg.YDot(:,iTr)= (YS-YS0)/Tau; % Y serves as initial estimate of YS
   
    % writing geometric data TriSeg to be used for output, not for solving
    TriSeg.VS(:,iTr)= VS     ; % final solution Rv-Sv-Lv
    TriSeg.YS(:,iTr)= YS     ; % volume of septal cap, only used for output
    TriSeg.AL(:,iTr)= A(:,1) ; % cross-sectional area LV
    TriSeg.AR(:,iTr)= A(:,2) ; % cross-sectional area RV
    TriSeg.ZL(:,iTr)= Z0(:,1); % source impedance LV
    TriSeg.ZR(:,iTr)= Z0(:,2); % source impedance RV
    TriSeg.pTransL(:,iTr)= -pTrans(:,1)-dpLo(:,1); % LV transmural pressure
    TriSeg.pTransR(:,iTr)=  pTrans(:,3)-dpLo(:,2); % RV transmural pressure
end

P.Wall  = Wall;
P.TriSeg= TriSeg;
end

%=== Additional function ============ 
function [dV,dY]=DvDp(VS,YS,VL,VR,Aw0,DADT)
% input: VS=septal cap volume, YS= junction radius
% [VL, VR]= [left,right] midwall enclosed volumes
% Aw0,DADT= [zero load wall area, compliance wall] for L- S- R- walls
VM   = [VS-VL,VS,VS+VR]; % cap volumes LSR
ARef = pi*YS.^2; % reference area for normalization
V    = VM./(ARef.*YS); %normalized cap volumes
A0   = Aw0./ARef     ; %normalized wall areas
dTdA = ARef./DADT    ; % normalized wall compliance

X    = VdPi2X(V); % solves 3rd order polynomial analytically
XX   = X.^2;
RR   = XX+1; %wall area (normalized to junction area (pi y^2)
RR2  = RR.^2;

dA = RR-A0; % wall area stretch increment
T  = max(0,dA.*dTdA);% buckle possibility

Av =  4*X  ./RR; % partial derivative dA/dVS, normalized
Ap = (1-XX)./RR; % partial derivative dA/dYS, normalized

Avv= 8*Ap./RR2; % 2nd order partial derivatives
Avp= -2*Av./RR2;
App= 3*V.*Av./RR2;

Fv= sum(T.*Av,2); %2*Tx, dE/dVS ~force
Fp= sum(T.*Ap,2); %-Ty, dE/dYS ~force

Fvv=sum(dTdA.*Av.^2 +T.*Avv,2); % 2nd derivative,'stiffness'
Fvp=sum(dTdA.*Av.*Ap+T.*Avp,2); % 2nd derivative,'stiffness' 
Fpp=sum(dTdA.*Ap.^2 +T.*App,2); % 2nd derivative,'stiffness'

Det= Fvv.*Fpp-Fvp.*Fvp; % determinant
dv= -( Fv.*Fpp-Fp.*Fvp)./Det; %solution of 2 parameters by matrix inversion
dp= -(-Fv.*Fvp+Fp.*Fvv)./Det;
dV= pi  * YS.^3 .* dv; % reversion of normalization, correction pi WB16
dY= 0.5 * YS    .* dp;

end

function X=VdPi2X(VdPi)
% Analytical solution of X as a function of cap volume V
% VdPi = normalized Vcap/(pi*Rjunction^3)
% assumption: cap junction radius is normalized to 1
% Solving 3rd order polynomial
% Mathematica: Solve[x^3 + 3x  - 2V == 0, x]
% Q= (V + Sqrt(V^2 + 1))^(1/3);  x= Q - 1/Q;
V = 3*VdPi;
Q = (V + sqrt(V.^2+1)).^(1/3);
X = Q - 1./Q;
end