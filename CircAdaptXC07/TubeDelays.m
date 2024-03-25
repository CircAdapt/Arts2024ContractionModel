function TubeDelays
% function TubeDelays
% Waves propagate on the basis of transmural pressures
% Reflections are determined by proximal and distal flow.
% Input: flows and pressures
% Output: wave amplitudes at both starting sides to fill the delay lines
% Zero-flow pressures at both sides, including estimations for the next
% time sample.
% Theo Arts, Maastricht University, June 8, 2018

global P
if P.Tube.n==0 || size(P.t,1)>1
    return
end

dt  = P.General.dt; % time step
t   = P.t; % current time
tEnd= size(P.Tube.pL,1)*P.General.dt; % maximum delay

Tube= P.Tube; % copy Tube record for local use
nt  = size(Tube.uP,1); % number of time samples in this beat
TauL= Tube.TauL; % row of delay time intervals L-wave
TauR= Tube.TauR; % row of delay time intervals R-wave

cL  = Tube.cL ; % wave velocity to left (backward)
cR  = Tube.cR ; % wave velocity to right (forward)
Att = Tube.Att; % Womersley attenuation (see TubeV2p)
Len = Tube.Len; % Tube length

pTrans= Tube.pTrans; % transmural tube center pressure
qProx = Tube.qProx ; % proximal inflow
qDist = Tube.qDist ; % distal outflow
pSProx= Tube.pSProx; % proximal zero-flow source pressure
pSDist= Tube.pSDist; % distal zero-flow source pressure
pL    = Tube.pL    ; % pL= pLWave*(ZR+ZL)/ZL
pR    = Tube.pR    ; % pR= pRWave*(ZR+ZL)/ZR
ZL    = Tube.ZL    ; % leftward wave impedance
ZR    = Tube.ZR    ; % rightward wave impedance

it  = round(t/dt+1)    ; % time index of Tube-related array
it2 = Mod1(it+[0;1],nt); % cyclically indexed pair it+[0,1] of time
q   = Tube.q(it2(1),:) ; % 'DC' tube flow
% New estimate of delay time taking while avoiding time reversal
TauL= TauL + dt*( 1 - TauL.*cL./Len ); %delay increment left wave
TauR= TauR + dt*( 1 - TauR.*cR./Len ); %delay increment right wave
TauL= max(dt,min(0.9*tEnd,TauL));% limitation for extremely long delay
TauR= max(dt,min(0.9*tEnd,TauR));% limitation for extremely long delay
% Modulo length time interval = tEnd calculation
tL  = mod(t-TauL,tEnd); % delay leftward wave
tR  = mod(t-TauR,tEnd); % delay rightward wave
% Initiation of pressures wave, relative to external tube pressure
pP  = Fu(tL, pL, dt); % delayed pL by interpolation, DC included
pD  = Fu(tR, pR, dt); % delayed pR by interpolation, DC included
% Womersley determined wave attenuation
AttL= exp(-TauL.*Att); % wave attenuation
AttR= exp(-TauR.*Att); % wave attenuation

pTrans= max(0,pTrans) ; % Transmural wave pressure cannot be negative
pP0   = pTrans - q.*ZR; % 'DC'-component is not attenuated
pD0   = pTrans + q.*ZL; % 'DC'-component is not attenuated
ZSum  = ZL+ZR;
% Initial left and right pressure waves
pL    = pSDist - qDist.*ZSum;% pL= pLWave*(ZR+ZL)/ZL
pR    = pSProx + qProx.*ZSum;% pR= pRWave*(ZR+ZL)/ZR

% pL and pR, entering the tube, are stored. At the time the wave
% arrives at the other side, the stored wave pressure is determined by
% linear interpolation, using the current time minus the delay time for
% the wave in the tube. The wave is attenuated according to Womersley,
% but the 'DC' flow and pressure component is not attenuated.

% 'DC' component of flow. If there is no wave and no pressure drop, but
% a flow q, the difference between proximal and distal zero-flow
% pressure equals pSDist-pSProx, and q=(pSDist-pSProx)/(ZR+ZL)

dq= 0.5*dt*((-pSProx+pSDist)./ZSum - q)./(TauL+TauR); %DC-flow increment
% Estimate of q, pL and pR for current and the first future sample.
% Maybe neeeded for very short time delays.
Tube.q(it2,:) = [1;1]*(q+dq); % estimate DC flow at current time sample
Tube.pL(it2,:)= [1;1]*pL   ; % time sample, and also used as estimate
Tube.pR(it2,:)= [1;1]*pR   ; % in upcoming time sample

Tube.TauL = TauL; % store new delay
Tube.TauR = TauR; % store new delay
% The dynamic wave is attenuated, but the 'DC' component is not
% pP and pD are estimated pressure wave amplitudes, arriving at the
% other side, i.e., after the wave delay.

Tube.uP(it2,:)=[1;1]*((AttL.*(pP-pP0)+pP0)./sqrt(ZR)); % also used for estimate
Tube.uD(it2,:)=[1;1]*((AttR.*(pD-pD0)+pD0)./sqrt(ZL)); % in upcoming time sample

P.Tube= Tube;

end

function j=Mod1(i,nt)
% generates array index by modulo interval 1:nt
j=mod(i-1,nt)+1;
end

function f=Fu(t,F,dt)
% Linear interpolation of matrix F along columns
% Matrix t, width the same as F, along column directions t-value,
% Serving as basis for F interpolation.
% Correspondence: i=1 -> t=0
% Output matrix f = interpolated with size(f)=size(t)
[nt,nc] = size(F);
F(end,:)= F(1,:)   ; % Circularity of 1st with last array element
nd      = size(t,1); % number of rows of t
it      = mod(t/dt,nt-1); % non-integer time sample shift
jt      = it+nt*(0:nc-1)+1; % single column index
f       = interp1(F(:),jt(:)); % single column interpolation
f= reshape(f,[nd,nc]); % reshaping to match size t-matrix and output
end
