function ValveqDot
% function ValveqDot
% Node pressure differences -> valve flow acceleration qDot
% Valves AvValves refer to ventricular valves with papillary muscles
% Improvement 04-22-2023: avoiding discontinuities in valve closure
% Better suited for DE-solution with non-fixed time step.
% Theo Arts, Maastricht University, Sep 26, 2023

global P

iNodeProx= P.Valve.iNodeProx;
iNodeDist= P.Valve.iNodeDist;
q        = P.Valve.q; % flow
nt       = size(q,1);
AOpen    = repmat(P.Valve.AOpen,[nt,1]); % open valve cross-section
ALeak    = repmat(P.Valve.ALeak,[nt,1]); % closed valve cross-section
Len      = P.Valve.Len; % effective length of flow channel
RhoB     = P.General.RhoB; % density of blood
TauV     = P.General.dt; % Decay time to avoid numerical closure artifacts

% allows AV-valve diastolic regurgitation
% may need further improvement
Ws   = P.Valve.AvWalls;
Vs   = P.Valve.AvValves;
T    = P.Wall.T(:,Ws); %wall tension related to pap. muscle
DADT = P.Wall.DADT(:,Ws); %wall stiffness
Aw0  = P.Wall.Aw0(:,Ws); %wall zero-stress area
Diast= tanh(30*max(0,T.*DADT./Aw0-0.1).^2); % diastole->1.0
ALeak(:,Vs)= 0.3*(AOpen(:,Vs)-ALeak(:,Vs)).*Diast+ALeak(:,Vs); 

pDrop = P.Node.p(:,iNodeProx)-P.Node.p(:,iNodeDist); % pressure drop

AMax  = max(AOpen,ALeak); % Valve cross-sectional area
AProx = P.Node.A(:,iNodeProx)-AMax; % prox summed area of node connections
ADist = P.Node.A(:,iNodeDist)-AMax; % dist summed area of node connections
AOpen = min(AOpen,max(AProx,ADist));% Avoid AOpen> all other A's
ALeak = min(ALeak,max(AProx,ADist));% Avoid AOpen> all other A's

Sq   = sign(q); % flow direction
Rhov3= 1e-3*RhoB; % energy treshold for valve opening/closure 
Avd2 = (AOpen-ALeak)/2;
Asd2 = (AOpen+ALeak)/2;
v    = 0.5*q./(Sq.*Avd2+Asd2);
z    = max(0,pDrop.*v./Rhov3).*Sq+2.0*Avd2./Asd2;
A    = Avd2.*tanh(z)+Asd2;

ADS= max(A,0.5*((ADist+AProx)-Sq.*(ADist-AProx)));%Downstream ADS,always >A
% Aux   = 0.5*(A.^-2 - ADS.^-2);
Aux   = 0.78*(1./A - 1./ADS).^2; % irreversible loss component only (X926)
R     = RhoB*abs(q).*Aux; % Bernouilli resistance
pDropB= q.*R; % Bernouilli pressure drop
Len1  = sqrt(Len.^2+A+(TauV*q./A).^2);%freq. limitation by length increase
L1    = 1.33*RhoB*Len1./A; % Valve inertia only
qDot  = (pDrop-pDropB)./L1; %non-lin Bernoulli solution

P.Valve.qDot= qDot; % flow derivative
end

