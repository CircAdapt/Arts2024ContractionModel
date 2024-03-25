function Adapt
%function Adapt
% Standard beat to beat (absence of) adaptation
% Blood pressure is controlled by change of circulating volume
% Theo Arts, Maastricht University, March 16, 2019

global P

% Systemic blood flow adjusted per ArtVen-element by resistance change
% Control of systemic pressure by adjustment of circulatory blood volume
FbFac  = P.General.AdaptFeedback; % Feedback constant. FbFac==0->no control
p0     = P.General.p0;
q0     = P.General.q0;
ErrStop= P.General.ErrStop; %Stop criterium simulation
In     = P.General.In ;
Out    = P.General.Out;
tEnd   = P.General.tEnd;
pNode  = P.Node.p;
iBaro  = P.Node.iBaro;
Name   = P.ArtVen.Name;
q      = P.ArtVen.q;
qRef   = P.ArtVen.qRef;% ArtVen flow at rest
q0Av   = P.ArtVen.q0Av;
p0Av   = P.ArtVen.p0Av;
kExc   = P.ArtVen.kExc;% coefficient to accomodate flow to exercise

% Finding systemic and pulmonary ArtVen's
jPu=contains(lower(Name),'pu');%Pulm. names contain 'pu'
jSy=~jPu;

%=== Estimation of TargetFlow in systemic ArtVen's for flow control
qSy    = mean(q(:,jSy)); % systemic ArtVen flow
qPu    = mean(q(:,jPu)); % pulmonary ArtVen flow
qRefSy = qRef(jSy); % Systemic reference flow at rest

%      AV press.-drop
q0Ref= sum(qRefSy); % total systemic flow at rest
qT   = qRefSy.*(q0/q0Ref).^kExc(jSy); % estimated target flow per ArtVen
qTSy = q0*qT./sum(qT); %target flows per ArtVen, so that total qSys=q0

% Vasodilation as simulated by a decrease of p0Av
p0AvSy= p0Av(jSy);
p0AvSy= p0AvSy.*abs(qSy./qTSy).^FbFac; 

% Systemic blood pressure control
FacpControl= (mean(pNode(:,iBaro))/p0)^(0.5*FbFac);
q0Av(:,jSy)=qTSy; % save Sy-Target flows in P.ArtVen.q0Av
p0Av(:,jSy)=p0AvSy/FacpControl; % flow autoregulation by resistance change

% Control of Ca conductivity to maintain duration of systole
PatchCaControl

% Flow in Sy and Pu to judge steady state
FlowVec=[sum(qSy),sum(qPu)]; % Systemic/Pu flow
disp(['Flow/q0 for Sys Pu: ', num2str(FlowVec/q0,'%8.4f')]);

% Storage of signal value to judge Adapt-convergence
VecV=[P.ArtVen.VAr,P.ArtVen.VVe,P.Tube.V,P.Chamber.V,P.TriSeg.V];
In =[In ;VecV(  1,:)];
Out=[Out;VecV(end,:)];
% Judging quality of steady state
if size(Out,1)>1
    ErrVec= 1000*log( Out(end,:)./In(end,:) );
    disp(['Stationarity error: ',num2str(round(norm(ErrVec)))] );
    %=== ERROR criterium on flow stationarity
    if norm(ErrVec)< ErrStop
        tEnd=0.5*tEnd;
    end
disp(' ');
end

P.ArtVen.q0Av  = q0Av; % Artven target flow
P.ArtVen.p0Av  = p0Av; % ArtVen resistance control
P.General.tEnd = tEnd; % End time of simulation
P.General.In   = In  ; % Storage of signals to judge stationarity
P.General.Out  = Out ; % Storage of signals to judge stationarity
P.General.FacpControl= FacpControl; % pressure control by volume change

P2SVar; % load adapted values in P.SVar for initiation of next beat
end

function PatchCaControl
global P
Ca =P.Patch.Ca;
dCa=P.Patch.CaS-P.Patch.CaD;
a  =0.2;
CaN=(Ca-P.Patch.CaD)./dCa; % normalized Ca
f=tanh(max(0,CaN)/a);
DutyCycle=sum(f)*(P.General.Dt/P.General.tCycle);% of Ca curve
Target=P.Patch.TauCa+14./P.Patch.YCaD; % target dutycyle
b=0.5; % feedback factor
P.Patch.FacYCa=min(1.1,max(0.5,P.Patch.FacYCa)).*(DutyCycle./Target).^b;
end
