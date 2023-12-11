function PNew
% function PNew
%
% Sets backbone of structure P
% Attribute names to ArtVen, TriSeg, Chambers, Valves, Tubes, 
% Walls and Patches
% Connections between elements are defined by strings stored in
% structure field P.Map
%
% Structure is built on information in P.Map
% ArtVen : artery-peripheral resistance-vein of organ or body part
% Chamber: cavity enclosed by a myocardial wall (atria)
% TriSeg : combination of two cavities with three walls
% Bag    : passive elastic bag around part of the circulation (pericardium)
% Node   : named connection point
% Wall   : muscular wall can contract, composed of 1 or more patches
% Patch  : contracting part of a wall, having specific mechanical properties
% Valve  : valve with inertia connects proximal to distal node, may leak
% Tube   : elastic wave-guiding tube connects proximal to distal node
%
% Theo Arts, Maastricht University, June 8, 2018

global P;
P=CreateP; % sets the tree structure with all necessary fields

%======= DEFINITION OF STRUCTURE BY ELEMENTS AND CONNECTIONS ==========
Map.ArtVenName = {'Ca','Br','Ce','Fe','Pu'}; % carotid, brachial, celiac,
% femoral and pulmonary ArtVen
Map.ChamberName= {'La','Ra'}; %atria
Map.TriSegName = {'v'}; %ventricular unit

% Valve and Tube connections
Map.ValveNodes=[... % proximal and distal node of valves
    {'Vc'  ,'Ra'  };...
    {'Ra'  ,'Rv'  };...
    {'Rv'  ,'PuAr'};...
    {'PuVe','La'  };...
    {'La'  ,'Lv'  };...
    {'Lv'  ,'Ao'  }];

Map.TubeNodes=[... % proximal and distal node of elastic tubes
    {'Ao'  ,'BrAr'};...
    {'BrAr','CaAr'};...DepPath
    {'BrAr','CeAr'};...
    {'CeAr','FeAr'};...
    {'Vc'  ,'BrVe'};...
    {'BrVe','CaVe'};...
    {'BrVe','CeVe'};...
    {'CeVe','FeVe'}];

% Bags pressurize parts of the circulation
Map.Bag.Name   ={'Peri'     ,'Thorax'}; % pericardium, thorax
Map.Bag.Chamber={{'La','Ra'},{}      }; % enclosed chambers
Map.Bag.TriSeg ={{'v'}      ,{}      }; % enclosed TriSeg
Map.Bag.ArtVen ={{}         ,{'Pu'}  }; % enclosed ArtVen
Map.Bag.Tube   ={{}         ,{'AoBrAr','BrArCeAr','VcBrVe','BrVeCeVe'}};
Map.Bag.Bag    ={{}         ,{'Peri'}}; % pericardium inside thorax

Map.Baro ={'BrAr'}; % pressure control node
Map.Pace = 'Ra1'; % leading pacemaker patch
Map.MultiPatch= {}; % defines Wall's split in Patch's
Map.AvDepPath = 'Ra1Sv1'; %Av-depolarization pathway
DepPath={...
    'Ra1','Ra1';...
    'Ra1','La1';...
    'Ra1','Sv1';...
    'Sv1','Lv1';...
    'Sv1','Rv1';...
    };
Map.DepPath   = DepPath;

% Transfer of structure information to P-structure
P.Map = Map;

Indexation; % mutual element relations expressed by indices

%==========================================
%============== Data filling ==============
%==========================================

% Species specific general information
P.General.q0           = 85e-6; % mean systemic flow
P.General.p0           = 12200; % mean systemic pressure
P.General.tCycle       = 0.85 ; % cycle time
P.General.Dt           = 0.002;

P.General.tCycleRef    = 0.85 ; % reference cycle time
P.General.RhoB         = 1050 ; % blood density
P.General.EtaB         = 0.004; % blood viscosity
P.General.TauAv        = 0.12 ; % AV depolarization delay
P.General.TauAvRef     = 0.13 ; % AV depolarization delay
P.General.TauAvExc     = 0.10 ; % AV depolarization delay
P.General.AdaptFunction= 'Adapt';
P.General.SaturationControl= 0  ;
P.General.AdaptFeedback    = 0.5; % feedback strength
P.General.nStoredBeats     = 3; % maximum number of stored beats in Ps
P.General.ErrStop          = 1.0; % Error to end the simulation

% Crude estimates of maximum pressures at exercise
pLv=19000   ; % peak Lv pressure
pRv= 7000   ; % peak Rv pressure
pLa= 2000   ; % peak pLa
pRa= 1100   ; % peak pRa
pPu= pLa+300; % mean pulmonary artery pressure
p0 = P.General.p0;

% General settings
P.General.FacpControl= 1    ;
P.General.tEnd       = 1.0; % end of simulation, finish with full beat

% Parameter settings per element type
AdaptationParameters
SarcomereProperties
MakeHeart(pLv,pRv,pLa,pRa)
ArtVenParameters(p0,pPu,pLa,pRa)
TubeParameters
BagParameters
ValveParameters
DepPathParameters

P.t=0;
TimingInit % Determines depolarization sequence
P2SVar;

save P P

end

% ================= Auxilary functions ====================
function AdaptationParameters

global P;

% ArtVen.Adapt and Tube.Adapt
Put({'ArtVen','Adapt'},'WallStress','All',[1;1]*500e3);
Put({'ArtVen','Adapt'},'vFlowMean' ,'All',[1;1]*0.46 ); % excercise level
Put({'Tube'  ,'Adapt'},'WallStress','All',500e3);
Put({'Tube'  ,'Adapt'},'vFlowMean' ,'All',0.50 );

% Patch/Sarcomere
PatchA={}; % atrial patches
for iP=1:P.Patch.n
    if P.Patch.Name{iP}(2)=='a'
        PatchA=[PatchA,P.Patch.Name(iP)]; % find atrial patches
    end
end
% Ventricle=default
Put({'Patch','Adapt'},'SfXbMaxT' ,'All' , 84000); % max active Sf
Put({'Patch','Adapt'},'SfXbMaxT' ,PatchA, 36000); % max active Sf atria
Put({'Patch','Adapt'},'SfXT'     ,'All' ,   840); % max Xb*Ecm Sf
Put({'Patch','Adapt'},'SfXT'     ,PatchA,  6400); % max Xb*Ecm Sf atria
Put({'Patch','Adapt'},'SfEcmMaxT','All' ,  2400); % ECM Sf at BE
Put({'Patch','Adapt'},'SfEcmMaxT',PatchA, 25000); % ECM Sf at BE atria
Put({'Patch','Adapt'},'hIsoT'    ,'All' ,  1.13); % ~(LSarc/2) at B.Ej.
end

function SarcomereProperties
global P;
PatchA={};
for iP=1:P.Patch.n
    if P.Patch.Name{iP}(2)=='a'
        PatchA=[PatchA,P.Patch.Name(iP)];
    end
end

% Ventricular: default
Put('Patch','Lsi'             ,'All' , 2.1  );
Put('Patch','Xb'              ,'All' , 0.1);
Put('Patch','TauRefrac'       ,'All' , 0.25 ); % refractory period

Put('Patch','kP','All',15.0); % passive stiffness
Put('Patch','SfP',{'La1','Ra1','Lv1','Sv1','Rv1'},[8 4 5 5 5]*1e2);
    % passive stress at reference

% Sarcomere structure, h.. refers to half sarc length
P.Patch.hAct = 1.050; % [um] thin filament length
P.Patch.hMyo = 0.825; % [um] myosin length including bare zone
P.Patch.hBare= 0.050; % [um] bare zone length
P.Patch.fD   = 0.33 ; % Fraction D-zone
P.Patch.LnbC =-1.30*[1 1 1 1 1]; % Fraction Ca sensitivity of C-zone
P.Patch.gTit = 1.82*[1 1 1 1 1]; % Titin-Xb dissociation by strain

% Ca-Xb properties
P.Patch.hRef = 0.0144*[1 1 1 1 1]; % [um] ref length myosin-Xb tension
P.Patch.hSe  = 0.02; % [um] half sarc length series elasticity
P.Patch.vMx  = [12*[1 1], 8*[1 1 1]]; % [um/s]
P.Patch.SfA  = 4500 *[1 1 1 1 1];% MechChem: scaling to number of Xb's

% Ca Pulse
P.Patch.TauCa  =[0.03*[1 1],0.18*[1 1 1]];
Put('Patch','CaS'  ,'All',0.25); %Systolic Ca-injection
Put('Patch','CaD'  ,'All',0.11); %Diastolic Ca-removal
Put('Patch','YCaS' ,'All', 4500*[1 1 1 1 1]);
Put('Patch','YCaD' ,'All', 1800*[1 1 1 1 1]);
P.Patch.FacYCa= [1 1 1 1 1]; % controls systolic duration by YCaS and YCaD
% rate dependency of Ca-pulse parameters
P.Patch.aTauCa=  0.6; % duration systolic Ca-pulse
P.Patch.aCaS  =  0.5; % systolic [Ca]-source
P.Patch.aCaD  =  0.2; % diastolic [Ca]-source
P.Patch.aYCaS = -1.8; % systolic pulse Ca conductivity
P.Patch.aYCaD =  0.0; % diastolic background Ca conductivity
P.Patch.avMx  = -1.0; % vMx increases with heart rate

end
%=============================

function MakeHeart(pLv,pRv,pLa,pRa)
%pLv,pRv = peak pressures Lv, Rv
%pLa,pRa= mean pressures La, Ra

global P;

VStroke= P.General.q0*P.General.tCycle; %Stroke volume
x=log(pLv/pRv);
z=1./(1+x.^2);
% Midwall areas
SumAw = 6.0*(3*VStroke)^(2/3); % ==L+S+R midwall wall area
LRdSAw= 5.0-z; % ==(L+R)/S
LdRAw = exp(-0.4*tanh(x)); % ==L/R
Aw=SumAw/(1+LRdSAw)*[LRdSAw/(1+1/LdRAw),1,LRdSAw/(1+LdRAw)];
% Solution LSR midwall areas

%Wall volumes
SumVw = 1.5*VStroke*(pLv+pRv)/15000; % ==L+R+S
LRdSVw= 1/(0.33-0.13*z); % ==(L+R)/S
LdRVw = exp(0.7*x); % == L/R
Vw=SumVw/(1+LRdSVw)*[LRdSVw/(1+1/LdRVw),1,LRdSVw/(1+LdRVw)];
% Solution LSR wall volumes

AAw=0.12*SumAw*[1,1]; % atrial midwall areas
AVw=0.35*[pLa,pRa].*SumVw/(pLv+pRv); % atrial wall volumes

PatchA={'La1','Ra1'};
PatchV={'Lv1','Sv1','Rv1'};
Put('Patch','ApRef',PatchV,Aw);
Put('Patch','VWall',PatchV,Vw);
Put('Patch','ApRef',PatchA,AAw);
Put('Patch','VWall',PatchA,AVw);

% Init values of heart cavity volume (state variables V)
P.Chamber.V= VStroke*[0.5,0.5];
P.TriSeg.VL= VStroke;
P.TriSeg.VR= VStroke;
P.TriSeg.V = 0.5*VStroke; % septal displacement volume
P.TriSeg.Y = 0.8*VStroke^(1/3); % septal boundary radius
end

function ArtVenParameters(p0,pPu,pLa,pRa)
global P
Sy  = {'Ca','Br','Ce','Fe'}; % systemic ArtVen
Pu  = 'Pu'; % Pulmonary ArtVen
q0  = P.General.q0; % total systemic flow
vFlowMean= P.ArtVen.Adapt.vFlowMean;

%   peripheral flow distribution over ArtVen vessel beds
Put('ArtVen','qRef' ,Sy, [0.15,0.12,0.57,0.16]*q0); % ArtVen flow targets
Put('ArtVen','qRef' ,Pu, q0); % Pu-flow=sytemic flow
kExc= [0 1.8 -0.5 1.8 1]; % ArtVen flow dependency with exercise
P.ArtVen.kExc= kExc; % ArtVen flow dependency with exercise
Put('ArtVen','p0' ,Sy,[p0 ; pRa]); % Sy- Art and Ven pressures
Put('ArtVen','p0' ,Pu,[pPu; pLa]); % Pu- Art en Ven pressure
p0Av=P.ArtVen.p0;
qRef=P.ArtVen.qRef;
Put('ArtVen','k','All',[[14,12,12,16,9];[14,12,12,16,9]]);
P.ArtVen.q0Av= qRef;
Put('ArtVen','Len','All',6.0*qRef.^(1/3)); % size ArtVen-bed
Put('ArtVen','Len',{'Br','Fe'},[0.3,0.5]); % Arm and leg are long
Len      = P.ArtVen.Len;

% Parameters for adaptation
P.ArtVen.ASmall       = 0.7e-6; % smaller blood vessels adapt shear
P.ArtVen.vImpact      = 3.0; % maximum impact velocity
P.ArtVen.AdaptFeedback= 0.2; % gain of adaptation feedback
qExc= [1;1]*(3.^kExc.*qRef);
A0  = qExc./vFlowMean;


% Dependent parameter values, init before adaptation
P.ArtVen.A0   = A0; % cross-sect.
P.ArtVen.AWall= diag([0.20;0.10])*A0; % A+V wall cross-section
P.ArtVen.p0Av = p0Av(1,:)-p0Av(2,:); % A+V reference pressure

% Initialization cavity volumes of Artven (= state variables)
V            = A0 .* [Len;Len];
P.ArtVen.VAr = V(1,:);
P.ArtVen.VVe = V(2,:);
end

function TubeParameters
global P
% Flow distribution calculated from ArtVen flows and shunt flows
FlowDistribution;  % mean flow distribution in ArtVen, Tubes, Valves
p0   = P.General.p0; % arterial pressure
Tube = P.Tube;
Name = Tube.Name;
Adapt= Tube.Adapt;
q    = Tube.q;
Len  = Tube.Len;
Row  = zeros(1 ,Tube.n);

vFlowMean = Adapt.vFlowMean;
WallStress= Adapt.WallStress;
vImpact   = 3.0;


% Sequence TubeNames:
aux = regexpi(Name,'Ve');
iArt=  cellfun('isempty',aux); % Arterial tubes
iVen= ~cellfun('isempty',aux); % Venous tubes
Aux = Get('ArtVen','p0','Br'); % get brachial Art and Ven pressure
pRa = Aux(2); % pRa= brachial venous pressure
A0  = 3*abs(q./vFlowMean); % Tube cross-section
Len(iArt)= [4,15,20,15]/100;
Len(iVen)= Len(iArt);
k(iArt)  = [10,10,10,14];
k(iVen)  = k(iArt);

% Parameter related to adaptation
Tube.ASmall       = 0.7e-6; % Diam=0.9 mm Ref Arts 2012
Tube.vImpact      = vImpact; % max impact velocity
Tube.AdaptFeedback= 0.3; % feedback rate of adaptation

% pressure distribution in tubes, representing vessel segments
Tube.p0(iArt) = p0;
Tube.p0(iVen) = pRa;
p0Tb =Tube.p0;

% additional tube properties
Tube.k        = k;
Tube.Len      = Len;
%Dependent initializations
Tube.A0   = A0;
Tube.AWall= A0.*(12*p0./WallStress+vImpact*0.02);
Tube.V    = Len.*A0;

% Memory allocation for Wave delay lines, needed for Tube function
Tube.pL = p0Tb;% Left pressure wave, non-delayed, circular storage
Tube.pR = p0Tb;
u0=sqrt(p0Tb*P.General.q0);
Tube.uP = u0;% Delayed and attenuated pressure signal
Tube.uD = u0;
Tube.TauL= Row+P.General.Dt;% delay
Tube.TauR= Row+P.General.Dt;

P.Tube= Tube;
end

function BagParameters
global P;
% Pericardium and  Thorax Bags
P.Bag.k     = [10,10]; % volume stiffness
P.Bag.pAdapt= [500,50]; % transmural bag pressure pressure
%Estmate heart volume
VWall = sum(P.Patch.VWall);
VCav  = sum([P.Chamber.V,P.TriSeg.VL,P.TriSeg.VR]);
VHeart= VCav+VWall;
P.Bag.VRef= [1.4*VHeart,2.7*VHeart]; % Bag target volumes
end

function ValveParameters

A=Get('Tube','A0','AoBrAr'); % Aortic cross-section

Put('Valve','q'    ,'All',0.0);
Put('Valve','AOpen','All',A  ); % AOpen in arteries == aortic cross-section
Put('Valve','ALeak','All',A*1e-6); % ALeak is small but >0
Put('Valve','Len'  ,'All',sqrt(A)); % default valve-channel length
% Mitral and Tricuspid valve are larger
Put('Valve','AOpen',{'RaRv','LaLv'},...
    1.5* Get('Valve','AOpen',{'RaRv','LaLv'}) );
% vene-atrial orifices are always open
Put('Valve','ALeak',{'VcRa','PuVeLa'},...
    Get('Valve','AOpen',{'VcRa','PuVeLa'}) );

% Wall: ApDead (non-contractile area)=valve orifice area
ALv=sum(Get('Valve','AOpen',{'LvAo'  ,'LaLv'}));
ARv=sum(Get('Valve','AOpen',{'RvPuAr','RaRv'}));
ALa=sum(Get('Valve','AOpen',{'PuVeLa','LaLv'}));
ARa=sum(Get('Valve','AOpen',{'VcRa'  ,'RaRv'}));
Put('Wall','ApDead',{'Lv','Rv','La','Ra'},[ALv,ARv,ALa,ARa]);
end

function DepPathParameters
global P
Put('DepPath','Delay','All',...
    [P.General.tCycle,0.01,P.General.TauAv,0.005,0.005]);
end

function FlowDistribution
% calculates flow distribution through ArtVens, Valves and Tubes
% Flow ~=0 are assumed to be known. All other flows are unknown. If
% a sufficient number of flows is known, equations on steady state flow
% distribution are solved in a least squares sense

global P

nNode  = P.Node.n;
nValve = P.Valve.n;
nTube  = P.Tube.n;
nArtVen= P.ArtVen.n;

Ma =zeros(nNode,nArtVen);
Mt =zeros(nNode,nTube);
Mv =zeros(nNode,nValve);

for ia=1:nArtVen
    iP=P.ArtVen.iNode(ia);
    iD=iP+1;
    Ma(iP,ia)=Ma(iP,ia)-1;
    Ma(iD,ia)=Ma(iD,ia)+1;
end
for ia=1:nTube
    iP=P.Tube.iNodeProx(ia);
    iD=P.Tube.iNodeDist(ia);
    Mt(iP,ia)=Mt(iP,ia)-1;
    Mt(iD,ia)=Mt(iD,ia)+1;
end
for ia=1:nValve
    iP=P.Valve.iNodeProx(ia);
    iD=P.Valve.iNodeDist(ia);
    Mv(iP,ia)=Mv(iP,ia)-1;
    Mv(iD,ia)=Mv(iD,ia)+1;
end
M  = [Ma,Mt,Mv];
Rga= 1:nArtVen;
Rgt= nArtVen+(1:nTube);
Rgv= nArtVen+nTube+(1:nValve);
nM = size(M,2);
q  = [P.ArtVen.qRef,P.Tube.q(1,:),P.Valve.q(1,:)];
Rg0= find(q~=0); %known flows q
Rg1= setdiff(1:nM,Rg0); % unknown flows q
q0 = q(Rg0);
q1 = -pinv(M(:,Rg1))*M(:,Rg0)*q0';
q(Rg1)       = q1;
P.ArtVen.qRef= q(Rga);
P.Tube.q     = q(Rgt);
P.Valve.q    = q(Rgv);
end


function P=CreateP
% Creates empty structure P with fields
P=[];

FieldsP={
    'General'
    'ArtVen'
    'Chamber'
    'TriSeg'
    'Valve'
    'Tube'
    'Node'
    'Wall'
    'Patch'
    'Bag'
    'DepPath'
    'Map'
    'SVar'
    'SVarDot'
    't'
    'tDot'
    };

FieldsGeneral={
    'q0'
    'p0'
    'tCycle'
    'tCycleRef'
    'tEnd'
    'Dt'
    'dt'
    'BeatNr'
    'FacpControl'
    'TauAv'
    'TauAvRef'
    'TauAvExc'
    'RhoB'
    'EtaB'
    'AdaptFunction'
    'In'
    'Out'
    'AdaptFeedback'
    'nStoredBeats'
    'SaturationControl'
    'ErrStop'
    };

FieldsArtVen={
    'Name'
    'n'
    'iNode'
    'k'
    'Len'
    'A0'
    'p0'
    'AWall'
    'p0Av'
    'qRef'
    'kExc'
    'q0Av'
    'ASmall'
    'vImpact'
    'AdaptFeedback'
    'VAr'
    'VVe'
    'VArDot'
    'VVeDot'
    'q'
    'qAr'
    'qVe'
    'pTransAr'
    'pTransVe'
    'ZAr'
    'ZVe'
    'AAr'
    'AVe'
    'pExt'
    'Adapt'
    };

FieldsChamber={
    'Name'
    'n'
    'iWall'
    'iNode'
    'V'
    'VDot'
    'A'
    'Z'
    'pTrans'
    'pExt'
    };

FieldsTriSeg={
    'Name'
    'n'
    'iWall'
    'iNode'
    'VL'
    'VR'
    'VLDot'
    'VRDot'
    'V'
    'Y'
    'VDot'
    'YDot'
    'VS'
    'YS'
    'AL'
    'AR'
    'ZL'
    'ZR'
    'pTransL'
    'pTransR'
    'pExt'
    };

FieldsWall={
    'Name'
    'n'
    'nPatch'
    'iPatch'
    'VWall'
    'Aw0'
    'ApDead'
    'DADT'
    'T'
    'Cm'
    'Aw'
    'pTrans'
    'Wall2Patch'
    };

FieldsPatch={
    'Name'
    'n'
    'VWall'
    'ApRef'
    'SfP'
    'kP'
    'SfA'
    'TauRefrac'
    'iPace'
    'iDepPathProx'
    'tDep'
    'tEnd'
    'hAct'
    'hMyo'
    'hBare'
    'hRef'
    'gTit'
    'LnbC'
    'fD'
    'hSe'
    'vMx'
    'TauCa'
    'CaS'
    'CaD'
    'YCaS'
    'YCaD'
    'aTauCa'
    'avMx'
    'aCaS'
    'aCaD'
    'aYCaS'
    'aYCaD'
    'FacYCa'
    'Adapt'
    'Lsi'
    'LsiDot'
    'Xb'
    'XbDot'
    'Ef'
    'Ls'
    'SfEcm'
    'Sf'
    'DSfDEf'
    'T'
    'DADT'
    'Ap0'
    'Ap'
    'Ca'
    };

FieldsNode={
    'Name'
    'n'
    'iBaro'
    'q'
    'p'
    'Y'
    'A'
    };

FieldsBag={
    'Name'
    'n'
    'iTube'
    'iChamber'
    'iTriSeg'
    'iArtVen'
    'iBag'
    'OK'
    'VRef'
    'k'
    'pAdapt'
    'V'
    'pTrans'
    'pExt'
    'p'
    'Ch'
    'Tr'
    'Av'
    'Tb'
    'Bg'
    };

FieldsValve={
    'Name'
    'n'
    'iNodeProx'
    'iNodeDist'
    'AOpen'
    'ALeak'
    'Len'
    'q'
    'qDot'
    'AvValves'
    'AvWalls'
    'Valve2NodeProx'
    'Valve2NodeDist'
    };

FieldsDepPath={
    'Name'
    'n'
    'iPatchProx'
    'iPatchDist'
    'Delay'
    'iCycle'
    'iAv'
    };

FieldsTube={
    'Name'
    'n'
    'iNodeProx'
    'iNodeDist'
    'k'
    'Len'
    'A0'
    'p0'
    'AWall'
    'dt'
    'ASmall'
    'vImpact'
    'AdaptFeedback'
    'V'
    'VDot'
    'q'
    'A'
    'pTrans'
    'pSProx'
    'pSDist'
    'pExt'
    'qProx'
    'qDist'
    'ZL'
    'ZR'
    'pL'
    'pR'
    'uP'
    'uD'
    'TauL'
    'TauR'
    'cL'
    'cR'
    'Att'
    'Tube2NodeProx'
    'Tube2NodeDist'
    'Adapt'
    };

FieldsMap={
    'General'
    'ArtVenName'
    'ChamberName'
    'TriSegName'
    'ValveNodes'
    'TubeNodes'
    'Bag'
    'Baro'
    'Pace'
    'MultiPatch'
    'AvDepPath'
    'DepPath'
    };

Empty10=ones(1,0);
for i=1:length(FieldsP)
    fld=FieldsP{i};
    P.(fld)=[];
    FldStr=['Fields',fld];
    if exist(FldStr,'var')
        SubFldStr=eval(FldStr);
        for j=1:length(SubFldStr)
            P.(fld).(SubFldStr{j})=Empty10;
        end
    end
end

end

