function AdaptPatch(StrPatch,AdaptType)
% function PatchAdapt(StrPatch,AdaptType);
% Adapts geometry and stiffness of Patch to mechanical load
% StrPatch= array of Patch names, e.g. {'Lv1','Sv1'}
% AdaptType= {'WallVolume','WallArea','EcmStress'} = type of adaptation
% Theo Arts, Maastricht University, March 12, 2024

global P

iPatch= Get('Patch','Index',StrPatch); %get related Patch indices
Check=@(Str) sum(strcmp(Str,AdaptType))>0;
if Check('All')
    AdaptWallVolume=1; AdaptWallArea=1; AdaptEcmStress=1;
else
    AdaptWallVolume=Check('WallVolume');
    AdaptWallArea  =Check('WallArea'); %-> Myom area
    AdaptEcmStress =Check('EcmStress'); %-> EcmArea
end

iP0  = P.Wall.iPatch; % first index of patches connected to the wall
nP0  = P.Wall.nPatch; % number of patches (multi-patch) per wall
Vw0  = P.Patch.VWall;
ApRef= P.Patch.ApRef;
SfP0 = P.Patch.SfP;

% Find L, S and R patches
iW = Get('Wall','Index',{'Lv','Sv','Rv'});
iPL= iP0(iW(1))+(1:nP0(iW(1)))-1; % indices LvFW patches
iPS= iP0(iW(2))+(1:nP0(iW(2)))-1; % indices Septal patches
iPR= iP0(iW(3))+(1:nP0(iW(3)))-1; % indices RvFW patches
iV = [iPL,iPS,iPR]; % ventriculaire patches

iV2W=[iPL*0+1,iPS*0+2,iPR*0+3];

% Patch boundaries adapt using geometric wall properties
At= P.Wall.Aw(:,iW); % wall area
Ct= P.Wall.Cm(:,iW); % wall curvature
Tt= P.Wall.T(:,iW); % wall tension
ht= P.Wall.VWall(iW)./At; % wall tension
zt = At.*Ct.^2/(4*pi); % wall fraction on spherical surface
T=mean(Tt);
z=mean(zt.*Tt)./T;
C=mean(Ct.*Tt)./T;
h=mean(ht.*Tt)./T;
s = 2*z-1; % vertical coordinate of wall vector in junction
c = sign(C).*sqrt(max(0,1-s.^2)); % horizontal coordinate

% Change of wall vector orientation to obtain physiological shape of the
% TriSeg geometry
aa= 0.4; % ++0.4++ correction parameter for size of Lv, Sv and Rv wall segments
c = aa*c;
k = sqrt(s.^2+c.^2);
s = s./k;
c = c./k; % new wall vector set

% Center of gravity of wall vectors
cMn= mean(c);
sMn= mean(s);
Aux= h.^4;
wh = Aux/sum(Aux); % weight per wall segment by wall tension

q  = wh.*(cMn*c+sMn*s); % wall growth factor
% q  = q-h*sum(q); % set average wall volume change to zero
qW=2*q(iV2W);

%===============
%This part has to be copied to get the sensed adaptation deviations
DLnSns= Sens(); % column of patch-averaged sensed signals
%End copy
%===============
DLnSns(1,iV)=DLnSns(1,iV)+1.0*qW; % geometric adjustment for TriSeg
DLnSns(2,iV)=DLnSns(2,iV)-0.0*qW; % geometric adjustment for TriSeg
DLnSns(3,iV)=DLnSns(3,iV);

dEffdSns=[...
   -1.2436   -0.1872    2.2845
    0.2741   -0.0448   -1.3219
    0.1824    0.7808  -16.1288
    ];
dEff= -dEffdSns*DLnSns;

% Clipping of Fac around 1 with range +/-Clip
a= 0.2; % gain of adaptation feedback. 1.0 represents best fit matrix
ClipFac= @(x,Clip) exp(Clip*tanh(log(x)/Clip));
Clip=0.10;
FacVWall  = ClipFac(exp(a*(dEff(1,:))),Clip);%adjust wall volume
FacApRef  = ClipFac(exp(a*(dEff(2,:))),Clip);%adjust wall areas
FacSfP    = ClipFac(exp(a* dEff(3,:)),Clip); %adjust passive stress

%=== Carrying out adaptation

if AdaptWallVolume
    Vw0= Vw0 .* FacVWall;
end

if AdaptWallArea
    ApRef= ApRef .* FacApRef;
end

if AdaptEcmStress
   SfP0= SfP0 .* FacSfP  ;
end

P.Patch.VWall(iPatch)= Vw0(iPatch);
P.Patch.ApRef(iPatch)= ApRef(iPatch);
P.Patch.SfP(iPatch)  = SfP0(iPatch);

% Display adaptation process
disp(['Patch adaptation [Volume Area Stiffness] x1000: ',num2str(1e3*[...
    std(log(FacVWall(:))),...
    std(log(FacApRef(:)))...
    std(log(FacSfP(:)))...
    ],'%7.0f')]);

end

% Function Sens reads target of sensed signals and the deviation from the
% target values

function DLnS=Sens()
global P
Patch=P.Patch;
% Adaptation targets
XbMaxT   = Patch.Adapt.XbMaxT; % Adaptation targets
SfXT     = Patch.Adapt.SfXT  ;
hMeanT   = Patch.Adapt.hMeanT;

LnSnsTarget  = log([XbMaxT; SfXT; hMeanT]);

% Row vector of sensed signals, used for adaptation
hi   = Patch.Lsi/2;
SfEcm= Patch.SfEcm;
Xb   = Patch.Xb;

XbMax = Mx(Xb);
SfX   = mean(SfEcm.*Xb)./mean(Xb);
hMean = mean(Xb.*hi)./mean(Xb);

LnS = log([XbMax; SfX; hMean]);
DLnS= LnS-LnSnsTarget;
end

function mx=Mx(y) % determines maximum value
y=sort(y,'descend');
s=2*y(1,:)-y(2,:)-y(3,:);
v=y(2,:)-y(3,:);
mx=y(1,:)+0.125*v.^2./s;
end
