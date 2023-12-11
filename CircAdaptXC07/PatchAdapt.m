function PatchAdapt(StrPatch,AdaptType)
% function PatchAdapt(StrPatch,AdaptType);
% Adapts geometry and stiffness of Patch to mechanical load
% StrPatch= array of Patch names, e.g. {'Lv1','Sv1'}
% AdaptType= {'WallVolume','WallArea','EcmStress'} = type of adaptation
% Theo Arts, Maastricht University, Nov 25, 2021

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
Aw0  = P.Wall.Aw0;
Vw0  = P.Patch.VWall;
ApRef= P.Patch.ApRef;
SfP0 = P.Patch.SfP;
FacT = ones(1,P.Patch.n);

% insertion ========
% The applied adaptation rules are not sufficient to guarantee physiologic
% size ratio of ventricular septum and left and right atrium. 
% Therefore, size of septum and both atria are stabilized by additional
% adaptation precautions

% Fixation of Septal/LvFree wall area by adjustment of SfEcmBeT
SdL0=0.6; % determines target of the ratio AwRef septal/LvFW
iW = Get('Wall','Index',{'Lv','Sv'});
iPL= iP0(iW(1))+(1:nP0(iW(1)))-1; % indices LvFW patches
iPS= iP0(iW(2))+(1:nP0(iW(2)))-1; % indices Septal patches
ALS= mean(Aw0(:,iW)); % mean wall area LvFW and Septum
SdL= ALS(2)/ALS(1); % Ratio ASep/ALvFW
AL =SdL0/SdL;
AS =1/AL;
FacT(iPL)=FacT(iPL)*AL; % FacT= multiplication factor for SfEcmBeT
FacT(iPS)=FacT(iPS)*AS;

% Fixation of size La and Ra by adjustment of SfEcmBeT
adv=0.24; %target of the ratio AwRef atrial/total ventricular
iW = Get('Wall','Index',{'La','Ra','Lv','Sv','Rv'});
iPL= iP0(iW(1))+(1:nP0(iW(1)))-1; % indices La patches
iPR= iP0(iW(2))+(1:nP0(iW(2)))-1; % indices Ra patches
A  = mean(Aw0(:,iW)); % mean wall areas
Av = sum(A(3:5)); % ventricular walls
aLa= A(1)/(adv*Av); % left  atrial adaptation factor
aRa= A(2)/(adv*Av); % right atrial adaptation factor
FacT(iPL)=FacT(iPL)*aLa;
FacT(iPR)=FacT(iPR)*aRa;
% end insertion =======

%===============
%This part has to be copied for determinimng the adaptation coefficients
DLnSns= Sens(); % column of patch-averaged sensed signals
%End copy
%===============

% Matrix dEffdSns is obtained by 'FindPatchAdaptCoef'
dEffdSns=[...
    0.3508   -1.2265    0.2985   -1.7354
    0.8239   -0.7688    0.3246   -3.2350
    0.4406    2.0946   -1.2691  -12.3534
    ];
dEff= -dEffdSns*DLnSns;

% Clipping of Fac around 1 with range +/-Clip
a= 0.1; % gain of adaptation feedback. 1.0 represents best fit matrix
ClipFac= @(x,Clip) exp(Clip*tanh(log(x)/Clip));
Clip=0.10;
FacVWall  = ClipFac(exp(a*(dEff(1,:)-log(FacT))),Clip);%adjust wall volume
FacApRef  = ClipFac(exp(a*(dEff(2,:)-log(FacT))),Clip);%adjust wall areas
FacSfP    = ClipFac(exp(a*dEff(3,:)),Clip); %adjust passive stress

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

function DLnS=Sens
global P
% Adaptation targets
SfXbMaxT   = P.Patch.Adapt.SfXbMaxT ; % Adaptation targets
SfXT       = P.Patch.Adapt.SfXT     ;
SfEcmMaxT  = P.Patch.Adapt.SfEcmMaxT;
hIsoT      = P.Patch.Adapt.hIsoT;
LnSnsTarget= log([SfXbMaxT; SfXT; SfEcmMaxT; hIsoT]);

% Row vector of sensed signals, used for adaptation
hi   = P.Patch.Lsi/2;
SfEcm= P.Patch.SfEcm;
SfXb = P.Patch.Sf-SfEcm;
SfX  = mean(SfEcm.*SfXb)./mean(SfEcm+SfXb);
Xb   = P.Patch.Xb;
SfEcmMax= Max(SfEcm);
SfXbMax = Max(SfXb);

XbThr2  = Max(Xb)/2;
hIso    = LookUp(Xb,hi,XbThr2);
LnS     = log([SfXbMax; SfX; SfEcmMax; hIso]);
DLnS    = LnS-LnSnsTarget;
end

function MaxF= Max(F)
% function Maxf= Max(F)
% Determines maximum value of signals, not exactly, but insensitive to
% spikes and less sensitive to HF components in the signals
% F   : matrix, represents sampled signals in columns
% MaxF: row vector of (approximated) maximum values per column
%
p1= 0.1; % parameter= fractional threshold around max thr=max(y)*(1-p1)
y =-sort(-F); % descending order
nx= size(y,1); % number of x samples (column)
xC= (1:nx)'-0.5; % column of x values
ym= y(round(size(y,1)/2),:); % y-value reference far from peak
Lx= ( y > y(1,:)-p1*(y(1,:)-ym) ); % near peak region in matrix

N   = sum(Lx); % number of near peak samples per column
x   = repmat(xC,[1,size(y,2)]); % x-value matrix
xLx2= x.^2.*Lx; % near peak x-values
yLx = y.*Lx; % near peak y-values
X2  = sum(xLx2)./N;     % row vectors
X4  = sum(xLx2.^2)./N;  % "
Y   = sum(yLx)./N;      % "
X2Y = sum(xLx2.*yLx)./N;% "
MaxF=(X4.*Y-X2.*X2Y)./(X4-X2.^2);% row vector of maxima
end

function g0=LookUp(f,g,f0)
% Finds rising slope crossings of function f with value level f0
% At these crossings, the value g0 of function g is looked up
% f,g: matrices, column= function samples, row= different functions
% Both matrices have the same size
% f0=row vector of f function levels at crossings
p1   = 3.0; % Normalization, Increase p1 -> narrow weighing fu around zero
Df   = f-f0; % function crosses zero at lookup level
StdDf= std(Df); % amplitude normalization reference
Df   = p1*Df./StdDf; % normalized Df
Rise = (conv2(Df([end,1:end],:),[1;-1],'valid')>0); % boolean= Rise of Df
W    = Rise.*cos(min(abs(Df),pi/2)).^2; % weighing function around lookup
g0   = sum(g.*W)./sum(W); % resulting value of function g
end
