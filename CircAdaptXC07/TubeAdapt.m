function TubeAdapt(StrTube,AdaptType)
% function TubeAdapt(StrTube,AdaptType);
% Adaptation of Diameter and Wall thickness of Tube to
% pressure and flow.
% StrTube= array of Tube names, e.g. {'AoBr','BrCe'}
% AdaptType= {'Diameter', 'WallVolume'} indicates type of adaptation
% Theo Arts, Maastricht University, April 2, 2014

global P

Tube    = P.Tube;
Adapt   = Tube.Adapt;
if Tube.n==0; return; end; % if no Tube exists

% Determine Tube indices
iTube= Get('Tube','Index',StrTube); %get related tube indices
% Determine type(s) of adaptation
Check=@(Str) sum(strcmp(Str,AdaptType))>0;
if Check('All') % if all types of adaptation
    AdaptDiameter= 1; AdaptWallVolume= 1;
else
    AdaptDiameter   = Check('Diameter')  ; % find Diameter adaptation type
    AdaptWallVolume = Check('WallVolume'); % find Wall adaptation type
end

% reading P-structure information
ASmall  = Tube.ASmall;
Feedback= Tube.AdaptFeedback; % feedback factor
AWall   = Tube.AWall(iTube); % wall cross-section
A0      = Tube.A0(iTube); % reference cross-section
p0      = Tube.p0(iTube); % reference transmuaral pressure
% signals depending on time
A       = Tube.A(:,iTube); % cross-section
q       = (Tube.qProx(:,iTube)+Tube.qDist(:,iTube))/2; % mean flow
pTrans  = Tube.pTrans(:,iTube); % transmural pressure
k       = Tube.k(iTube); % stiffness coefficient
Z       = (Tube.ZL(:,iTube)+Tube.ZR(:,iTube))/2; % mean wave impedance
vImpact = Tube.vImpact; % assumed body impact velocity
% targets of adaptation
vMeanT = Adapt.vFlowMean(iTube)./(1+sqrt(ASmall./A0));% target velocity
WallStressT= Adapt.WallStress(:,iTube); % target wall stress

AMax   = Max(A);
AMean  = mean(A);
ZAMax  = Max(A.*Z);
pMax   = Max(pTrans) + ZAMax*vImpact; % max. pressure
pMean  = mean(pTrans); % effective pressure load
WallStress= 3*pMax.*(0.5+AMax./AWall); % maximum wall stress
vMean  = max(mean(max(q,0)),mean(max(-q,0)))./AMean; %changed S620

% Correction of wall shear rate for small blood vessels
FacvFlow     = vMean./vMeanT; % mean velocity/target value
FacWallStress= WallStress./WallStressT; % Wall stress/target
Facp0        = pMean./p0; % mean pressure/target value

%=== Carrying out adaptation
FacV=ones(size(A0));
if AdaptDiameter; % adapts diameter to flow
    FacV = FacvFlow.^Feedback; % volume stabilization during adaptation
    A0   = A0 .* FacV; % set working cross-section
    p0   = p0 .* Facp0.^(3*Feedback./k); % set working pressure
end

if AdaptWallVolume % adapts wall thickness to pressure
    AWall= AWall .* FacWallStress.^Feedback;
end
%---
Tube.V(end,iTube)= Tube.V(end,iTube).*FacV; 
Tube.A0(iTube)   = A0; % working cross-section
Tube.p0(iTube)   = p0; % working pressure
Tube.AWall(iTube)= AWall; % wall cross-section

% % Display data on adaptation step
% aux=[Tube.Name;repmat({' '},[1,Tube.n])];
% disp(['Deviation of ',aux{:},':'])
% disp(['vFlow     : ', num2str(round(1000*log(FacvFlow    )),'%+6.3d')]);
% disp(['WallStress: ', num2str(round(1000*log(FacWallStress)),'%+6.3d')]);
disp(['Tube   adaptation [vFlow WallStress]     x1000: ',...
    num2str(1e3*[std(log(FacvFlow(:))),...
    std(log(FacWallStress(:)))],'%7.0f') ]);

P.Tube=Tube;

end

function v=Max0(M)
nt=size(M,1);
M=-sort(-M);
it=1+floor(sqrt(nt));
v=M(it,:);
end

function maxF= Max(F)
% function maxf= Max(F)
% Determines maximum value of signals, not exactly, but insensitive to
% spikes and less sensitive to HF components in the signals
% F   : matrix, represents sampled signals in columns
% MaxF: row vector of (approximated) maximum values per column
% 
s=sort(F);
ds=conv2(s,[1;0;-1] ,'valid');%column derivative sorted f
ys=conv2(s,[1;2;1]/4,'valid');%value sorted f
nds=size(ds,1); % number of samples after convolution
g=(1:nds)'; g=g.*flipud(g); % weighing function =(x*(endx-x)), boundaries=0
gds= g.*ds; % s multiplied by weighing function
a1 = (max(gds)==gds); % find steepest part in sorted f
a2 = cumsum(a1); % find multiple maxima
thr= ys(a1 & a2==1)'; % select the first maximum only, use as threshold
a3 = ys>thr; % select values above threshold
n  = sum(a3); % number of super-threshold values per column
y1 = sum(ys.*a3); % sum of those values
y2 = sum(ys.^2.*a3); % sum of squares
mn = y1./n; % mean of super-threshold values
sd = sqrt(y2./n-mn.^2); % sd of super-threshold values
maxF= mn+sqrt(1.25)*sd; % estimated maximum value
end