function ArtVenAdapt(StrAv,AdaptType)
% function ArtVenAdapt(StrAv,AdaptType);
% Adaptation of Diameter and Wall thickness of Art and Ven to
% pressure and flow.
% StrAv= cell array of ArtVen names, e.g. {'Br','Pu'}
% AdaptType= {'Diameter', 'WallVolume'} indicates type of adaptation
% Theo Arts, Maastricht University, Jun 1, 2018

global P

ArtVen  = P.ArtVen;
Adapt   = ArtVen.Adapt;

% Determine indices of the adapting ArtVen's
iAv= Get('ArtVen','Index',StrAv);
% Determine type(s) of adaptation
Check=@(Str) sum(strcmp(Str,AdaptType))>0;
if Check('All') % if all types of adaptation
    AdaptDiameter= 1; AdaptWallVolume= 1;
else
    AdaptDiameter   = Check('Diameter')  ; % find Diameter adaptation type
    AdaptWallVolume = Check('WallVolume'); % find Wall adaptation type
end

% Reading P-structure information
% Parameters, ordered in matrix [2,nAv]
ASmall  = ArtVen.ASmall; % for smaller cross-section adaptation by shear
Feedback= ArtVen.AdaptFeedback;
AWall   = ArtVen.AWall(:,iAv); % wall cross-section
A0      = ArtVen.A0(:,iAv); % lumen cross-section working point
p0      = ArtVen.p0(:,iAv); % pTrans working point
% signals depending on time in matrix nt*(2*nAv)
A       = [ArtVen.AAr,ArtVen.AVe]; % Vessel cross-section
q       = ArtVen.q(:,[iAv,iAv]); % mean arterial and venous flow
p       = [ArtVen.pTransAr,ArtVen.pTransVe]; % vessel pTrans
k       = ArtVen.k(:,iAv); % stiffness coefficient
Z       = [ArtVen.ZAr,ArtVen.ZVe]; % Wave impedance
vImpact = ArtVen.vImpact; % Max impact velocity of whole body
% Target values
vMeanT  = Adapt.vFlowMean(:,iAv) ./ (1+sqrt(ASmall./A0));
% Mean target flow velocity, corrected for small vessel effect
WallStressT= Adapt.WallStress(:,iAv); % Wall stress

% adaptation target values
nAv  = numel(iAv);
szT  = [nAv,2]; % size of transposed matrix
AMax = reshape(Max(A)   ,szT)'; % max cross-section Ar and Ve connections
AMean= reshape(mean(A)  ,szT)'; % mean A
ZAMax= reshape(Max(A.*Z),szT)'; % max A*Z
pMax = reshape(Max(p)   ,szT)'+ ZAMax*vImpact; % max pressure + impact
pMean= reshape(mean(p)  ,szT)'; % mean pressure
WallStress= 3 * pMax .* (0.5+AMax./AWall); % max vessel wall stress
vMean= reshape(mean(abs(q)),szT)' ./AMean; % mean flow velocity

FacvFlow     = vMean./vMeanT; % flow velocity correction factor
FacWallStress= WallStress./WallStressT; % wall stress correction factor
Facp0        = pMean./p0; % mean pressure/target value

if AdaptDiameter
    A0= A0 .* FacvFlow.^Feedback;% Feedback 20% for optimum control
    p0= p0.*Facp0.^(3*Feedback./k);% working pressure increment
end

if AdaptWallVolume % adapts wall thickness to pressure
    AWall= AWall .* FacWallStress.^Feedback;
end
P.ArtVen.A0(:,iAv)   = A0   ; % working cross-section
P.ArtVen.p0(:,iAv)   = p0   ; % working pressure
P.ArtVen.AWall(:,iAv)= AWall; % wall cross-section

% Display data on adaptation step
% aux=[P.ArtVen.Name;repmat({' '},[1,P.ArtVen.n])];
% disp(['Deviation of ',aux{:},':'])
% disp(['vFlow     : ',num2str(round(1e3*log(FacvFlow(:)'     )),'%+6.3d')]);
% disp(['WallStress: ',num2str(round(1e3*log(FacWallStress(:)')),'%+6.3d')]);
disp(['ArtVen adaptation [vFlow WallStress]     x1000: ',...
    num2str(1e3*[std(log(FacvFlow(:))),...
    std(log(FacWallStress(:)))],'%7.0f') ]);
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
a1 = max(gds)==gds; % find steepest part in sorted f
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