function pNodeVDot
% function pNodeVDot
% Input: Source pressures and impedances of the outlets of the elements
% ArtVen, Chamber, TriSeg, Tube, and Valve
% The elements are connected to the nodes
% Output: pressures in Nodes and flows from nodes to elements
% Includes collapsibility of tubes
% Theo Arts, Maastricht University, March 7, 2024

global P

% counteracting collapse of cavities by negative pressure at small VCav
Facq    = P.General.FacpControl; % pressure control by volume change
nt      = numel(P.t);

iCh2Nd        = P.Chamber.iNode;
YCh           = 1./P.Chamber.Z  ; % conductivity of Chamber
ACh           = P.Chamber.A     ; % flow cross-section of chamber
pSCh          = P.Chamber.pTrans; % zero-flow transmural pressure
pExtCh        = P.Chamber.pExt  ; % external pressure

iTr2NdL       = P.TriSeg.iNode  ; % node index related to L cavity
YTrL          = 1./P.TriSeg.ZL  ; % conductivity of left ventricle
YTrR          = 1./P.TriSeg.ZR  ; % conductivity of right ventricle
ATrL          = P.TriSeg.AL     ; % flow cross-section of LV
ATrR          = P.TriSeg.AR     ; % flow cross-section of RV
pSTrL         = P.TriSeg.pTransL; % zero-flow left transmural pressure
pSTrR         = P.TriSeg.pTransR; % zero-flow right transmural pressure
pExtTr        = P.TriSeg.pExt   ; % external pressure

Valve2NodeProx= P.Valve.Valve2NodeProx; % linking matrix
Valve2NodeDist= P.Valve.Valve2NodeDist; % linking matrix
qValve        = P.Valve.q    ; % valve flow
AOpen         = P.Valve.AOpen; % open state cross-section
ALeak         = P.Valve.ALeak; % closed state cross-section

nAv           = P.ArtVen.n       ; % number of ArtVen elements
iAv2NdAr      = P.ArtVen.iNode   ; % indices of related nodes
LenAv         = P.ArtVen.Len     ; % reperesentative large vessel length
qAv           = P.ArtVen.q       ; % microcirculatory art-ven flow
YAr           = 1./P.ArtVen.ZAr  ; % arterial wave conductivity
YVe           = 1./P.ArtVen.ZVe  ; % venous wave conductivity
pTransAr      = P.ArtVen.pTransAr; % zero-flow prox transmural pressure
pTransVe      = P.ArtVen.pTransVe; % zero-flow dist transmural pressure
VAr           = P.ArtVen.VAr     ; % volume arterial compartment
VVe           = P.ArtVen.VVe     ; % volume venous compartment
pExtAv        = P.ArtVen.pExt    ; % external pressure

nTb           = P.Tube.n; % number of Tube elements
iTb2NdP       = P.Tube.iNodeProx    ; % proximal node indices
iTb2NdD       = P.Tube.iNodeDist    ; % distal node indices
Tube2NodeProx = P.Tube.Tube2NodeProx; % linking matrix
Tube2NodeDist = P.Tube.Tube2NodeDist; % linking matrix
YTbP          = 1./P.Tube.ZR ; % R-wave conductivity per tube
YTbD          = 1./P.Tube.ZL ; % L-wave conductivity per tube
pSTbP         = P.Tube.pSProx; % zero-flow prox transmuaral pressure
pSTbD         = P.Tube.pSDist; % zero-flow dist transmuaral pressure
ATb           = P.Tube.A     ; % length averaged cross-section (fu(t))
pExtTb        = P.Tube.pExt  ; % external pressure

nNd           = P.Node.n; % number of nodes

% Memory allocation
YNode=zeros(nt,nNd);
qNode=zeros(nt,nNd);
ANode=zeros(nt,nNd);

% Chamber
pSCh   = pSCh + pExtCh ; % absolute source pressure
qSCh   = pSCh .* YCh   ; % source flow to node
YNode(:,iCh2Nd)=YNode(:,iCh2Nd)+YCh ; % added Chamber related conductivity
qNode(:,iCh2Nd)=qNode(:,iCh2Nd)+qSCh; % added Chamber related inflow
ANode(:,iCh2Nd)=ANode(:,iCh2Nd)+ACh ; % added Chamber related cross-section

% TriSeg
iTr2NdR = iTr2NdL+1; % right node indices
iTr2Nd  = [iTr2NdL,iTr2NdR];
pSTrL   = pSTrL+pExtTr ; % absolute source pressure
pSTrR   = pSTrR+pExtTr ; % absolute source pressure
pSTr    = [pSTrL,pSTrR]; % absolute source pressure
YTr     = [YTrL,YTrR]; % conductivity
ATr     = [ATrL,ATrR]; % cross-sectional area
qSTr    = pSTr .* YTr; % source flow to node
YNode(:,iTr2Nd)=YNode(:,iTr2Nd)+YTr ; % added TriSeg related conductivity
qNode(:,iTr2Nd)=qNode(:,iTr2Nd)+qSTr; % added TriSeg related inflow
ANode(:,iTr2Nd)=ANode(:,iTr2Nd)+ATr ; % added TriSeg related cross-section

%==== Valve flow contribution to node pressure
DqValve= qValve*(-Valve2NodeProx+Valve2NodeDist); %Added flow into node
DAValve= repmat(max(AOpen,ALeak) * ...
    (Valve2NodeProx+Valve2NodeDist),[nt,1]); %Added flow area around node
qNode= qNode+DqValve;
ANode= ANode+DAValve;

%==== AV pressure ArtVen to node pressure and flow
iAv2NdVe= iAv2NdAr+1; % distal node indices
iAv2Nd  = [iAv2NdAr,iAv2NdVe]; % node indices
AAv     = [VAr,VVe]./[LenAv,LenAv]; % cross-section
YAv     = [YAr,YVe] ; % source conductivity
pSAr    = pTransAr+pExtAv-qAv./YAr; % absolute source pressure, effect qAv
pSVe    = pTransVe+pExtAv+qAv./YVe; % absolute source pressure, effect qAv
pSAv    = [pSAr,pSVe] ; % absolute source pressure
qSAv    = [pSAr,pSVe] .* YAv; % source flow to node
YNode(:,iAv2Nd)=YNode(:,iAv2Nd)+YAv ; % added ArtVen related conductivity
qNode(:,iAv2Nd)=qNode(:,iAv2Nd)+qSAv; % added ArtVen related inflow
ANode(:,iAv2Nd)=ANode(:,iAv2Nd)+AAv ; % added ArtVen related cross-section


%===== Tube contributions to node pressure and flow
iTb2Nd= [iTb2NdP,iTb2NdD];
YTb   = [YTbP,YTbD] ; % source conductivity
pSTb  = [pSTbP+pExtTb,pSTbD+pExtTb] ; % absolute source pressure
qSTb  = pSTb .* YTb; % source flow to node
Tb2Nd = [Tube2NodeProx;Tube2NodeDist]; % Tube -> matrix -> Node P+D
YNode = YNode + YTb *Tb2Nd; % added Tube related conductivity
qNode = qNode + qSTb*Tb2Nd; % added Tube related inflow
ANode = ANode + ATb * (Tube2NodeProx+Tube2NodeDist); % added cross-section

% If no waterfall:
pNode= qNode ./ YNode; % node pressure by solving NodeInflow Dq=0

% Detect waterfall conditions in ArtVen, insert pressure drop
pAvNd= pNode(:,iAv2Nd); % artven related node pressures
YAvNd= YNode(:,iAv2Nd); % node conductivity
dpAv = dpWf(pAvNd,[pExtAv,pExtAv],pSAv,YAv,YAvNd); %Waterfall pressure drop
dqAv = -YAv.*dpAv; % waterfall node inflow increment per node
pSAv = pSAv-dpAv ; % waterfall pressure correction
qNode(:,iAv2Nd)= qNode(:,iAv2Nd) + dqAv; % waterfall node inflow correction
RgAr = 1:nAv  ; % indices Prox
RgVe = RgAr+nAv; % indices Dist
YAr  = YAv(:,RgAr) ; % arterial conductivity
YVe  = YAv(:,RgVe) ; % venous conductivity
pSAr = pSAv(:,RgAr); % corrected value zero-flow art pressure
pSVe = pSAv(:,RgVe); % corrected value zero-flow ven pressure

% Detect waterfall conditions in Tube, insert pressure drop
pTbNd= pNode(:,iTb2Nd); % node pressures
YTbNd= YNode(:,iTb2Nd); % node conductivity
dpTb = dpWf(pTbNd,[pExtTb,pExtTb],pSTb,YTb,YTbNd); %Waterfall pressure drop
dqTb = -YTb.*dpTb ; % waterfall node inflow increment per node
pSTb = pSTb-dpTb  ; % waterfall pressure correction
qNode= qNode + dqTb * Tb2Nd; % waterfall node inflow correction
RgP  = 1:nTb      ; % indices Prox
RgD  = RgP+nTb    ; % indices Dist
YTbP = YTb(:,RgP) ; % proximal tube conductivity
YTbD = YTb(:,RgD) ; % distal tube conductivity
pSTbP= pSTb(:,RgP); % corrected value zero-flow prox pressure
pSTbD= pSTb(:,RgD); % corrected value zero-flow dist pressure

% Node pressure taking into account waterfalls
pNode = qNode ./ YNode; % Waterfall corrected node pressures

% ArtVen inflow and outflow
pAvNdAr= pNode(:,iAv2NdAr); %pressure Prox-node
pAvNdVe= pNode(:,iAv2NdVe); %pressure Dist-node
qAr    = (pAvNdAr-pSAr).*YAr; % prox tube flow
qVe    = (pSVe-pAvNdVe).*YVe; % dist tube flow

% Tube inflow and outflow
pTbNdP= pNode(:,iTb2NdP); %pressure Prox-node
pTbNdD= pNode(:,iTb2NdD); %pressure Dist-node
qTbP  = (pTbNdP-pSTbP).*YTbP; % prox tube flow
qTbD  = (pSTbD-pTbNdD).*YTbD; % dist tube flow

% d/dt Cavity volume
P.Chamber.VDot= (pNode(:,iCh2Nd) -pSCh ).*YCh ;
P.TriSeg.VLDot= (pNode(:,iTr2NdL)-pSTrL).*YTrL;
P.TriSeg.VRDot= (pNode(:,iTr2NdR)-pSTrR).*YTrR;

P.ArtVen.qAr  = qAr; % ArtVen arterial inflow
P.ArtVen.qVe  = qVe; % ArtVen venous outflow
P.Tube.qProx  = qTbP; % Tube inflow
P.Tube.qDist  = qTbD; % Tube outflow
P.Tube.VDot   = qTbP-qTbD    ; % Tube volume change
P.ArtVen.VArDot= qAr-qAv*Facq; % Art volume change
P.ArtVen.VVeDot= qAv/Facq-qVe; % Ven volume change
P.Node.p      = pNode; % Node pressure
P.Node.Y      = YNode; % Total Node conductivity
P.Node.q      = qNode; % zero pressure node inflow
P.Node.A      = ANode; % total cross-section of connections to node
P.Chamber.p   = pSCh ; % absolute chamber pressures
P.TriSeg.pL   = pSTrL; % absolute TriSeg L cavity pressure
P.TriSeg.pR   = pSTrR; % absolute TriSeg R cavity pressure

end

function dp=dpWf(pN,pE,pS,YS,YN)
% {pN,pE,pS}= pressure {Node, External, Source}
% {YS,YN}= conductivity {Y waterfall, Total Y-Node having YS included}
% dp= decrease of source pressure to simulate waterfall, i.e. limitation
% of flow if outflow pressure pN < external pressure pE
% always dp>0
% if inflow>0, pS<pN, no waterfall
% else if pS<pE, no flow, i.e., dp=pS-pN
% else waterfall: dp=max(0,-PS./(YN./YS-1+eps))
% eps: eps>0, but small, used to avoid zero-division
%=====================
eps= 0.01 ; % safety to avoid zero division
PS = pS-pE; % subtract external pressure
PN = pN-pE; % subtract external pressure
dp = (PS>PN).*max(0,-PS./((PS>0).*(YN./YS-2+eps)+1)-PN);
% Waterfall pressure gradient
end
