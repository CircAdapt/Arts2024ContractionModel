function AdaptExc
% function AdaptExc
% Starts with
% Standard beat to beat (absence of) adaptation
% Blood pressure is controlled by change of circulating volume
% followed by specific exercise adaptation
% Theo Arts, Maastricht University, Mar 5, 2024

global P

Adapt;
AdaptE; %SPECIFIC ADAPTATION

P2SVar; % load adapted values in P.SVar for initiation of next beat
end


function AdaptE
%function AdaptExc;
% Specific adaptation to exercise conditions
% Control of:
% Blood vessel wall thickness by wall stress
% Heart wall thickness and cavity volume by myocardial tissue load
% Bag pressure
% Theo Arts, Maastricht University, Feb 20, 2013

%=== Begin adapt special
%=== Actions of adaptation at excercise
%=== Adapt ArtVen wall thickness and Patches
AdaptArtVen('All','All');
AdaptTube('All','All');
AdaptPatch('All','All');
AdaptBag;
%=== end actions of adaptation

end

%==== Specific adaptation functions ==========================

function AdaptBag %  Bag adaptation to pAdapt reference
global P
if P.Bag.n==0; return; end
pMax      = max(P.Bag.pTrans);
Fac_pMax  = max(0.1,pMax./P.Bag.pAdapt);
P.Bag.VRef= P.Bag.VRef.*(Fac_pMax.^(0.3 ./ P.Bag.k));
disp(['Bag Adaptation ',[P.Bag.Name{:}],'                x1000: ',...
    num2str(1000*log(Fac_pMax),'%7.0f')]);

end
