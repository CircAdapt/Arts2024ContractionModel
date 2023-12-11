% Script DisplayWave
% Theo Arts, Maastricht University, Oct 13, 2012
% Displays wave pressures and flow velocities, based on structure P

pAr={'Lv','Ao','BrAr','CaAr','CeAr','FeAr'};
pVe={'Ra','Vc','BrVe','CeVe','FeVe'};
TubeArtFlow={'AoBrAr','BrArCeAr','CeArFeAr'};
TubeVenFlow={'VcBrVe','BrVeCeVe','CeVeFeVe'};
ArtFlowDist={'CaAr','BrAr','CeAr','FeAr'};

pThorax=Get('Bag','p','Thorax');
pPeri  =Get('Bag','p','Peri');
t=P.t;

figure(2); M=Get('Node','p',pAr);
plot(t,M);
axis([min(t),max(t),min(M(:)),max(M(:))])
legend(pAr)
title('Arterial Pressures')

figure(3); M=[Get('Node','p',pVe),pThorax,pPeri];
plot(t,M);
axis([min(t),max(t),min(M(:)),max(M(:))])
legend([pVe,{'Thorax','Peri'}])
title('Venous Pressures')

M1= [Get('Tube','qProx',TubeArtFlow),Get('Tube','qDist',TubeArtFlow)];
A= Get('Tube','A',TubeArtFlow);
M=M1./[A,A];
figure(4);
plot(t,M);
axis([min(t),max(t),min(M(:)),max(M(:))])
legend([TubeArtFlow,TubeArtFlow])
title('Prox and Dist Artery Velocity')

M1= [-Get('Tube','qProx',TubeVenFlow),-Get('Tube','qDist',TubeVenFlow)];
A= Get('Tube','A',TubeVenFlow);
M= M1./[A,A];
figure(5);
plot(t,M);
axis([min(t),max(t),min(M(:)),max(M(:))])
legend([TubeVenFlow,TubeVenFlow])
title('Prox and Dist Venous Velocity')

Len= Get('ArtVen','Len','All');
nt=length(P.t);
AAr= P.ArtVen.VAr./Len;
M  = P.ArtVen.qAr./AAr;
figure(6);
plot(t,M);
axis([min(t),max(t),min(M(:)),max(M(:))])
legend(P.ArtVen.Name)
title('Arterial Flow Velocity into Periphery')

