function SVar2P
% function SVar2P
% Transfer of scaled state variables in P.SVar to P-structure (SI-units)
% Inverse of P2SVar
% Theo Arts, Maastricht University, Jun 5, 2018

global P
% Scalings
qRef= P.General.q0;
VRef= qRef*P.General.tCycle;

nAv= P.ArtVen.n;
nTb= P.Tube.n;
nCh= P.Chamber.n;
nTr= P.TriSeg.n;
nPa= P.Patch.n;
nVa= P.Valve.n;
a  = cumsum([0,1,nAv,nAv,nTb,nVa,nCh,nTr,nTr,nTr,nTr,nPa,nPa]);
iB = a(1:end-1)+1; iE=a(2:end); % successive begin and end indices

i=1;
P.t         =              P.SVar(:,iB(i):iE(i))  ; i=i+1;
P.ArtVen.VAr=         exp( P.SVar(:,iB(i):iE(i)) ); i=i+1;
P.ArtVen.VVe=         exp( P.SVar(:,iB(i):iE(i)) ); i=i+1;
P.Tube.V    =         exp( P.SVar(:,iB(i):iE(i)) ); i=i+1;
P.Valve.q   = qRef * sinh( P.SVar(:,iB(i):iE(i)) ); i=i+1;
P.Chamber.V =         exp( P.SVar(:,iB(i):iE(i)) ); i=i+1;
P.TriSeg.VL =         exp( P.SVar(:,iB(i):iE(i)) ); i=i+1;
P.TriSeg.VR =         exp( P.SVar(:,iB(i):iE(i)) ); i=i+1;
P.TriSeg.V  = VRef * sinh( P.SVar(:,iB(i):iE(i)) ); i=i+1;
P.TriSeg.Y  =         exp( P.SVar(:,iB(i):iE(i)) ); i=i+1;
P.Patch.Xb  =              P.SVar(:,iB(i):iE(i))  ; i=i+1;
P.Patch.Lsi =              P.SVar(:,iB(i):iE(i))  ; i=i+1;

end
