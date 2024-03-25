function pBag
% function pBag
% Bags enclose Chambers, Trisegs, Walls, Tubes and other Bags
% Transmural pressures (pTrans) are added to absolute pressures (p)
% Bags may be nested
% Theo Arts, Maastricht University, June 8, 2018

global P
Bag  = P.Bag;

if P.Bag.n>0
    Ch= Bag.Ch; % Linking matrices
    Tr= Bag.Tr;
    Av= Bag.Av;
    Tb= Bag.Tb;
    Bg= Bag.Bg;
    
    % volumes of enclosed elements
    VWall   = P.Wall.VWall;
    VChamber= P.Chamber.V+VWall(P.Chamber.iWall);
    iW      = P.TriSeg.iWall;
    Aux     = VWall(iW)+VWall(iW+1)+VWall(iW+2);
    VTriSeg = P.TriSeg.VL+P.TriSeg.VR+Aux;
    VTube   = P.Tube.V + P.Tube.AWall .* P.Tube.Len;
    VArtVen = P.ArtVen.VAr +P.ArtVen.VVe;
    % Bag pressure volume relation
    VRef    = P.Bag.VRef  ;
    pAdapt  = P.Bag.pAdapt;
    k       = P.Bag.k;
    % Total enclosed volume per bag
    VBag    = VChamber*Ch+VTriSeg*Tr+VTube*Tb+VArtVen*Av;
    pTrans  = (VBag./VRef).^k.*pAdapt;

    % Pressure increments due to Bag-pTrans
    % These increments are external pressures on Cavity, Tube and Bag
    DpCh= pTrans * Ch'; %Chamber
    DpTr= pTrans * Tr'; %TriSeg
    DpTb= pTrans * Tb'; %Tube
    DpAv= pTrans * Av'; %ArtVen
    DpBg= pTrans * Bg'; %ArtVen
    
    % Bag output
    P.Bag.V     = VBag;
    P.Bag.pTrans= pTrans;
    P.Bag.p     = pTrans + DpBg;
else % no bags
    DpCh= 0;
    DpTr= 0; 
    DpTb= 0;
    DpAv= 0;
    DpBg= 0;
end

% External pressures
P.Chamber.pExt= DpCh;
P.TriSeg.pExt = DpTr;
P.Tube.pExt   = DpTb;
P.ArtVen.pExt = DpAv;
P.Bag.pExt    = DpBg;

end
