function DisplayPs
% function DisplayPs
% Theo Arts, Maastricht University, Dec 13, 2017
global P

load Ps
p=[];q=[];V=[];Sarc=[];
for j=1:numel(Ps)
    P=Ps(j);
    pj=Get('Node','p',{'La','Ra','Lv','Rv','Ao','PuAr'});
    qj=Get('Valve','q',{'VcRa','RaRv','RvPuAr','PuVeLa','LaLv','LvAo'});
    Vj=[P.Chamber.V,P.TriSeg.VL,P.TriSeg.VR];
    Sarcj=[P.Patch.Ef,P.Patch.Sf];
    p=[p;pj];
    q=[q;qj];
    V=[V;Vj];
    Sarc=[Sarc;Sarcj];
end
CircDisplay
t=(0:size(p,1)-1)'*P.General.dt;
figure(2); plot(t,p)
figure(3); plot(t,V)
figure(4); plot(V,p(:,1:4))
figure(5); plot(t,q)
np=size(Sarc,2)/2;
figure(6); plot(Sarc(:,1:np),Sarc(:,np+(1:np)));
