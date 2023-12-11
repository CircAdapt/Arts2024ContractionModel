function Indexation
% function Indexation
% Sets backbone of structure P
% Naming the elements ArtVen, TriSeg, Chamber, Valve, Tube, Node,
% Wall, Patch and DepPath
% Connections between elements are defined by strings stored in
% structure field P.Map
%
% The P-Structure is built on the information in P.Map
% ArtVen represents artery-peripheral resistance-vein of organ or body part
% Chamber represents a cavity enclosed by a myocardial wall (atria)
% TriSeg represents combination of two cavities with three walls
% Bag represents passive elastic bag, encapsulating part of the
% circulation, like the pericardium
% Node: named connection point
% Wall: contractile muscular wall, composed of 1 or more patches
% Patch: contractile part of a wall, having specific mechanical properties
% Valve: valve with inertia connects proximal to distal node, may leak
% Tube: elastic tube connects proximal to distal node
% DepPath: depolarization pathway between patches, acting as nodes
%
% Theo Arts, Maastricht University, Nov 27, 2021

global P;
Map= P.Map; % contains complete information about element composition 
% and mutual connections of Heart and Circulation

%==== Creating name-strings for elements =============

% ArtVen->Node
ArtVenName = Map.ArtVenName; % getting ArtVen names
Aux        =[strcat(ArtVenName,{'Ar'});strcat(ArtVenName,{'Ve'})];
NodeName   = Aux(:)'; % names of arterial and venous coupled cavities
ArtVen2Node= Aux(1,:)'; % refers to arterial node

% Chamber->[Node,Wall]
ChamberName  = Map.ChamberName; % getting chamber names
NodeName     = [NodeName,ChamberName]; % append related cavities
Chamber2Node = ChamberName;
Chamber2Wall = ChamberName;
WallName     = ChamberName; % chamber related walls

% TriSeg->[Node,Wall]
TriSegName   = Map.TriSegName; % getting TriSeg name(s)
Aux          = [strcat({'L'},TriSegName);strcat({'R'},TriSegName)];
NodeName     = [NodeName,Aux(:)']; % Append TriSeg cavities to Node record
TriSeg2Node  = Aux(1,:);
Aux=[strcat({'L'},TriSegName);...
    strcat({'S'},TriSegName);...
    strcat({'R'},TriSegName)];
WallName     = [WallName,Aux(:)']; % Name of L, S and R wall
TriSeg2Wall  = Aux(1);

[NodeName,ValveName,Valve2NodeProx,Valve2NodeDist]=...
    DNode(NodeName,Map.ValveNodes);

[NodeName,TubeName,Tube2NodeProx  ,Tube2NodeDist ]=...
    DNode(NodeName,Map.TubeNodes);

% Wall->Patch, Wall.nPatch
Wall2Patch= strcat(WallName,'1'); % Patches are named to wall and numbered
PatchName = Wall2Patch;
WallnPatch= ones(1,length(WallName)); % 1st patch of each wall
% Walls composed with MultiPatch
MultiPatchWall= Map.MultiPatch; % search multi-patched walls
nMultiPatch   = length(MultiPatchWall);
for i=1:nMultiPatch
    NameMulti= {};
    Str      = MultiPatchWall{i};
    LtrStr   = isletter(Str);%logical letter positions
    NameStr  = Str(LtrStr); %name of multipatched wall
    nPatch   = str2double(Str(~LtrStr)); %number of patches in multipatch
    iWall    = strcmp( NameStr     ,WallName )*(1:numel(WallName))';
    iPatch   = strcmp([NameStr,'1'],PatchName)*(1:numel(PatchName))';
    for k=2:nPatch
        NameMulti=[NameMulti,[NameStr,num2str(k)]];
    end
    PatchName=[PatchName(1:iPatch),NameMulti,PatchName(iPatch+1:end)];
    WallnPatch(iWall)= nPatch; % store number of patches in wall record
end

% DepPath, Patch~Node, DepPath connects Patches directionally
[PatchName,DepPathName,DepPath2PatchProx,DepPath2PatchDist ]=...
    DNode(PatchName, Map.DepPath);

%==== END Creating name-strings for the elements =============

%==== indexations, numbering of the name to enhance calulation speed =====
% Naming and counting of elements, determined by name strings
P.ArtVen = Naming('ArtVen' ,ArtVenName );
P.Chamber= Naming('Chamber',ChamberName);
P.TriSeg = Naming('TriSeg' ,TriSegName );
P.Valve  = Naming('Valve'  ,ValveName  );
P.Tube   = Naming('Tube'   ,TubeName   );
P.Node   = Naming('Node'   ,NodeName   );
P.Wall   = Naming('Wall'   ,WallName   );
P.Patch  = Naming('Patch'  ,PatchName  );
P.DepPath= Naming('DepPath',DepPathName);

% indices determine mutual relations between elements
P.ArtVen.iNode      = Get('Node'  ,'Index',ArtVen2Node   );
P.Chamber.iWall     = Get('Wall'  ,'Index',Chamber2Wall  );
P.Chamber.iNode     = Get('Node'  ,'Index',Chamber2Node  );
P.TriSeg.iWall      = Get('Wall'  ,'Index',TriSeg2Wall   );
P.TriSeg.iNode      = Get('Node'  ,'Index',TriSeg2Node   );
P.Valve.iNodeProx   = Get('Node'  ,'Index',Valve2NodeProx);
P.Valve.iNodeDist   = Get('Node'  ,'Index',Valve2NodeDist);
P.Tube.iNodeProx    = Get('Node'  ,'Index',Tube2NodeProx );
P.Tube.iNodeDist    = Get('Node'  ,'Index',Tube2NodeDist );
P.Wall.nPatch       = WallnPatch;
P.Wall.iPatch       = Get('Patch' ,'Index',Wall2Patch    );
P.DepPath.iPatchProx= Get('Patch' ,'Index',DepPath2PatchProx );
P.DepPath.iPatchDist= Get('Patch' ,'Index',DepPath2PatchDist );

% Indexation of Bag-content
BagIndexation

% Baroreceptor and pacemaker location
P.Node.iBaro    = Get('Node','Index',P.Map.Baro); %Baro receptor
P.Patch.iPace   = Get('Patch','Index',P.Map.Pace); % Pacing patch
P.DepPath.iCycle= find(P.DepPath.iPatchProx==P.Patch.iPace & ...
    P.DepPath.iPatchDist==P.Patch.iPace); % Circular DepPath=pace maker
P.DepPath.iAv   = Get('DepPath','Index',P.Map.AvDepPath); % AV DepPath

% Creates index-based sparce matrices to define inter-element links
LinkMatrix

end


% ======== AUXILARY FUNCTIONS ============

function Element= Naming(ElementType,ElementName)
% Naming
% Naming and number of elements stored in P-structure
global P;
Element=P.(ElementType);
Element.Name= ElementName;
Element.n   = length(Element.Name);
end

function [NodeName,ElementName,Prox,Dist]=DNode(NodeName,StrNode)
% From data in (P.)Map not yet used node names nodes are added to already
% used nodes
% In:
% NodedName= name of all already used nodes
% StrNode  = cell array with rows of prox and dist node (maybe known) names
% Out:
% NodeName   : Nodes with additions
% ElementName: Created element name
% Prox,Dist  : Created proximal and distal noded name

Prox= StrNode(:,1)';
Dist= StrNode(:,2)';
ElementName=strcat(Prox,Dist);
dNodeName= setdiff([Prox,Dist],NodeName);
% remove double names 
NodeName = [NodeName,dNodeName];
end

function BagIndexation
% function BagIndexation
% Bags may encapsulate ArtVen-, Tube-, Chamber-, TriSeg- and Bag- elements
% Nesting of bags requires attention
% Elements connections are represented by matrices.
global P
MapBag      = P.Map.Bag; %reading structure definition in P.Map.Bag
BagName     = MapBag.Name;
P.Bag.Name  = BagName; % Bag Name is needed for Bag identification by 'Get'
Bag         = P.Bag;
Bag.n       = length(BagName);
Bag.OK      = zeros(size(BagName)); % keeps track of nesting of Bags
CellVec     = cell(1,Bag.n);
Bag.iChamber= CellVec;
Bag.iTriSeg = CellVec;
Bag.iArtVen = CellVec;
Bag.iTube   = CellVec;
Bag.iBag    = CellVec;
for ib=1:Bag.n % fill in indices
   Bag.iChamber{ib}= Get('Chamber','Index',MapBag.Chamber{ib});
   Bag.iTriSeg{ib} = Get('TriSeg' ,'Index',MapBag.TriSeg{ib} );
   Bag.iArtVen{ib} = Get('ArtVen' ,'Index',MapBag.ArtVen{ib} );
   Bag.iTube{ib}   = Get('Tube'   ,'Index',MapBag.Tube{ib}   );
   Bag.iBag{ib}    = Get('Bag'    ,'Index',MapBag.Bag{ib}    );
end

Bag.OK=cellfun(@isempty,Bag.iBag); % check on nesting of bags
% find bags in bag:
while ~all(Bag.OK) % label OK marks if taken care for nesting
    Rg=find(~Bag.OK); % marked not OK
    for i=1:numel(Rg) % indices contained Bags
        ib1=Rg(i); % Outer bag
        ib2=cell2mat(Bag.iBag(ib1)); % inner Bag(s) of ib1
        if all(Bag.OK(ib2)) % in inner bags OK
            % merge inner bag to outer bag
            Bag.OK(ib1)=1; % this Bag becomes OK
            for j=1:numel(ib2); % add encapsulated elements of inner bag
                ib=ib2(j);
                Bag.iChamber{ib1}= [Bag.iChamber{ib1}, Bag.iChamber{ib}];
                Bag.iTriSeg{ib1} = [Bag.iTriSeg{ib1} , Bag.iTriSeg{ib} ];
                Bag.iArtVen{ib1} = [Bag.iArtVen{ib1} , Bag.iArtVen{ib} ];
                Bag.iTube{ib1}   = [Bag.iTube{ib1}   , Bag.iTube{ib}   ];
            end
        end
    end
end

P.Bag= Bag;
end

function LinkMatrix
% function LinkMatrix
% Transfers linking information from strings into indices.
% Links may be defined by sparse matrices or cell arrays
% Matrices are constructed, based on information, described by strings
% refering to element names, created by the beginning of 'indexation'
% Main purpose: use of indices speeds up calculation, because string
% reading is slow
% Theo Arts, Maastricht University, April 10, 2018

global P
nWall    = P.Wall.n;
nPatch   = P.Patch.n;
nChamber = P.Chamber.n;
nTriSeg  = P.TriSeg.n;
nArtVen  = P.ArtVen.n;
nValve   = P.Valve.n;
nTube    = P.Tube.n;
nNode    = P.Node.n;
nBag     = P.Bag.n;
nDepPath = P.DepPath.n;

% sparse matrix for Valve <-> Node connections
% robust to nValve==0
Aux= sparse(zeros(nValve,nNode));
P.Valve.Valve2NodeProx=Aux;
P.Valve.Valve2NodeDist=Aux;
Rg= 1:nValve;
AuxProx= Rg+(P.Valve.iNodeProx(Rg)-1)*nValve;
AuxDist= Rg+(P.Valve.iNodeDist(Rg)-1)*nValve;
P.Valve.Valve2NodeProx(AuxProx)=+1;
P.Valve.Valve2NodeDist(AuxDist)=+1;

% sparse matrix for Tube <-> Node connections
% robust to nTube==0
Aux= sparse(zeros(nTube,nNode));
P.Tube.Tube2NodeProx=Aux;
P.Tube.Tube2NodeDist=Aux;
Rg= 1:nTube;
AuxProx=Rg+(P.Tube.iNodeProx(Rg)-1)*nTube;
AuxDist=Rg+(P.Tube.iNodeDist(Rg)-1)*nTube;
P.Tube.Tube2NodeProx(AuxProx)=+1;
P.Tube.Tube2NodeDist(AuxDist)=+1;

% sparse matrix for Patch <-> Wall connections
Aux= sparse(nPatch,nWall);
j=0;
for i=1:nWall
    j=P.Wall.iPatch(i)-1;
    Aux(j+(1:P.Wall.nPatch(i)),i)=1;
end
P.Wall.Wall2Patch=Aux';

% Valve-Wall connections (papillary muscles)
P.Valve.AvValves=Get('Valve','Index',{'LaLv','RaRv'});
P.Valve.AvWalls =Get('Wall','Index',{'Lv','Rv'});

% Bag indexings
% sparse matrices for volume contributions to Bag by:
% Wall-VW, Tube-VT concatenated to VCWT
% sparse matrices for pressure contributions from Bag to:
% Chamber-Ch, TriSeg-Tr, ArtVen-Av, Tube-Tb, Bag-Bg
Ch=sparse(nChamber,nBag);
Tr=sparse(nTriSeg,nBag);
Av=sparse(nArtVen,nBag);
Tb=sparse(nTube,nBag);
Bg=sparse(nBag,nBag);

for iBag=1:nBag
    Ch(P.Bag.iChamber{iBag},iBag)=1; %Chamber link matrix
    Tr(P.Bag.iTriSeg{iBag} ,iBag)=1; %TriSeg link matrix
    Av(P.Bag.iArtVen{iBag} ,iBag)=1; %ArtVen link matrix
    Tb(P.Bag.iTube{iBag}   ,iBag)=1; %Tube link matrix
    Bg(P.Bag.iBag{iBag}    ,iBag)=1; %Bag-Bag link matrix
end
P.Bag.Ch= Ch;
P.Bag.Tr= Tr;
P.Bag.Av= Av;
P.Bag.Tb= Tb;
P.Bag.Bg= Bg;

% linking from patch to entrance depolarization pathways
iPatchProx= P.DepPath.iPatchProx; % proximal patches
Aux= cell(1,nPatch); % cell array refers to pathway indices per patch
for iD=1:nDepPath
    Aux{iPatchProx(iD)}=[Aux{iPatchProx(iD)},iD];
end
P.Patch.iDepPathProx= Aux; % Patch -> beginning of new DepPath

end

