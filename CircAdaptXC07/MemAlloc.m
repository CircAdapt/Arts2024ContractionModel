function MemAlloc
% function MemAlloc
% Allocates memory for matrices with column length= number of time points
% These allocations are necessary because elements in matrices get a value
% by gradual substitutions, using indices for locating them in the matrix.
% Theo Arts, Maastricht University, June 8, 2018

global P;

nt=size(P.t,1); % number of time points

% test on change of column length, only then MemAlloc is needed
if isfield(P,'SVarDot')
    ExecuteMemAlloc = (size(P.SVar,1)~=size(P.SVarDot,1));
else
    ExecuteMemAlloc = 1;
end
    
if ExecuteMemAlloc
    
    nWall  =P.Wall.n; % number of walls
    nTriSeg=P.TriSeg.n; % number of TriSeg objects
    
    % Matrix constructions
    MatWall  = zeros(nt,nWall);
    MatTriSeg= zeros(nt,nTriSeg);
    
    % Allocating matrices for the walls
    P.Wall.Aw0   = MatWall; %unstressed mid-wall area
    P.Wall.DADT  = MatWall; %area compliance
    P.Wall.T     = MatWall; %tension
    P.Wall.Cm    = MatWall; %curvature
    P.Wall.Aw    = MatWall; %actual mid-wall area
    P.Wall.pTrans= MatWall; %transmural pressure

    % Allocating matrices for the TriSeg objects
    P.TriSeg.VS     = MatTriSeg; % volume displacement septal wall
    P.TriSeg.YS     = MatTriSeg; % radius junction ring
    P.TriSeg.VDot   = MatTriSeg; % SVar-Dot
    P.TriSeg.YDot   = MatTriSeg; % SVar-Dot
    P.TriSeg.AL     = MatTriSeg; % LV-cross-section
    P.TriSeg.AR     = MatTriSeg; % RV-cross-section
    P.TriSeg.ZL     = MatTriSeg; % LV source impedance
    P.TriSeg.ZR     = MatTriSeg; % RV source impedance
    P.TriSeg.pTransL= MatTriSeg; % LV transmural pressure
    P.TriSeg.pTransR= MatTriSeg; % RV transmural pressure
end
end

