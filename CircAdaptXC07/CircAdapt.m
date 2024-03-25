function CircAdapt
% function CircAdapt
% Core of the CircAdapt model
% Simulates a series of beats based on de contents of structure P
% Theo Arts, Maastricht University, April 26, 2014

global P
P.General.dt= P.General.Dt;%init
Indexation; % Indexation to define connections between elements
% and matrices defining connections between elements

% Initialization P.SVar (= state variables) to last state
P2SVar;
nBeat= P.General.nStoredBeats;% max nr of stored P-structures

% clearing In-Out record to judge steady state
P.General.In     = [];
P.General.Out    = []; % reset of storage input/output per beat

%==== simulation of successive cardiac cycles
P.General.BeatNr= 0;
Ps=[]; % init of array of saved P-structures
while 0<P.General.tEnd % new beat
    TimingInit; % Setting depolarization sequence
    % Also setting pacemaker to tCycle
    % P.Patch.tEnd marks upcoming time interval length (~tCycle)
    % Also, Tube delay-signals are shifted in time

    tEnd     = P.Patch.tEnd; % end time of current beat set by DepPath
    dt       = P.General.dt; % used time step
    SVar     = P.SVar(end,:); % SVar= row vector, boundary cond. t=0
    SVar(1)  = 0; % setting P.t=0 for each beat
        
    % just simulated beat
    disp(['Time to go= ', num2str(P.General.tEnd)]); pause(0.01);
    P.SVar = OdeCA('SVarDot',dt,tEnd,SVar); % ODE-solver
    % OdeCA is a simple ODE solver, taylored to CircAdapt
    CircDisplay; % displays time course of 1-beat hemodynamics
    
    Ps   = [Ps(max(1,end-nBeat+2):end),P]; % maximized number of stored bts
    if nBeat>0; save Ps Ps; end % saving the last beats
    
    feval(P.General.AdaptFunction); % specific adaptation function
    P.General.tEnd= P.General.tEnd-P.t(end); % remaining simulation time    
end
end
