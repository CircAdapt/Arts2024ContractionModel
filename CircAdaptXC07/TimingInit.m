function TimingInit
% function TimingInit
% Sets patch depolarization times for upcoming beat
% Theo Arts, Maastricht University, Oct 24, 2018

global P

iPace       = P.Patch.iPace; % index of pacing patch
TauRefrac   = P.Patch.TauRefrac; % refractory periods of patch
nP          = P.Patch.n; % number of patches
tDep        = P.Patch.tDep; % pair of depolarization times per Patch
iDepPathProx= P.Patch.iDepPathProx; % index to DepPath
Delay       = P.DepPath.Delay; % delays in the depolarization pathways 
iPatchProx  = P.DepPath.iPatchProx; 
iPatchDist  = P.DepPath.iPatchDist;
tCycle      = P.General.tCycle; % cycle time
tCycleRef   = P.General.tCycleRef; % reference cycle time
TauAvRef    = P.General.TauAvRef; % AV-delay at rest
TauAvExc    = P.General.TauAvExc; % AV-delay with exercise
iCycle      = P.DepPath.iCycle; % tCycle-related pacemaker circular DepPath
iAv         = P.DepPath.iAv; % AV-delay related DepPath

Much = 100; % serves as 'infinite' time
nk   = 1  ; % if nk==2 -> recalculation tDep (needed for 1st beat)

if numel(tDep)==0 % if there is no history of previous beat
    Dep      = zeros(1,nP); % boolean: depolarization occurred
    t0       = Dep-Much;
    t1       = Dep+Much;
    t        = 0; % time (progressing from dep to dep
    t0(iPace)= t; % depolarization starts with pacing-patch
    nk       = 2; % extra cycle to estimate steady state dep history
else % rolling tDep-memory by one step
    tA= tDep(1,:); % old depolarization time
    t0= tDep(2,:); % last depolarization, occurred in previous beat
    t1= 0*t0+Much; % expected new, yet unknown depolarization times
end

t= t0(iPace); % starting time, maybe non-zero in case of a fast circle
% Depolarization pathways
for kk=1:nk % if no history of previous beat: nk=2, else: nk=1
    t0   = t0-t; % time shift
    Dep  = 0*t0; %Boolean detecting depolarization in upcoming beat
    Dep(iPace)= 1; % first dep
    t=0;
    tIn = t0(iPatchProx); % starting time of DepPath
    tOut= tIn+Delay; % ending time of DepPath
    RgD = find( t>=tIn & t<tOut); % active delays
    EndBeat= 0; % Boolean marking end of beat interaval
    tA=t0; % known last depolatization times
    while ~EndBeat
        tOutA = tOut(RgD); % tOut for active delays
        RgP   = iPatchDist(RgD); % Patches distal to active delays
        tDepA = t0(RgP); % last depolarizations of distal patches
        tRefrA= tDepA+TauRefrac(RgP); % refractory time of distal patches
        t1A=t1(RgP); % expected depolarization time
        
        B1= (tOutA>tRefrA) & (tOutA<t1A); % if a faster pathway is found
        while any(B1) % iteration until all relevant pathways are included
            t1A(B1)=tOutA(B1);
            B1= (tOutA>tRefrA) & (tOutA<t1A);
        end
        t1(RgP)   = t1A; % candidate time of earliest depolarization
        t         = min(t1A); % setting time to current depolarization
        iPace     = RgP(t1A==t); % Depolarization of Patch
        t1(iPace) = Much; % expected next depolarization time
        EndBeat   = any(Dep(iPace)); %if 2nd depolarization occurs
        Dep(iPace)= 1; % record depolarization on patch
        t0(iPace) = t; % new depolarization time
        Rg        = cell2mat(iDepPathProx(iPace)); % just activated delays
        tIn(Rg)   = t; % upcoming activated delays
        iD        = RgD(tOut(RgD)>t); % remainder of active delays
        iD1       = Rg; % added delays
        tOut(Rg)  = tIn(Rg)+Delay(Rg); % pathway delayed times
        RgD       = unique([iD,iD1]) ; % active pathways
    end
end
tEnd         = t;
TauAv        = TauAvRef+2*(tCycle/tCycleRef-1)*(TauAvRef-TauAvExc);
Delay(iCycle)= tCycle; % Delay(Ra1->Ra1)
Delay(iAv   )= TauAv ; % Delay(Ra1-Sv1)

% Shifting/preparing Tube-delayed signals
nt  = ceil(tEnd/P.General.Dt); % number of time points upcoming beat
dt  = tEnd/(nt-1); % set integer number of time steps per cycle
P.General.dt= dt; %slightly changed dt to get integer number of t-steps

t1   = P.t; % previous beat
pL1  = P.Tube.pL; % Wave signals previous beat may enter current beat
pR1  = P.Tube.pR;
uP1  = P.Tube.uP;
uD1  = P.Tube.uD;
q1   = P.Tube.q ;

t2   = (0:nt-1)'*dt; % time samples next beat
t12  = t1-t1(end)+tEnd; % time shift previous->current beat

% Delayed signals interpolated in time,relevant with change of tCycle or Dt
P.Tube.pL= interp0(t12,pL1,t2);%left  wave used as delay line
P.Tube.pR= interp0(t12,pR1,t2);%right wave used as delay line
P.Tube.uP= interp0(t12,uP1,t2);%proximal zero flow pressure
P.Tube.uD= interp0(t12,uD1,t2);%distal zero flow pressure
P.Tube.q = interp0(t12,q1 ,t2);%mean tube flow
P.Tube.dt= dt; %sampling interval Tube delay lines
% end Tube signal preparation

P.DepPath.Delay= Delay; % delays per DepPath
P.Patch.tDep   = [tA;t0]; % time of depolarization per Patch
P.Patch.tEnd   = tEnd; % end of current depolarization time span
end

function p2=interp0(t1,p1,t2)
% converts time scale of previous beat to that of current beat
% if interpolation out of interval, boundary point is used
if numel(t1)==1
    t1=[0;t1];
    p1=p1([1,1],:);
end
t1(1)=min(0,t1(1));
p2=interp1(t1,p1,t2);
end
