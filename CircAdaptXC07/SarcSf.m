function SarcSf
% function SarcSf
% Linearizes the relation between myofiber strain Ef and myofiber stress Sf
% for all Patches.
% Input: sarcomere length 2hi, Xb-level
% Output: sarcomere stiffness dSf/dEf and Sf, extrapolated to the L=0 state
% Linearization Sf=Sfi+dSf/dEf (Ef-Efi)
% The contractile mechanism is based on the MechChem model
% Ref: Dupuis, PLoS 2016 and Arts, JMCC 2024
% Theo Arts, Maastricht University, March 7, 2024

global P

Sarc= P.Patch;

% Ca-Xb properties
hSe = Sarc.hSe  ; % [um] half sarc length series elasticity
kP  = Sarc.kP   ; % exponent ECM stiffness
SfP = Sarc.SfP  ; % ECM stress scaling
SfA = Sarc.SfA  ; % active stress scaling
hi  = Sarc.Lsi/2; % half zero-Xb load sarcomere length (state variable)
Xb  = Sarc.Xb   ; % Xb amount (state variable)
Efi = log(hi)   ; % natural strain 
EXb = Xb.*hi.^2.*SfA/hSe; % Xb Cauchy stiffness

% passive structures
% tension linearization through L=1, isometric
dEf  = hSe./hi; % series elastic strain with L=1
SfEcm= exp(kP.*(Efi+dEf)).*SfP; % ECM stress
EEcm = SfEcm.*kP; % ECM stiffenss
SfEcm= SfEcm-EEcm.*dEf; % SfEcm, exprapolated to condition L=0

DSfDEf= EXb+EEcm; %summation of Xb and ECM stiffnesses
Sf    = SfEcm; % Sf for L=0 (extrapolation)

P.Patch.SfEcm = SfEcm  ; % passive stress ECM for L=0
P.Patch.Sf    = Sf     ; % total stress for L=0
P.Patch.DSfDEf= DSfDEf ; % sarcomere stiffness

end
