function SarcSf
% function SarcSf
% Linearizes the relation between myofiber strain Ef and myofiber stress Sf
% for all Patches.
% Linearization parameters are:
% 1) stress Sf for strain Efi=ln(Lsi/2) and
% 2) the derivative dSf/dEf
% The sarcomere is embedded in Patch
% The contractile mechanism is based on the MechChem model
% MechChem thick filament/ C-D difference
% Ref: Dupuis, PLoS 2016 and Arts 2022
% Theo Arts, Maastricht University, Jun 21, 2022

global P

Sarc= P.Patch;

% Ca-Xb properties
hSe = Sarc.hSe  ; % [um] half sarc length series elasticity
kP  = Sarc.kP   ; % exponent ECM stiffness
SfP = Sarc.SfP  ; % ECM stress scaling
SfA = Sarc.SfA  ; % active stress scaling
hi  = Sarc.Lsi/2; % half zero-load sarcomere length (state variable)
Xb  = Sarc.Xb   ; % Xb amount (state variable)
Efi = log(hi)   ; % natural strain 
EXb = Xb.*hi.^2.*SfA/hSe; % Xb Cauchy stiffness

% passive structures
% tension linearization through L=1, isometric
dEf  = hSe./hi; % series elastic strain with L=1
SfEcm= exp(kP.*(Efi+dEf)).*SfP; % ECM stress
EEcm = SfEcm.*kP; % ECM stiffenss
SfEcm= SfEcm-EEcm.*dEf; % SfEcm for L=0

DSfDEf= EXb+EEcm; %summation of Xb and ECM stiffnesses
Sf    = SfEcm; % Sf for L=0 (extrapolation)

Sarc.SfEcm = SfEcm  ; % passive stress ECM for L=0
Sarc.Sf    = Sf     ; % total stress for L=0
Sarc.DSfDEf= DSfDEf ; % sarcomere stiffness
P.Patch = Sarc;

end
