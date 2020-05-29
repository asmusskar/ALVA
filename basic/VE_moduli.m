function alva = VE_moduli(D0,Dinf,tauD,nD,alva)
%--------------------------------------------------------------------------
% DESCRIPTION:
% This function determines the creep compliance modulus and relaxation
% modulus at selected time increments

% INPUT PARAMETERS
% ti   : Time increments [s]
% tt   : Total travel time [s]
% D0   : Instantaneous or glassy compliance, [1/MPa]
% Dinf : Long time equilibrium or rubbery compliance [1/MPa]
% tauD : Retardation time of the material response [s]
% nD   : Slope of the creep curve in the transient region [-]
%--------------------------------------------------------------------------

ti   = alva.ti;
tt   = alva.tt;
time = zeros(ti,1);
Et   = zeros(ti,1);
E0   = 1/D0;
Einf = 1/Dinf;

for i=1:ti
    if i==1
        time(i) = 0;
    else
        time(i) = tt/10^(ti-i);
    end
    % Dt      = Dinf +((D0-Dinf)/(1+(time(i)/tauD)^nD));
    Et(i)   = Einf*(1+(time(i)/tauD)^nD)/((time(i)/tauD)^nD+(Einf/E0));
end

alva.Et   = Et; 
alva.time = time;