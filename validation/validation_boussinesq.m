function boussinesq = validation_boussinesq(alva)
%--------------------------------------------------------------------------
% DESCRIPTION:
% This function calculates the analytical response for an elastic half
% -space/Boussinesq solution [1].

% INPUT PARAMETERS:
Xd     = alva.Xd;     % Output coordinates (x,y,z)
a      = alva.a;      % Load radii
q      = alva.q;      % Load pressure
E      = alva.E;      % Layer Young's moduli
nu     = alva.nu;     % Layer Poisson's ratios
% -------------------------------------------------------------------------
% References
% -------------------------------------------------------------------------
% [1] Boussinesq, J. (1885). Application des potentiels a l'équilibre et 
% du mouvement des solides élastiques..., volume 4. Gauthier-Villars.
% -------------------------------------------------------------------------

z = Xd(:,3);

sigma_z = zeros(length(z),1); dz = sigma_z;
for i = 1:length(z)
    sigma_z(i) = q*(1-(1/(sqrt(1+(a/z(i))^2)^3)));
    dz(i)      = ((1-nu(1)^2)*2*q*a/E(1))*(sqrt(1+(z(i)/a)^2)-z(i)/a)...
                *(1+(z(i)/a)/(2*(1-nu(1))*sqrt(1+(z(i)/a)^2)));
end

sigma_z = sigma_z;
dz      = dz;

boussinesq = [sigma_z dz];
