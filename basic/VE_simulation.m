function evres = VE_simulation(letres,alva)
%--------------------------------------------------------------------------
% DESCRIPTION:
% This function construct a matrix of response at different times
% (i.e., for different moduli) used for simulating a moving load

% INPUT PARAMETERS
% dt    : Time increments moving vehicle 
% Xr    : Refined mesh
% Xd    : Output coordinates (x,y,z)
% ti    : Time increments complience curve 
% letres: Linear elastic response in knots
%--------------------------------------------------------------------------

dt   = alva.dt;
Xr   = alva.Xr;
Xd   = alva.Xd;
ti   = alva.ti;
time = alva.time;

% Initialize parameters for refined mesh
rp = zeros(length(Xr),1);             
rr = zeros(length(Xr),length(ti));  
xs  = Xd(:,1);
xr  = Xr(:,1);

% Find response for refined mesh using cubic spline interpolation
for j=1:ti
    rs = letres(:,j);
    for i=1:length(Xr)
        rp(i) = spline(xs,rs,xr(i));
    end
    rr(:,j) = rp;
end

% Mirror response for -x0 to x0
% Response at time t=0...t=tt, i.e. [R(t0), R(t1), R(t2)...R(tt)]
rr =[flipud(rr(2:1:end,:)); rr];   

% Start interpolating beween all the time increments
% -> resulting in a row no.: position / collumn no.:time matrix 
r = zeros(length(rr),length(rr));
for i=1:length(rr)
    t = (i-1)*dt;
    if i==1
        r(:,1)=rr(:,1);
    elseif i > 1 && i < length(rr)
        idt2   = min(find(time > t));
        idt1   = idt2-1;
        r1     = rr(:,idt1);
        r2     = rr(:,idt2);
        r(:,i) = r1 + ((r2-r1)/(time(idt2)-time(idt1)))*(t-time(idt1));
    elseif i==length(rr)
        r(:,i) = rr(:,end);
    end
end

% Simulate moving load: Start point = -x0, Point of interest = 0, 
% i.e. load point = 0 and response point = x0 

% Pseudo-code
% t0: R(E(0),x(0))
% t1: R(E(0),x(dx)) +R(E(dt),x(0)) -R(E(0),x(0))
% t2: R(E(0),x(2dx))+R(E(dt),x(dx))-R(E(0),x(dx))+R(E(2dt),x(0))-R(E(dt),x(0))

evres = zeros(length(rr),1);      % response at evaluation point
for i=1:length(rr)
    k=-(length(rr)-1)+(i-1);
    if i==1
        evres(i)=sum(diag(r,k));
    else
        rp       = diag(r,k);
        rn       = diag(r,k-1);
        evres(i) = sum(rp)-sum(rn);
    end
    kk(i)=k;
end