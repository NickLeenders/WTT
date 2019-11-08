% WTT Assignment 2

%% Q2
Vnom = 690;
Pnom = 4.2E+6;
nPoles = 40;
Ke = 7.8019;
%Ls = 1.2E-3;
Ls = 1.2E-3;
Rs= 4E-3;
Kmech = 373963.3497;

RPM = linspace(6, 15, 20);
omegaMech = 2*pi*RPM/60;
omegaEl = omegaMech*nPoles;
Z = Rs+omegaEl*Ls*1i;

Pmech = Kmech*omegaMech.^3;
dVals = asin(2.*Pmech.*Ls./(Ke.^2.*omegaEl))./2;
Ig = Ke.*sin(dVals)./Ls;
Vg = Ke.*omegaEl.*cos(dVals)-Ig.*Rs;
PlossPhase = Rs.*Ig.^2;

Pout = 3.*Pmech-3*PlossPhase;

% figure()
% plot(RPM, Pout)
% xlabel('RPM')
% ylabel('Generator output power')
% fprintf('Generator active power output at 15 RPM: Pout = %0.3e\n',Pout(end))

%% Q1:
omegaGrid = 100*pi;
omegaGen = 40*pi;

L1 = 2E-6;
R1 = 2E-3;
Z1 = R1+L1*omegaGrid;
Z2 = Z1;
Rc = 1.0;
Lc = 5E-3;
Cc = 1E-6;
ZcPrime = (omegaGrid*Lc*1j+Rc);
Zc = 1j/(omegaGrid*Cc);
Vlow = 690;
Vhigh = 33000;
alpha = Vlow/Vhigh;

Z1prime = Z1*(1/alpha)^2;
Z2prime = Z1prime;
Ztot = Z1prime+Z2prime+ZcPrime;

syms V1 I1 Igrid Vgrid Ic Vc complex

SgridVals = zeros(2, length(Pout));
IgridVals = zeros(2, length(Pout));
VgridVals = zeros(2, length(Pout));
V1Vals = zeros(2, length(Pout));
I1Vals = zeros(2, length(Pout));
for i = 1:length(Pout)
    iPout = Pout(i)/3;
    equations = cell(5,1);
    equations{1, 1} = V1 == Vgrid+(Igrid+Vgrid./Zc).*Ztot;
    equations{2, 1} = I1 == Igrid+Vgrid./Zc;
    
    % Boundary conditions
    equations{3, 1} = iPout == I1.*V1;
    equations{4, 1} = Vgrid == 33000/sqrt(3);
    equations{5, 1} = imag(conj(I1.*V1)) == 0;
    
    Solution = solve(equations{:});
    
    % Removing the complex conjugate here makes the reactive power
    % negative. Maybe this is the way to do it?
    Sgrid = 3*conj(Solution.Igrid).*Solution.Vgrid;
    SgridVals(:,i) = double(Sgrid);
    IgridVals(:,i) = double(Solution.Igrid);
    VgridVals(:,i) = double(Solution.Vgrid);
    V1Vals(:,i)    = double(Solution.V1);
    I1Vals(:,i)    = double(Solution.I1);
end
%plot(RPM, abs(real(SgridVals)))
Pgrid = real(SgridVals);
Qgrid = abs(imag(SgridVals));
%Qcap = -abs(I1Vals-IgridVals).^2*omegaGrid*Cc;

PF = Pgrid./sqrt(Pgrid.^2+Qgrid.^2);
Ploss1 = Ztot.*I1Vals(1, :).^2;
% Note: Apparently, the explanation for RPMs up to ~6 producing higher
% reactive power than active power is that for these low RPMs, any power
% produced is used to charge up the capacitor.
figure()
plot(RPM, Pgrid(1,:))
hold on
plot(RPM, Qgrid(1,:))
grid on
xlabel('RPM','interpreter','latex','FontSize',14)
ylabel('$$P_{poc} [W], Q_{poc} [Var]$$','interpreter','latex','FontSize',14)
legend('Active power [W]','Reactive power [Var]','interpreter','latex','Location','NorthWest','FontSize',12)

figure()
plot(RPM, Pgrid(1,:))
grid on
xlabel('RPM','interpreter','latex','FontSize',14)
ylabel('$$P_{poc} [W]$$','interpreter','latex','FontSize',14)


figure()
plot(RPM, Qgrid(1,:))
grid on
xlabel('RPM','interpreter','latex','FontSize',14)
ylabel('$$Q_{poc} [Var]$$','interpreter','latex','FontSize',14)


%% Q3: Calculate and plot the efficiency
Eff = Pgrid(1,:)./(Pout)*100;

figure()
plot(RPM, Eff)
grid on
xlabel('RPM','interpreter','latex','FontSize',14)
ylabel('$$Efficiency [\%]$$','interpreter','latex','FontSize',14)

figure()
plot(RPM, 3*Pmech)
hold on
plot(RPM, Pgrid(1,:))
grid on
xlabel('RPM','interpreter','latex','FontSize',14)
ylabel('Power [W]','interpreter','latex','FontSize',14)
legend('Kinetic power','POC power','interpreter','latex','Location','NorthWest','FontSize',12)


%% Q4: What is the phase angle between the step up transformer LV voltage 
% and POC voltage.

% Since the POC voltage angle is taken as reference, this angle will
% entirely depend on the angle of V1 with respect to reference (?).

phiV = 180/pi*atan(imag(V1Vals)./real(V1Vals));
figure()
plot(RPM, phiV(1, :))
grid on
xlabel('RPM','interpreter','latex','FontSize',14)
ylabel('$$\phi_{V} [^{o}]$$','interpreter','latex','FontSize',14)


%% Q5: Power factor at POC = 1. Would the system operate with higher or lower
% losses? Discuss.
for i = 1:length(Pout)
    iPout = Pout(i)/3;
    equations = cell(6,1);
    equations{1, 1} = V1 == Vgrid+(Igrid+Vgrid./Zc).*Ztot;
    equations{2, 1} = I1 == Igrid+Vgrid./Zc;
    
    % Boundary conditions
    equations{3, 1} = iPout == real(I1.*V1);
    equations{4, 1} = Vgrid == 33000/sqrt(3);
    equations{5, 1} = imag(Igrid) == 0;
    equations{6, 1} = imag(Vgrid) == 0;
    
    Solution = solve(equations{:});
    
    % Removing the complex conjugate here makes the reactive power
    % negative. Maybe this is the way to do it?
    Sgrid5 = 3*conj(Solution.Igrid).*Solution.Vgrid;
    SgridVals5(:,i) = double(Sgrid5);
    IgridVals5(:,i) = double(Solution.Igrid);
    VgridVals5(:,i) = double(Solution.Vgrid);
    V1Vals5(:,i)    = double(Solution.V1);
    I1Vals5(:,i)    = double(Solution.I1);
end

Pgrid5 = real(SgridVals5);
Qgrid5 = abs(imag(SgridVals5));
Ploss5 = Ztot.*I1Vals5(2, :).^2;
PF5 = Pgrid5./sqrt(Pgrid5.^2+Qgrid5.^2);

Eff5 = Pgrid5(2,:)./(Pout)*100;

figure()
plot(RPM, Eff5)
grid on
xlabel('RPM','interpreter','latex','FontSize',14)
ylabel('$$Efficiency [\%]$$','interpreter','latex','FontSize',14)

figure()
plot(RPM, real(Ploss1-Ploss5))
grid on
xlabel('RPM','interpreter','latex','FontSize',14)
ylabel('$$Efficiency [\%]$$','interpreter','latex','FontSize',14)

