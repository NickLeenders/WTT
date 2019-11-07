% WTT Assignment 2

%% Q2
Vnom = 690;
Pnom = 4.2E+6;
nPoles = 40;
Ke = 7.8019;
Ls = 1.2E-3;
Rs= 4E-3;
Kmech = 373963.3497;

RPM = linspace(3, 15, 20);
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
Zc = -1j/(omegaGrid*Cc);
Vlow = 690;
Vhigh = 33000;
alpha = Vlow/Vhigh;

Z1prime = Z1*(1/alpha)^2;
Z2prime = Z1prime;
Ztot = Z1prime+Z2prime+ZcPrime;



syms V1 I1 Igrid Vgrid Ic Vc complex
% equations = cell(6,1);
% 
% equations{1, 1} = I1-Ic+Igrid == 0;
% equations{2, 1} = -I1*Ztot-Ic*Zc+V1 == 0;
% equations{3, 1} = -I1*Ztot+V1-Vgrid == 0;
% equations{4, 1} = I1*V1 == Pout;
% equations{5, 1} = Vgrid == 33000;
% equations{6, 1} = Ic == Vc/Zc;
% 
% Solution = solve(equations{:});
% 
% Sgrid = conj(Solution.Igrid).*Solution.Vgrid;
% Sgrid = subs(Sgrid)
SgridVals = zeros(2, length(Pout));
IgridVals = zeros(2, length(Pout));
VgridVals = zeros(2, length(Pout));
V1Vals = zeros(2, length(Pout));
I1Vals = zeros(2, length(Pout));
for i = 1:length(Pout)
    iPout = Pout(i);
    equations = cell(7,1);
    equations{1, 1} = I1-Ic+Igrid == 0;
    equations{2, 1} = -I1.*Ztot-Ic.*Zc+V1 == 0;
    equations{3, 1} = -I1.*Ztot+V1-Vgrid == 0;
    equations{4, 1} = I1.*V1 == iPout;
    equations{5, 1} = Vgrid == 33000;
    equations{6, 1} = Ic == Vc/Zc;
    equations{7, 1} = real(conj(Igrid).*Vgrid)<iPout;
    
    Solution = solve(equations{:});
    Sgrid = conj(Solution.Igrid).*Solution.Vgrid;
    SgridVals(:,i) = eval(subs(Sgrid));
%     IgridVals(:,i) = eval(subs(Igrid));
%     VgridVals(:,i) = eval(subs(Vgrid));
%     V1Vals(:,i)    = eval(subs(V1));
%     I1Vals(:,i)    = eval(subs(I1));
end
%plot(RPM, abs(real(SgridVals)))
Pgrid = abs(real(SgridVals));
Qgrid = abs(imag(SgridVals));

PF = Pgrid./sqrt(Pgrid.^2+Qgrid.^2);
figure()
plot(RPM, Pgrid(1,:))
hold on
plot(RPM, Qgrid(1,:))




