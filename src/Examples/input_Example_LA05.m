%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example LA 05: LA of spherical shell segment under radially outwards edge load.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% [1] - Flügge W., Stresses in Shells, 2nd edition, Springer-Verlag, Berlin Heidelberg, Germany, 1973
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.

%%%%%%%%%%%%%%%%%
% Preliminaries %
%%%%%%%%%%%%%%%%%
clearvars, clc, close all
addpath('..','../AQUINAS_Analysis','../AQUINAS_Objects','../Maple_Source');

%%%%%%%%%%
% Inputs %
%%%%%%%%%%

% Note of caution! : The analytical solution provided here is only valid for small values of kappa, as obtained from eq.(6.20) of [1].
% For further details see notes on limitations of analytical solution, provided by eqs.(6.26) of [1].

% Geometry
a = 381.0; % [mm] - Sphere meridional radius
t = 25.4; % [mm] - Sphere thickness
alpha = 10*pi/180; % [rad] Half angle of spherical shell section (see Fig 3. of Chapter 6 of [1])

% Material
MatName = 'Steel'; % Material name
E = 2.0e5; % [N/mm2] - Young's modulus of steel
nu = 0.3; % [-] - Poisson ratio of steel

% Loading
H = 1; % [N/mm] Radial edge load, acting at both exposed edges of the spherical shell section (see Fig 3. of Chapter 6 of [1])

% Shell segment FE discretization
numelem = 100;

% Plot controls
FS = 20; % [-] - Font Size (pts)
FSL = 18; % [-] - Font Size for Legends (pts)
MS = 3; % [-] - Markser size (pts)
LW = 3; % [-] - Line Width (pts)


%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%

M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)));

S1 = AQUINAS_Segment_Object('type','Ellipse','material','Steel','thickness',t,'rzbot',[a  0],'rztop',[a*sin(pi/2-alpha) a*cos(pi/2-alpha)],'geom',[0 0 a a],'els',numelem);

C1 = AQUINAS_Constraint_Object('rcoord',a,'zcoord',0,'dofs',{'w'});
C2 = AQUINAS_Constraint_Object('rcoord',a,'zcoord',0,'dofs',{'b'});

F1 = AQUINAS_Nodal_Force_Object('type','Fu','rcoord',a*sin(pi/2-alpha),'zcoord',a*cos(pi/2-alpha),'magnitude',H);

SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',6);

A1 = AQUINAS_Analysis_Object('type','LA');

O1 = AQUINAS_Output_Mode_Object('mode','VarArgOut');

LA_out = AQUINAS_Protocol(M1,S1,C1,C2,F1,SOL1,A1,O1);


%%%%%%%%%%%%%%%%%%%%%%%
% Analytical solution %
%%%%%%%%%%%%%%%%%%%%%%%
% The following bending theory solution can be found in chapter 6.2.1 of [1],
% while the example itself can be found in pages 332-333 of [1].

% Please note a difference in notation below. In AQUINAS, 'u' and 'w' are the
% global displacements parallel to the r and z axes respectively. According to
% Flugge's notation in [1], however, 'w' is the displacement normal to the shell's
% wall, while 'v' is the meridional displacement. Consequently, in the below the
% meaning of 'w' and 'v' is the OPPOSITE of that used in AQUINAS for 'w' and 'u'.

numstan = 100*ceil(alpha*180/pi) + 1; % Number of Analytical Stations, necessary for desired accuracy of displacement field

% Membrane theory solution
Nphi_m = @(phi) -H*((sin(pi/2-alpha))^2)*(sin(alpha))./((sin(phi)).^2); % Meridional membrane stress resultant
Ntheta_m = @(phi) -Nphi_m(phi); % Circumferential membrane stress resultant
EpsPhi_m = @(phi) (Nphi_m(phi)-nu*Ntheta_m(phi))/(E*t); % Meridional midsurface strain for membrane stress state only
EpsTheta_m = @(phi) (Ntheta_m(phi)-nu*Nphi_m(phi))/(E*t); % Circumferential midsurface strain for membrane stress state only
SigmaPhi_m = @(phi) Nphi_m(phi)/t; % Meridional membrane stress
SigmaTheta_m = @(phi) Ntheta_m(phi)/t; % Circumferential membrane stress

% Bending Theory solution
kappa4 = 3*(1-(nu^2))*((a^2)/(t^2)) - (nu^2)/4;  kappa2 = sqrt(kappa4);  kappa = sqrt(kappa2); % kappa value, eq. 6.20 of [1]
D = E*t/(1-nu^2); % Membrane stiffness
K = (t^2)*D/12; % Bending stiffness

% Hypergeometric series solution and derivatives, eq. 6.25 of [1]
z2 = @(phi) cos(phi) + 5*((cos(phi)).^3)/6 + (95-4*kappa4)*((cos(phi)).^5)/120 + (3895-260*kappa4)*((cos(phi)).^7)/5760;
z4 = @(phi) (kappa2/3)*((cos(phi).^3) + (6/5)*(cos(phi).^5) + (1079-4*kappa4)*((cos(phi)).^7)/840 + (10063-68*kappa4)*((cos(phi)).^9)/7560);
z2p = @(phi) -sin(phi).*(1 + 5*((cos(phi)).^2)/2 + (95-4*kappa4)*((cos(phi)).^4)/24 + (3895-260*kappa4)*((cos(phi)).^6)/720);
z4p = @(phi) -sin(phi)*(kappa2/3).*(3*(cos(phi).^2) + (6*(cos(phi).^4)) + (1079-4*kappa4)*((cos(phi)).^7)/120 + (10063-68*kappa4)*((cos(phi)).^8)/840);

phi1 = pi/2-alpha; % Complementary angle to alpha, see Fig. 6.3 of [1]

% Factors for the top edge boundary condition, needed for the evaluation of C2 and C4 constants (only those two needed due to symmetry of the problem)
BC1C2 = z2(phi1)*sin(phi1);
BC1C4 = z4(phi1)*sin(phi1);
BC2C2 = ((1+nu)*sin(phi1)*(2*kappa2*z4p(phi1) - nu*z2p(phi1))+cos(phi1)*(2*kappa2*z4(phi1) - nu*z2(phi1)))*K/(D*a*(1-nu));
BC2C4 = -((1+nu)*sin(phi1)*(2*kappa2*z2p (phi1) + nu*z4p(phi1))+cos(phi1)*(2*kappa2*z2(phi1) + nu*z4(phi1)))*K/(D*a*(1-nu));

% Forming and solving matrix system to find C2 and C4
BCfactors = [BC1C2   BC1C4
             BC2C2   BC2C4];
V = [H*cos(alpha) 0]';
C = BCfactors\V; C2 = C(1); C4 = C(2);


Nphi_b = @(phi) -cos(phi).*(C2*z2(phi)+C4*z4(phi)); % Meridional membrane stress resultant, according to eq. 6.26b
Ntheta_b = @(phi) -sin(phi).*(C2*z2p(phi)+C4*z4p(phi)) - cos(phi).*(C2*z2(phi)+C4*z4(phi)); % Circumferential membrane stress resultant, according to eq. 6.26b
Mphi_b = @(phi) (K/(D*(1-nu^2)*a))*sin(phi).*(C2*(2*kappa2*z4p(phi)-nu*z2p(phi))-C4*(2*kappa2*z2p(phi)+nu*z4p(phi))) + (K/(D*(1-nu)*a))*cos(phi).*(C2*(2*kappa2*z4(phi)-nu*z2(phi))-C4*(2*kappa2*z2(phi)+nu*z4(phi))); % Meridional bending stress resultant, according to eq. 6.26e
Mtheta_b = @(phi) (K*nu/(D*(1-nu^2)*a))*sin(phi).*(C2*(2*kappa2*z4p(phi)-nu*z2p(phi))-C4*(2*kappa2*z2p(phi)-nu*z4p(phi))) + (K/(D*(1-nu)*a))*cos(phi).*(C2*(2*kappa2*z4(phi)-nu*z2(phi))-C4*(2*kappa2*z2(phi)+nu*z4(phi)));  % Circumferential bending stress resultant, according to eq. 6.26e
% The following stresses-strains are obtained from general shell theory relations
SigmaPhi_bi = @(phi) Nphi_b(phi)/t + 6.0*Mphi_b(phi)/(t*t); % Complete meridional (axial) membrane stress - inner surface
SigmaPhi_bm = @(phi) Nphi_b(phi)/t; % Complete meridional (axial) membrane stress - midfsurface
SigmaPhi_bo = @(phi) Nphi_b(phi)/t - 6.0*Mphi_b(phi)/(t*t); % Complete meridional (axial) membrane stress - outer surface
SigmaTheta_bi = @(phi) Ntheta_b(phi)/t + 6.0*Mtheta_b(phi)/(t*t); % Complete circumferential membrane stress - inner surface
SigmaTheta_bm = @(phi) Ntheta_b(phi)/t; % Complete circumferential membrane stress - midsurface
SigmaTheta_bo = @(phi) Ntheta_b(phi)/t - 6.0*Mtheta_b(phi)/(t*t); % Complete circumferential membrane stress - outer surface
KappaPhi = @(phi) (Mphi_b(phi)-nu*Mtheta_b(phi))/(K*(1-nu^2)); % Meridional curvature
KappaTheta = @(phi) (Mtheta_b(phi)-nu*Mphi_b(phi))/(K*(1-nu^2)); % Circumferential curvature
EpsPhi_bi = @(phi) ((Nphi_b(phi)-nu*Ntheta_b(phi))/(D*(1-nu^2))) + 0.5*t*KappaPhi(phi); % Complete meridional (axial) strain - inner surface
EpsPhi_bm = @(phi) (Nphi_b(phi)-nu*Ntheta_b(phi))/(D*(1-nu^2)); % Complete meridional (axial) strain - midsurface
EpsPhi_bo = @(phi) ((Nphi_b(phi)-nu*Ntheta_b(phi))/(D*(1-nu^2))) - 0.5*t*KappaPhi(phi); % Complete meridional (axial) strain - outer surface
EpsTheta_bi = @(phi) ((Ntheta_b(phi)-nu*Nphi_b(phi))/(D*(1-nu^2))) + 0.5*t*KappaTheta(phi); % Complete circumferential strain - inner surface
EpsTheta_bm = @(phi) (Ntheta_b(phi)-nu*Nphi_b(phi))/(D*(1-nu^2)); % Complete circumferential strain - midsurface
EpsTheta_bo = @(phi) ((Ntheta_b(phi)-nu*Nphi_b(phi))/(D*(1-nu^2))) - 0.5*t*KappaTheta(phi); % Complete circumferential strain - outer surface

% Numerical scheme for the calculation of the midsurface displacements,
% that capitilizes on the boundary conditions at phi=pi/2 (v = 0 there)
v_prev = 0.0;
dvds = 0.0;
vg = zeros(numstan,1);  wg = vg;
for i=1:numstan
    iphi = pi/2 + (i-1)*alpha/(numstan-1);
    vbn = v_prev;
    wbn = a*EpsTheta_bm(iphi) - v_prev*cos(iphi)/sin(iphi); % evaluation of normal displacement at current numerical station, according to eqs. 6.4 of [1], simplified for spherical shell case (r1=r2=r)
    dvds = EpsPhi_bm(iphi) - wbn/a; % evaluation of meridional displacement derivative at current numerical station, according to eqs. 6.4 of [1], simplified for spherical shell case (r1=r2=r)
    v_prev = vbn + dvds*a*alpha/(numstan-1); % numerical integration of meridional displacements, based on the finite difference method, assuming a finite length of ds=r*alpha/(numstan-1)
    vg(i) = cos(iphi)*vbn + sin(iphi)*wbn; % transformation for the evaluation of global radial displacement v (local to global transformation), for comparison with AQUINAS
    wg(i) = -(cos(iphi)*wbn - sin(iphi)*vbn); % transformation for the evaluation of global axial displacement w (local to global transformation), for comparison with AQUINAS
end
wg = flipud(wg); vg = flipud(vg);

%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%
%% Radial displacement plot
figure('Name','Radial displacement plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(2,1);

% % Shell bending theory solution
plot(vg, 1:-(1/(numstan-1)):0,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
lgnd{1} = 'Bending theory displacement';

% AQUINAS solution
plot(LA_out.DOFs.u{1}, 1:-(1/numelem):0,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',MS);
plot(LA_out.DOFs.u{1}, 1:-(1/numelem):0,...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','w','MarkerSize',0.8*MS);
lgnd{2} = 'AQUINAS FE solution';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Radial displacement of the wall [$mm$]','interpreter','latex','fontsize',FS);
ylabel('Normalized phi angle from midpoint of spherical zone [-]','interpreter','latex','fontsize',FS);
set(gca, 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Axial displacement plot
figure('Name','Axial displacement plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(2,1);

% Shell bending theory solution
plot(wg, 1:-(1/(numstan-1)):0,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
lgnd{1} = 'Bending theory displacement';

% AQUINAS solution
plot(LA_out.DOFs.w{1}, 1:-(1/numelem):0,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',MS);
plot(LA_out.DOFs.w{1}, 1:-(1/numelem):0,...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','w','MarkerSize',0.8*MS);
lgnd{2} = 'AQUINAS FE solution';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Axial displacement of the wall [$mm$]','interpreter','latex','fontsize',FS);
ylabel('Normalized phi angle from midpoint of spherical zone [-]','interpreter','latex','fontsize',FS);
set(gca, 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Meridional strain plot
figure('Name','Meridional strain plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(7,1);

% Shell membrane theory solution
plot(EpsPhi_m((pi/2-alpha):(alpha/(numstan-1)):pi/2), 1:-(1/(numstan-1)):0,...
    'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',LW); grid on; hold on;
lgnd{1} = 'Membrane theory midsurface strain';

% Shell bending theory solution
plot(EpsPhi_bi((pi/2-alpha):(alpha/(numstan-1)):pi/2), 1:-(1/(numstan-1)):0,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(EpsPhi_bm((pi/2-alpha):(alpha/(numstan-1)):pi/2), 1:-(1/(numstan-1)):0,...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(EpsPhi_bo((pi/2-alpha):(alpha/(numstan-1)):pi/2), 1:-(1/(numstan-1)):0,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{2} = 'Bending theory strain (inner surface)';
lgnd{3} = 'Bending theory strain (midsurface)';
lgnd{4} = 'Bending theory strain (outer surface)';

% AQUINAS solution
plot(LA_out.Strains.Epsilon.Eps_phi.i{1}, 1:-(1/numelem):0,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_phi.m{1}, 1:-(1/numelem):0,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_phi.o{1}, 1:-(1/numelem):0,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_phi.i{1}, 1:-(1/numelem):0,...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Strains.Epsilon.Eps_phi.m{1}, 1:-(1/numelem):0,...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(LA_out.Strains.Epsilon.Eps_phi.o{1}, 1:-(1/numelem):0,...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{5} = 'AQUINAS FE solution (inner surface)';
lgnd{6} = 'AQUINAS FE solution (midsurface)';
lgnd{7} = 'AQUINAS FE solution (outer surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Meridional strains in the wall [$-$]','interpreter','latex','fontsize',FS);
ylabel('Normalized phi angle from midpoint of spherical zone [-]','interpreter','latex','fontsize',FS);
set(gca, 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Circumferential strain plot
figure('Name','Circumferential strain plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(7,1);

% Shell membrane theory solution
plot(EpsTheta_m((pi/2-alpha):(alpha/(numstan-1)):pi/2), 1:-(1/(numstan-1)):0,...
    'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',LW); grid on; hold on;
lgnd{1} = 'Membrane theory midsurface strain';

% Shell bending theory solution
plot(EpsTheta_bi((pi/2-alpha):(alpha/(numstan-1)):pi/2), 1:-(1/(numstan-1)):0,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(EpsTheta_bm((pi/2-alpha):(alpha/(numstan-1)):pi/2), 1:-(1/(numstan-1)):0,...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(EpsTheta_bo((pi/2-alpha):(alpha/(numstan-1)):pi/2), 1:-(1/(numstan-1)):0,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{2} = 'Bending theory strain (inner surface)';
lgnd{3} = 'Bending theory strain (midsurface)';
lgnd{4} = 'Bending theory strain (outer surface)';

% AQUINAS solution
plot(LA_out.Strains.Epsilon.Eps_theta.i{1}, 1:-(1/numelem):0,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_theta.m{1}, 1:-(1/numelem):0,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_theta.o{1}, 1:-(1/numelem):0,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_theta.i{1}, 1:-(1/numelem):0,...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Strains.Epsilon.Eps_theta.m{1}, 1:-(1/numelem):0,...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(LA_out.Strains.Epsilon.Eps_theta.o{1}, 1:-(1/numelem):0,...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{5} = 'AQUINAS FE solution (inner surface)';
lgnd{6} = 'AQUINAS FE solution (midsurface)';
lgnd{7} = 'AQUINAS FE solution (outer surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Circumferential strains in the wall [$-$]','interpreter','latex','fontsize',FS);
ylabel('Normalized phi angle from midpoint of spherical zone [-]','interpreter','latex','fontsize',FS);
set(gca, 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Meridional stress plot
figure('Name','Meridional stress plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(7,1);

% Shell membrane theory solution
plot(SigmaPhi_m((pi/2-alpha):(alpha/(numstan-1)):pi/2), 1:-(1/(numstan-1)):0,...
    'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',LW); grid on; hold on;
lgnd{1} = 'Membrane theory midsurface stress';

% Shell bending theory solution
plot(SigmaPhi_bi((pi/2-alpha):(alpha/(numstan-1)):pi/2), 1:-(1/(numstan-1)):0,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(SigmaPhi_bm((pi/2-alpha):(alpha/(numstan-1)):pi/2), 1:-(1/(numstan-1)):0,...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(SigmaPhi_bo((pi/2-alpha):(alpha/(numstan-1)):pi/2), 1:-(1/(numstan-1)):0,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{2} = 'Bending theory stress (inner surface)';
lgnd{3} = 'Bending theory stress (midsurface)';
lgnd{4} = 'Bending theory stress (outer surface)';

% AQUINAS solution
plot(LA_out.Stresses.Sig_phi.i{1}, 1:-(1/numelem):0,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stresses.Sig_phi.m{1}, 1:-(1/numelem):0,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(LA_out.Stresses.Sig_phi.o{1}, 1:-(1/numelem):0,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stresses.Sig_phi.i{1}, 1:-(1/numelem):0,...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Stresses.Sig_phi.m{1}, 1:-(1/numelem):0,...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(LA_out.Stresses.Sig_phi.o{1}, 1:-(1/numelem):0,...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{5} = 'AQUINAS FE solution (inner surface)';
lgnd{6} = 'AQUINAS FE solution (midsurface)';
lgnd{7} = 'AQUINAS FE solution (outer surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Meridional stresses in the wall [$N/mm^2$]','interpreter','latex','fontsize',FS);
ylabel('Normalized phi angle from midpoint of spherical zone [-]','interpreter','latex','fontsize',FS);
set(gca, 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Circumferential stress plot
figure('Name','Circumferential stress plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(7,1);

% Shell membrane theory solution
plot(SigmaTheta_m((pi/2-alpha):(alpha/(numstan-1)):pi/2), 1:-(1/(numstan-1)):0,...
    'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',LW); grid on; hold on;
lgnd{1} = 'Membrane theory midsurface stress';

% Shell bending theory solution
plot(SigmaTheta_bi((pi/2-alpha):(alpha/(numstan-1)):pi/2), 1:-(1/(numstan-1)):0,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(SigmaTheta_bm((pi/2-alpha):(alpha/(numstan-1)):pi/2), 1:-(1/(numstan-1)):0,...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(SigmaTheta_bo((pi/2-alpha):(alpha/(numstan-1)):pi/2), 1:-(1/(numstan-1)):0,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{2} = 'Bending theory stress (inner surface)';
lgnd{3} = 'Bending theory stress (midsurface)';
lgnd{4} = 'Bending theory stress (outer surface)';

% AQUINAS solution
plot(LA_out.Stresses.Sig_theta.i{1}, 1:-(1/numelem):0,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stresses.Sig_theta.m{1}, 1:-(1/numelem):0,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(LA_out.Stresses.Sig_theta.o{1}, 1:-(1/numelem):0,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stresses.Sig_theta.i{1}, 1:-(1/numelem):0,...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Stresses.Sig_theta.m{1}, 1:-(1/numelem):0,...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(LA_out.Stresses.Sig_theta.o{1}, 1:-(1/numelem):0,...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{5} = 'AQUINAS FE solution (inner surface)';
lgnd{6} = 'AQUINAS FE solution (midsurface)';
lgnd{7} = 'AQUINAS FE solution (outer surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Circumferential stresses in the wall [$N/mm^2$]','interpreter','latex','fontsize',FS);
ylabel('Normalized phi angle from midpoint of spherical zone [-]','interpreter','latex','fontsize',FS);
set(gca, 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Membrane stress resultants
figure('Name','Membrane stress resultant plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(6,1);

% Shell membrane theory solution
plot(Nphi_m((pi/2-alpha):(alpha/(numstan-1)):pi/2), 1:-(1/(numstan-1)):0,...
    'Color',[1.0 0.7 0.7],'LineStyle','--','LineWidth',LW); grid on; hold on;
plot(Ntheta_m((pi/2-alpha):(alpha/(numstan-1)):pi/2), 1:-(1/(numstan-1)):0,...
    'Color',[0.7 0.7 1.0],'LineStyle','--','LineWidth',LW);
lgnd{1} = 'Membrane theory meridional MSR';
lgnd{2} = 'Membrane theory circumferential MSR';

% Shell bending theory solution
plot(Nphi_b((pi/2-alpha):(alpha/(numstan-1)):pi/2), 1:-(1/(numstan-1)):0,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); hold on;
plot(Ntheta_b((pi/2-alpha):(alpha/(numstan-1)):pi/2), 1:-(1/(numstan-1)):0,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); hold on;
lgnd{3} = 'Bending theory meridional MSR';
lgnd{4} = 'Bending theory circumferential MSR';

% AQUINAS solution
plot(LA_out.Stress_Reslts.Membrane.Nphi{1}, 1:-(1/numelem):0,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stress_Reslts.Membrane.Ntheta{1}, 1:-(1/numelem):0,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stress_Reslts.Membrane.Nphi{1}, 1:-(1/numelem):0,...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Stress_Reslts.Membrane.Ntheta{1}, 1:-(1/numelem):0,...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{5} = 'AQUINAS FE solution meridional MSR';
lgnd{6} = 'AQUINAS FE solution circumferential MSR';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Membrane stress resultants (MSR) in the wall [$N/mm$]','interpreter','latex','fontsize',FS);
ylabel('Normalized phi angle from midpoint of spherical zone [-]','interpreter','latex','fontsize',FS);
set(gca, 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Bending stress resultants
figure('Name','Bending stress resultant plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(4,1);

% Shell bending theory solution
plot(Mphi_b((pi/2-alpha):(alpha/(numstan-1)):pi/2), 1:-(1/(numstan-1)):0,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
plot(Mtheta_b((pi/2-alpha):(alpha/(numstan-1)):pi/2), 1:-(1/(numstan-1)):0,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{1} = 'Bending theory meridional BSR';
lgnd{2} = 'Bending theory circumferential BSR';

% AQUINAS solution
plot(-LA_out.Stress_Reslts.Bending.Mphi{1}, 1:-(1/numelem):0,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(-LA_out.Stress_Reslts.Bending.Mtheta{1}, 1:-(1/numelem):0,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(-LA_out.Stress_Reslts.Bending.Mphi{1}, 1:-(1/numelem):0,...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(-LA_out.Stress_Reslts.Bending.Mtheta{1}, 1:-(1/numelem):0,...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{3} = 'AQUINAS FE solution meridional BSR';
lgnd{4} = 'AQUINAS FE solution circumferential BSR';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Bending stress resultants (BSR) in the wall [$Nmm/mm$]','interpreter','latex','fontsize',FS);
ylabel('Normalized phi angle from midpoint of spherical zone [-]','interpreter','latex','fontsize',FS);
set(gca, 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


% BSD 3-Clause License
%  
% Copyright (c) 2023, Mr Achilleas Filippidis and Dr Adam Jan Sadowski of 
% Imperial College London. All rights reserved.
%  
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%  
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
%  
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
%  
% 3. Neither the name of the copyright holder, nor of Imperial College, nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
%  
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%  
% THE USE OF THE MATLAB LANGUAGE DOES NOT IMPLY ENDORSEMENT BY MATHWORKS.