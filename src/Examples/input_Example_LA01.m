%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example LA 01: LA of a thin cylindrical shell under combined uniform axial
% compression and uniform internal pressure. A BC1r / C1 condition is
% assumed at the base, and a BC2f / S3 condition at the top.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% [1] Sadowski A.J. & Rotter J.M. (2012) "Cylindrical shell bending theory
% for orthotropic shells under general axisymmetric pressure distributions"
% Engineering Structures, 42, 258-265.
% http://dx.doi.org/10.1016/j.engstruct.2012.04.024
% N.B. the analytical solution to this classical problem may be found in
% several classical texts (Timoshenko, Flugge etc.) but the above
% reference is possibly the most accessible.
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
% Geometry
r = 500.0; % [mm] - Cylinder radius
L = 500.0; % [mm] - Cylinder length / height
t = 5.0; % [mm] - Cylinder thickness

% Material
MatName = 'Steel'; % Material name
E = 2.0e5; % [N/mm2] - Young's modulus of steel
nu = 0.3; % [-] - Poisson ratio of steel

% Loading
N = 1000.0; % [N/mm] - Distributed axisymmetric axial load along top edge (acting vertically downwards)
pn = 1.0; % [N/mm2] - Distributed normal pressure along inner surface (acting radially outwards)

% Plot controls
FS = 20; % [-] - Font Size (pts)
FSL = 18; % [-] - Font Size for Legends (pts)
MS = 3; % [-] - Markser size (pts)
LW = 3; % [-] - Line Width (pts)


%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%
M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)));

S1 = AQUINAS_Segment_Object('type','Cone','material',MatName,'thickness',t,'rzbot',[r 0.0],'rztop',[r L],'els',100,'meshtype','E','g',0.5);

% BC1r condition at base of cylinder
C1 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'u'});
C2 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'w'});
C3 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'b'});

% BC2f condition at top of cylinder
C4 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',L,'dofs',{'u'});

F1 = AQUINAS_Nodal_Force_Object('type','Fw','rcoord',r,'zcoord',L,'magnitude',-N);

P1 = AQUINAS_Distributed_Pressure_Object('segments',{S1},'type','pn','functionHandle',@() pn);

SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',6);

A1 = AQUINAS_Analysis_Object('type','LA');

O1 = AQUINAS_Output_Mode_Object('mode','VarArgOut');

LA_out = AQUINAS_Protocol(M1,S1,C1,C2,C3,C4,F1,P1,SOL1,A1,O1);


%%%%%%%%%%%%%%%%%%%%%%%
% Analytical solution %
%%%%%%%%%%%%%%%%%%%%%%%
% Please note a difference in notation below. In AQUINAS, u and w are the
% global displacements parallel to the r and z axes respectively. In usual
% notation used in classical cylindrical shell theory, however, w is the
% displacement normal to the wall and u is the axial displacement.
% Consequently, in the below the meaning of 'w' and 'u' is the OPPOSITE of
% that used in AQUINAS.

% Membrane theory solution
NZm = @(z) -N*ones(size(z)); % Meridional (axial) membrane stress resultant
NThm = @(z) pn*r*ones(size(z)); % Circumferential membrane stress resultant
wm = @(z) (r/(E*t))*(NThm(z) - nu*NZm(z)); % Radial midsurface displacement for membrane stress state only
wm_p = @(z) zeros(size(z)); % First derivative of wm w.r.t. z
wm_pp = @(z) zeros(size(z)); % Second derivative of wm w.r.t. z
um = @(z) -z/(E*t).*(nu*NThm(z) - NZm(z)); % Meridional (axial) midsurface displacement for membrane stress state only
um_p = @(z) -(nu*NThm(z) - NZm(z))/(E*t); % First derivative of um w.r.t. z
EpsZ_m = @(z) um_p(z); % Meridional (axial) midsurface strain for membrane stress state only
EpsTh_m = @(z) wm(z)/r; % Circumferential midsurface strain for membrane stress state only
SigmaZ_m = @(z) NZm(z)/t; % Meridional (axial) membrane stress
SigmaTh_m = @(z) NThm(z)/t; % Circumferential membrane stress

% Linear bending half-wavelength
lambda = pi*sqrt(r*t)*(3.0*(1.0 - nu*nu))^(-0.25);

% Coefficients (their derivatives and anti-derivatives) of the radial
% displacement field
c1 = @(z) exp(-pi*z/lambda).*cos(pi*z/lambda);
c2 = @(z) exp(-pi*z/lambda).*sin(pi*z/lambda);
c3 = @(z) exp(pi*z/lambda).*cos(pi*z/lambda);
c4 = @(z) exp(pi*z/lambda).*sin(pi*z/lambda);
c1p = @(z) -(pi/lambda)*exp(-pi*z/lambda).*(cos(pi*z/lambda)+sin(pi*z/lambda));
c2p = @(z) (pi/lambda)*exp(-pi*z/lambda).*(cos(pi*z/lambda)-sin(pi*z/lambda));
c3p = @(z) (pi/lambda)*exp(pi*z/lambda).*(cos(pi*z/lambda)-sin(pi*z/lambda));
c4p = @(z) (pi/lambda)*exp(pi*z/lambda).*(cos(pi*z/lambda)+sin(pi*z/lambda));
c1pp = @(z) 2.0*(pi/lambda)^2*exp(-pi*z/lambda).*sin(pi*z/lambda);
c2pp = @(z) -2.0*(pi/lambda)^2*exp(-pi*z/lambda).*cos(pi*z/lambda);
c3pp = @(z) -2.0*(pi/lambda)^2*exp(pi*z/lambda).*sin(pi*z/lambda);
c4pp = @(z) 2.0*(pi/lambda)^2*exp(pi*z/lambda).*cos(pi*z/lambda);
c1ppp = @(z) 2.0*(pi/lambda)^3*exp(-pi*z/lambda).*(cos(pi*z/lambda)-sin(pi*z/lambda));
c2ppp = @(z) 2.0*(pi/lambda)^3*exp(-pi*z/lambda).*(cos(pi*z/lambda)+sin(pi*z/lambda));
c3ppp = @(z) -2.0*(pi/lambda)^3*exp(pi*z/lambda).*(cos(pi*z/lambda)+sin(pi*z/lambda));
c4ppp = @(z) 2.0*(pi/lambda)^3*exp(pi*z/lambda).*(cos(pi*z/lambda)-sin(pi*z/lambda));
c1i = @(z) 0.5*lambda/pi*(1-exp(-pi*z/lambda).*(cos(pi*z/lambda)-sin(pi*z/lambda)));
c2i = @(z) 0.5*lambda/pi*(1-exp(-pi*z/lambda).*(cos(pi*z/lambda)+sin(pi*z/lambda)));
c3i = @(z) 0.5*lambda/pi*(-1+exp(pi*z/lambda).*(sin(pi*z/lambda)+cos(pi*z/lambda)));
c4i = @(z) 0.5*lambda/pi*(1+exp(pi*z/lambda).*(sin(pi*z/lambda)-cos(pi*z/lambda)));

% Construction and solution of matrix system
a1_0 = c1(0.0); a2_0 = c2(0.0); a3_0 = c3(0.0); a4_0 = c4(0.0);
a1p_0 = c1p(0.0); a2p_0 = c2p(0.0); a3p_0 = c3p(0.0); a4p_0 = c4p(0.0);
a1_h = c1(L); a2_h = c2(L); a3_h = c3(L); a4_h = c4(L);
a1pp_h = c1pp(L); a2pp_h = c2pp(L); a3pp_h = c3pp(L); a4pp_h = c4pp(L);

% Forming matrix system
M = [a1_0    a2_0    a3_0    a4_0;...
     a1p_0   a2p_0   a3p_0   a4p_0;...
     a1_h    a2_h    a3_h    a4_h;...
     a1pp_h  a2pp_h  a3pp_h  a4pp_h];
V = -[wm(0) wm_p(0) wm(L) wm_pp(L)]';
disp('Ignore any warning relating to matrix ill-conditioning - it does not relate to AQUINAS.');
As = M\V; A1 = As(1); A2 = As(2); A3 = As(3); A4 = As(4);

% Bending theory solution
C = E*t/(1.0-nu*nu); % Membrane stiffnesss
D = E*t*t*t/(12.0*(1.0 - nu*nu)); % Bending stiffness
wb = @(z) A1*c1(z) + A2*c2(z) + A3*c3(z) + A4*c4(z) + wm(z); % Complete radial displacement
wb_pp = @(z) A1*c1pp(z) + A2*c2pp(z) + A3*c3pp(z) + A4*c4pp(z) + wm_pp(z); % Second derivative of wb w.r.t. z
ub = @(z) -z*N/C - nu/r*(A1*c1i(z) + A2*c2i(z) + A3*c3i(z) + A4*c4i(z) + z.*wm(z)); % Complete meridional (axial) displacement
ub_p = @(z) -N/C - nu/r*wb(z); % First derivative of ub w.r.t. z
KappaZ = @(z) wb_pp(z); % Meridional curvature
EpsZ_bi = @(z) ub_p(z) + 0.5*t*KappaZ(z); % Complete meridional (axial) strain - inner surface
EpsZ_bm = @(z) ub_p(z); % Complete meridional (axial) strain - midsurface
EpsZ_bo = @(z) ub_p(z) - 0.5*t*KappaZ(z); % Complete meridional (axial) strain - outer surface
EpsTh_bi = @(z) wb(z)/r; % Complete circumferential strain - inner surface
EpsTh_bm = @(z) wb(z)/r; % Complete circumferential strain - midsurface
EpsTh_bo = @(z) wb(z)/r; % Complete circumferential strain - outer surface
NZ_b = @(z) C*(ub_p(z) + nu/r*wb(z)); % Complete meridional (axial) membrane stress resultant
NTh_b = @(z) C*(wb(z)/r + nu*ub_p(z)); % Complete circumferential membrane stress resultant
MZ_b = @(z) D*wb_pp(z); % Meridional bending moment stress resultant
MTh_b = @(z) nu*MZ_b(z); % Circumferential bending moment stress resultant
SigmaZ_bi = @(z) NZ_b(z)/t + 6.0*MZ_b(z)/(t*t); % Complete meridional (axial) membrane stress - inner surface
SigmaZ_bm = @(z) NZ_b(z)/t; % Complete meridional (axial) membrane stress - midfsurface
SigmaZ_bo = @(z) NZ_b(z)/t - 6.0*MZ_b(z)/(t*t); % Complete meridional (axial) membrane stress - outer surface
SigmaTh_bi = @(z) NTh_b(z)/t + 6.0*MTh_b(z)/(t*t); % Complete circumferential membrane stress - inner surface
SigmaTh_bm = @(z) NTh_b(z)/t; % Complete circumferential membrane stress - midsurface
SigmaTh_bo = @(z) NTh_b(z)/t - 6.0*MTh_b(z)/(t*t); % Complete circumferential membrane stress - outer surface


%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%
%% Radial displacement plot
figure('Name','Radial displacement plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(3,1);

% Shell membrane theory solution
plot(wm(0:L/1000:L), 0:L/1000:L,...
    'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',LW); grid on; hold on;
lgnd{1} = 'Membrane theory displacement';

% Shell bending theory solution
plot(wb(0:L/1000:L), 0:L/1000:L,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{2} = 'Bending theory displacement';

% AQUINAS solution
plot(LA_out.DOFs.u{1}, LA_out.Shell_Geom.z{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',MS);
plot(LA_out.DOFs.u{1}, LA_out.Shell_Geom.z{1},...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','w','MarkerSize',0.8*MS);
lgnd{3} = 'AQUINAS FE solution';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Radial displacement of the wall [$mm$]','interpreter','latex','fontsize',FS);
ylabel('Cylinder vertical coordinate [$mm$]','interpreter','latex','fontsize',FS);
set(gca,'ticklabelinterpreter','latex','fontsize',FS,'tickdir','out');


%% Axial / meridional displacement plot
figure('Name','Axial displacement plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(3,1);

% Shell membrane theory solution
plot(um(0:L/1000:L), 0:L/1000:L,...
    'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',LW); grid on; hold on;
lgnd{1} = 'Membrane theory displacement';

% Shell bending theory solution
plot(ub(0:L/1000:L), 0:L/1000:L,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{2} = 'Bending theory displacement';

% AQUINAS solution
plot(LA_out.DOFs.w{1}, LA_out.Shell_Geom.z{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',MS);
plot(LA_out.DOFs.w{1}, LA_out.Shell_Geom.z{1},...
    'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor','none','MarkerSize',0.8*MS);
lgnd{3} = 'AQUINAS FE solution';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Axial / meridional displacement of the wall [$mm$]','interpreter','latex','fontsize',FS);
ylabel('Cylinder vertical coordinate [$mm$]','interpreter','latex','fontsize',FS);
set(gca,'ticklabelinterpreter','latex','fontsize',FS,'tickdir','out');


%% Axial / meridional strain plot
figure('Name','Axial strain plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(7,1);

% Shell membrane theory solution
plot(EpsZ_m(0:L/1000:L), 0:L/1000:L,...
    'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',LW); grid on; hold on;
lgnd{1} = 'Membrane theory midsurface strain';

% Shell bending theory solution
plot(EpsZ_bi(0:L/1000:L), 0:L/1000:L,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(EpsZ_bm(0:L/1000:L), 0:L/1000:L,...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(EpsZ_bo(0:L/1000:L), 0:L/1000:L,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{2} = 'Bending theory strain (inner surface)';
lgnd{3} = 'Bending theory strain (midsurface)';
lgnd{4} = 'Bending theory strain (outer surface)';

% AQUINAS solution
plot(LA_out.Strains.Epsilon.Eps_phi.i{1}, LA_out.Shell_Geom.z{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_phi.m{1}, LA_out.Shell_Geom.z{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_phi.o{1}, LA_out.Shell_Geom.z{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_phi.i{1}, LA_out.Shell_Geom.z{1},...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Strains.Epsilon.Eps_phi.m{1}, LA_out.Shell_Geom.z{1},...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(LA_out.Strains.Epsilon.Eps_phi.o{1}, LA_out.Shell_Geom.z{1},...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{5} = 'AQUINAS FE solution (inner surface)';
lgnd{6} = 'AQUINAS FE solution (midsurface)';
lgnd{7} = 'AQUINAS FE solution (outer surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Axial / meridional strains in the wall [$-$]','interpreter','latex','fontsize',FS);
ylabel('Cylinder vertical coordinate [$mm$]','interpreter','latex','fontsize',FS);
set(gca,'ticklabelinterpreter','latex','fontsize',FS,'tickdir','out');


%% Circumferential strain plot
figure('Name','Circumferential strain plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(7,1);

% Shell membrane theory solution
plot(EpsTh_m(0:L/1000:L), 0:L/1000:L,...
    'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',LW); grid on; hold on;
lgnd{1} = 'Membrane theory midsurface strain';

% Shell bending theory solution
plot(EpsTh_bi(0:L/1000:L), 0:L/1000:L,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(EpsTh_bm(0:L/1000:L), 0:L/1000:L,...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(EpsTh_bo(0:L/1000:L), 0:L/1000:L,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{2} = 'Bending theory strain (inner surface)';
lgnd{3} = 'Bending theory strain (midsurface)';
lgnd{4} = 'Bending theory strain (outer surface)';

% AQUINAS solution
plot(LA_out.Strains.Epsilon.Eps_theta.i{1}, LA_out.Shell_Geom.z{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_theta.m{1}, LA_out.Shell_Geom.z{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_theta.o{1}, LA_out.Shell_Geom.z{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_theta.i{1}, LA_out.Shell_Geom.z{1},...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Strains.Epsilon.Eps_theta.m{1}, LA_out.Shell_Geom.z{1},...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(LA_out.Strains.Epsilon.Eps_theta.o{1}, LA_out.Shell_Geom.z{1},...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{5} = 'AQUINAS FE solution (inner surface)';
lgnd{6} = 'AQUINAS FE solution (midsurface)';
lgnd{7} = 'AQUINAS FE solution (outer surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Circumferential strains in the wall [$-$]','interpreter','latex','fontsize',FS);
ylabel('Cylinder vertical coordinate [$mm$]','interpreter','latex','fontsize',FS);
set(gca,'ticklabelinterpreter','latex','fontsize',FS,'tickdir','out');


%% Axial / meridional stress plot
figure('Name','Axial stress plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(7,1);

% Shell membrane theory solution
plot(SigmaZ_m(0:L/1000:L), 0:L/1000:L,...
    'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',LW); grid on; hold on;
lgnd{1} = 'Membrane theory midsurface stress';

% Shell bending theory solution
plot(SigmaZ_bi(0:L/1000:L), 0:L/1000:L,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(SigmaZ_bm(0:L/1000:L), 0:L/1000:L,...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(SigmaZ_bo(0:L/1000:L), 0:L/1000:L,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{2} = 'Bending theory stress (inner surface)';
lgnd{3} = 'Bending theory stress (midsurface)';
lgnd{4} = 'Bending theory stress (outer surface)';

% AQUINAS solution
plot(LA_out.Stresses.Sig_phi.i{1}, LA_out.Shell_Geom.z{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stresses.Sig_phi.m{1}, LA_out.Shell_Geom.z{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(LA_out.Stresses.Sig_phi.o{1}, LA_out.Shell_Geom.z{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stresses.Sig_phi.i{1}, LA_out.Shell_Geom.z{1},...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Stresses.Sig_phi.m{1}, LA_out.Shell_Geom.z{1},...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(LA_out.Stresses.Sig_phi.o{1}, LA_out.Shell_Geom.z{1},...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{5} = 'AQUINAS FE solution (inner surface)';
lgnd{6} = 'AQUINAS FE solution (midsurface)';
lgnd{7} = 'AQUINAS FE solution (outer surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Axial / meridional stresses in the wall [$N/mm^2$]','interpreter','latex','fontsize',FS);
ylabel('Cylinder vertical coordinate [$mm$]','interpreter','latex','fontsize',FS);
set(gca,'ticklabelinterpreter','latex','fontsize',FS,'tickdir','out');



%% Circumferential stress plot
figure('Name','Circumferential stress plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(7,1);

% Shell membrane theory solution
plot(SigmaTh_m(0:L/1000:L), 0:L/1000:L,...
    'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',LW); grid on; hold on;
lgnd{1} = 'Membrane theory midsurface stress';

% Shell bending theory solution
plot(SigmaTh_bi(0:L/1000:L), 0:L/1000:L,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(SigmaTh_bm(0:L/1000:L), 0:L/1000:L,...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(SigmaTh_bo(0:L/1000:L), 0:L/1000:L,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{2} = 'Bending theory stress (inner surface)';
lgnd{3} = 'Bending theory stress (midsurface)';
lgnd{4} = 'Bending theory stress (outer surface)';

% AQUINAS solution
plot(LA_out.Stresses.Sig_theta.i{1}, LA_out.Shell_Geom.z{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stresses.Sig_theta.m{1}, LA_out.Shell_Geom.z{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(LA_out.Stresses.Sig_theta.o{1}, LA_out.Shell_Geom.z{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stresses.Sig_theta.i{1}, LA_out.Shell_Geom.z{1},...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Stresses.Sig_theta.m{1}, LA_out.Shell_Geom.z{1},...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(LA_out.Stresses.Sig_theta.o{1}, LA_out.Shell_Geom.z{1},...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{5} = 'AQUINAS FE solution (inner surface)';
lgnd{6} = 'AQUINAS FE solution (midsurface)';
lgnd{7} = 'AQUINAS FE solution (outer surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Circumferential stresses in the wall [$N/mm^2$]','interpreter','latex','fontsize',FS);
ylabel('Cylinder vertical coordinate [$mm$]','interpreter','latex','fontsize',FS);
set(gca,'ticklabelinterpreter','latex','fontsize',FS,'tickdir','out');


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