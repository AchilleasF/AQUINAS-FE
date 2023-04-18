%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example LA 07: LA of a thin elliptical shell cap under uniform pressure,
% with an angle span of pi/2. The shell is closed at its top apex and has a
% meridional displacement restraint applied at its bottom edge.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% [1] Flügge W., Stresses in Shells, 2nd edition, Springer-Verlag, Berlin Heidelberg, Germany, 1973
% [2] Abaqus 6.13, Dassault Systèmes Simulia Corp, 2014.
% [3] Adam J. Sadowski, Ludovica Pototschnig, Petrina Constantinou - 'The ‘panel analysis’
%     technique in the computational study of axisymmetric thin-walled shell systems',
%     Finite Elements in Analysis and Design, Volume 152, 2018, Pages 55-68
%     https://doi.org/10.1016/j.finel.2018.07.004.
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

% Discretization (only for AQUINAS)
numelemAQ = 200;

% Plot controls
FS = 20; % [-] - Font Size (pts)
FSL = 18; % [-] - Font Size for Legends (pts)
MS = 3; % [-] - Markser size (pts)
LW = 3; % [-] - Line Width (pts)

%%%%%%%%%%%%%%%%%
% Locked Inputs %
%%%%%%%%%%%%%%%%%
% The following input variables SHOULD NOT BE ALTERED, any change to them will affect the problem definition for AQUINAS (and therefore its solution).
% The ABAQUS solution of the elliptical shell cap that will be plotted in this example is numerically hardcoded into this script (by reading a very specific external .csv file).
% Hence, any change to the following inputs will lead to divergence between the two solutions, beyond the purposes of the current comparison.
% The following input variables are only presented here, and then used to define AQUINAS's objects, to give a better understanding of the problem definition.

% Geometry
a = 1000.0; % [mm] - Radius of curvature along the primary axis of the elliptical shell
b = 500.0; % [mm] - Radius of curvature along the secondary axis of the elliptical shell
t = 1.0; % [mm] - Elliptical shell segment thickness

% Material
MatName = 'Steel'; % Material name
E = 2.0e5; % [N/mm2] - Young's modulus of steel
nu = 0.3; % [-] - Poisson ratio of steel

% Loading
p = 1.0; % [N/mm2] - Distributed normal pressure (positive when acting towards the outer surface of the shell)

%%%%%%%%%%%%%%%%%%%%%%%
% Analytical solution %
%%%%%%%%%%%%%%%%%%%%%%%
% The following bending theory solution can be found in chapter 2.2.2.2 of [1]

% Please note a difference in notation below. In AQUINAS, 'u' and 'w' are the
% global displacements parallel to the r and z axes respectively. According to
% Flugge's notation in [1], however, 'w' is the displacement normal to the shell's
% wall, while 'v' is the meridional displacement. Consequently, in the below the
% meaning of 'w' and 'v' is the OPPOSITE of that used in AQUINAS for 'w' and 'u'.

numstan = 200 + 1; % Number of Analytical Stations, necessary for desired accuracy of displacement field

% Membrane theory solution
Nphi_m = @(phi) ((p*a^2)/2)*(1./sqrt((a^2)*(sin(phi).^2)+(b^2)*(cos(phi).^2))); % Meridional membrane stress resultant
Ntheta_m = @(phi) ((p*a^2)/(2*b^2))*(((b^2)-(a^2-b^2)*(sin(phi).^2))./sqrt((a^2)*(sin(phi).^2)+(b^2)*(cos(phi).^2))); % Circumferential membrane stress resultant
SigmaPhi_m = @(phi) Nphi_m(phi)/t; % Meridional membrane stress
SigmaTheta_m = @(phi) Ntheta_m(phi)/t; % Circumferential membrane stress
EpsPhi_m = @(phi) (SigmaPhi_m(phi)-nu*SigmaTheta_m(phi))/E; % Meridional midsurface strain for membrane stress state only
EpsTheta_m = @(phi) (SigmaTheta_m(phi)-nu*SigmaPhi_m(phi))/E; % Circumferential midsurface strain for membrane stress state only


%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%

M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)));

S1 = AQUINAS_Segment_Object('type','Ellipse','material',MatName,'thickness',t,'rzbot',[a 0.0],'rztop',[0.0 b],'geom',[0 0 a b],'els',numelemAQ);

C1 = AQUINAS_Constraint_Object('rcoord',a,'zcoord',0.0,'dofs',{'w'});

P1 = AQUINAS_Distributed_Pressure_Object('segments',{S1},'type','pn','functionHandle',@() p);

SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',6);

A1 = AQUINAS_Analysis_Object('type','LA');

O1 = AQUINAS_Output_Mode_Object('mode','VarArgOut');

LA_out = AQUINAS_Protocol(M1,S1,C1,P1,SOL1,A1,O1);


%%%%%%%%%%%%%%%%%%%
% ABAQUS solution %
%%%%%%%%%%%%%%%%%%%
% Abaqus CAE [2] was chosen as a reliable FEA tool in order to obtain a solution for the problem of an elliptical cap under uniform pressure.

numelemAB = 300; % The number of elements used for the discretization, chosen so that a fine enough mesh is assigned to the boundary layer of the elliptical cap
%                (following the 'rule of thumb' that is usually adopted in such scenarios, of at least 10 elements in the boundary layer of the shell, see [3]).
%                The extent of the boundary layer was found through an approximation, using the corresponding formula for a sphere, since the radii of meridional
%                and circumferential curvature of the elliptical shell depend on the meridional coordinate.
%                The necessary mesh resolution that was evaluated using this process was then adopted for the entirety of the shell's meridian. In addition instead of the 156 elements
%                that the above 'rule of thumb' would demand to have enough elements in the boundary layer, almost double that number (300 elements) was used, since the scope of these examples is not to test ABAQUS's accuracy.

% ABAQUS data, copied from .csv file into the DataAB array
DataAB = readtable('input_Example_LA07.csv'); DataAB = table2array(DataAB);


Rcoord_AB = DataAB(:,1);
Zcoord_AB = DataAB(:,3);
Scoord_AB = DataAB(:,4);
ScoordNorm_AB = DataAB(:,5);
uAB = DataAB(:,6);
wAB = DataAB(:,7);
EpsilonPhi_ABi = DataAB(:,9);
EpsilonPhi_ABm = (DataAB(:,9)+DataAB(:,8))/2;
EpsilonPhi_ABo = DataAB(:,8);
EpsilonTheta_ABi = DataAB(:,11);
EpsilonTheta_ABm = (DataAB(:,11)+DataAB(:,10))/2;
EpsilonTheta_ABo = DataAB(:,10);
SigmaPhi_ABi = DataAB(:,13);
SigmaPhi_ABm = (DataAB(:,13)+DataAB(:,12))/2;
SigmaPhi_ABo = DataAB(:,12);
SigmaTheta_ABi = DataAB(:,15);
SigmaTheta_ABm = (DataAB(:,15)+DataAB(:,14))/2;
SigmaTheta_ABo = DataAB(:,14);
NphiAB = t*SigmaPhi_ABm;
NthetaAB = t*SigmaTheta_ABm;
MphiAB = (t^2)*(SigmaPhi_ABi-SigmaPhi_ABo)/12;
MthetaAB = (t^2)*(SigmaTheta_ABi-SigmaTheta_ABo)/12;

%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%

%% Radial displacement plot
figure('Name','Radial displacement plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(2,1);

% ABAQUS solution
plot(uAB, ScoordNorm_AB,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
lgnd{1} = 'ABAQUS FE solution';

% AQUINAS solution
plot(LA_out.DOFs.u{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',MS);
plot(LA_out.DOFs.u{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','w','MarkerSize',0.8*MS);
lgnd{2} = 'AQUINAS FE solution';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Radial displacement [$mm$]','interpreter','latex','fontsize',FS);
ylabel('Normalised arc-length distance from the apex [$-$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Axial displacement plot
figure('Name','Axial displacement plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(2,1);

% ABAQUS solution
plot(wAB, ScoordNorm_AB,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
lgnd{1} = 'ABAQUS FE solution';

% AQUINAS solution
plot(LA_out.DOFs.w{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',MS);
plot(LA_out.DOFs.w{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','w','MarkerSize',0.8*MS);
lgnd{2} = 'AQUINAS FE solution';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Axial displacement [$mm$]','interpreter','latex','fontsize',FS);
ylabel('Normalised arc-length distance from the apex [$-$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Meridional strain plot
figure('Name','Meridional strain plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(6,1);

% Shell membrane theory solution
plot(EpsPhi_m(0:(pi/2)/(numstan-1):pi/2), (0:(pi/2)/(numstan-1):pi/2)/(pi/2),...
    'Color',[0.7 0.7 0.7],'LineStyle','--','LineWidth',2*LW); grid on; hold on;
lgnd{1} = 'Membrane theory (midsurface)';

% ABAQUS solution
plot(EpsilonPhi_ABi, ScoordNorm_AB,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
plot(EpsilonPhi_ABm, ScoordNorm_AB,...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(EpsilonPhi_ABo, ScoordNorm_AB,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{2} = 'ABAQUS FE solution (inner surface)';
lgnd{3} = 'ABAQUS FE solution (midsurface)';
lgnd{4} = 'ABAQUS FE solution (outer surface)';

% AQUINAS solution
plot(LA_out.Strains.Epsilon.Eps_phi.i{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_phi.m{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_phi.o{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_phi.i{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Strains.Epsilon.Eps_phi.m{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(LA_out.Strains.Epsilon.Eps_phi.o{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{5} = 'AQUINAS FE solution (inner surface)';
lgnd{6} = 'AQUINAS FE solution (midsurface)';
lgnd{7} = 'AQUINAS FE solution (outer surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Meridional strains in the ellipse [$-$]','interpreter','latex','fontsize',FS);
ylabel('Normalised arc-length distance from the apex [$-$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Circumferential strain plot
figure('Name','Circumferential strain plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(6,1);

% Shell membrane theory solution
plot(EpsTheta_m(0:(pi/2)/(numstan-1):pi/2), (0:(pi/2)/(numstan-1):pi/2)/(pi/2),...
    'Color',[0.7 0.7 0.7],'LineStyle','--','LineWidth',2*LW); grid on; hold on;
lgnd{1} = 'Membrane theory (midsurface)';

% ABAQUS solution
plot(EpsilonTheta_ABi, ScoordNorm_AB,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(EpsilonTheta_ABm, ScoordNorm_AB,...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(EpsilonTheta_ABo, ScoordNorm_AB,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{2} = 'ABAQUS FE solution (inner surface)';
lgnd{3} = 'ABAQUS FE solution (midsurface)';
lgnd{4} = 'ABAQUS FE solution (outer surface)';

% AQUINAS solution
plot(LA_out.Strains.Epsilon.Eps_theta.i{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_theta.m{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_theta.o{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_theta.i{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Strains.Epsilon.Eps_theta.m{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(LA_out.Strains.Epsilon.Eps_theta.o{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{5} = 'AQUINAS FE solution (inner surface)';
lgnd{6} = 'AQUINAS FE solution (midsurface)';
lgnd{7} = 'AQUINAS FE solution (outer surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Circumferential strains in the ellipse [$-$]','interpreter','latex','fontsize',FS);
ylabel('Normalised arc-length distance from the apex [$-$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Meridional stress plot
figure('Name','Meridional stress plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(6,1);

% Shell membrane theory solution
plot(SigmaPhi_m(0:(pi/2)/(numstan-1):pi/2), (0:(pi/2)/(numstan-1):pi/2)/(pi/2),...
    'Color',[0.7 0.7 0.7],'LineStyle','--','LineWidth',2*LW); grid on; hold on;
lgnd{1} = 'Membrane theory (midsurface)';

% ABAQUS solution
plot(SigmaPhi_ABi, ScoordNorm_AB,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(SigmaPhi_ABm, ScoordNorm_AB,...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(SigmaPhi_ABo, ScoordNorm_AB,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{2} = 'ABAQUS FE solution (inner surface)';
lgnd{3} = 'ABAQUS FE solution (midsurface)';
lgnd{4} = 'ABAQUS FE solution (outer surface)';

% AQUINAS solution
plot(LA_out.Stresses.Sig_phi.i{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stresses.Sig_phi.m{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(LA_out.Stresses.Sig_phi.o{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stresses.Sig_phi.i{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Stresses.Sig_phi.m{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(LA_out.Stresses.Sig_phi.o{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{5} = 'AQUINAS FE solution (inner surface)';
lgnd{6} = 'AQUINAS FE solution (midsurface)';
lgnd{7} = 'AQUINAS FE solution (outer surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Meridional stresses in the ellipse [$N/mm^2$]','interpreter','latex','fontsize',FS);
ylabel('Normalised arc-length distance from the apex [$-$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Circumferential stress plot
figure('Name','Circumferential stress plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(6,1);

% Shell membrane theory solution
plot(SigmaTheta_m(0:(pi/2)/(numstan-1):pi/2), (0:(pi/2)/(numstan-1):pi/2)/(pi/2),...
    'Color',[0.7 0.7 0.7],'LineStyle','--','LineWidth',2*LW); grid on; hold on;
lgnd{1} = 'Membrane theory (midsurface)';

% ABAQUS solution
plot(SigmaTheta_ABi, ScoordNorm_AB,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(SigmaTheta_ABm, ScoordNorm_AB,...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(SigmaTheta_ABo, ScoordNorm_AB,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{2} = 'ABAQUS FE solution (inner surface)';
lgnd{3} = 'ABAQUS FE solution (midsurface)';
lgnd{4} = 'ABAQUS FE solution (outer surface)';

plot(LA_out.Stresses.Sig_theta.i{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stresses.Sig_theta.m{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(LA_out.Stresses.Sig_theta.o{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stresses.Sig_theta.i{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Stresses.Sig_theta.m{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(LA_out.Stresses.Sig_theta.o{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{5} = 'AQUINAS FE solution (inner surface)';
lgnd{6} = 'AQUINAS FE solution (midsurface)';
lgnd{7} = 'AQUINAS FE solution (outer surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Circumferential stresses in the ellipse [$N/mm^2$]','interpreter','latex','fontsize',FS);
ylabel('Normalised arc-length distance from the apex [$-$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Membrane stress resultant plot
figure('Name','Membrane stress resultant plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(6,1);

% Shell membrane theory solution
plot(Nphi_m(0:(pi/2)/(numstan-1):pi/2), (0:(pi/2)/(numstan-1):pi/2)/(pi/2),...
    'Color',[1.0 0.7 0.7],'LineStyle','--','LineWidth',2*LW); grid on; hold on;
plot(Ntheta_m(0:(pi/2)/(numstan-1):pi/2), (0:(pi/2)/(numstan-1):pi/2)/(pi/2),...
    'Color',[0.7 0.7 1.0],'LineStyle','--','LineWidth',2*LW);
lgnd{1} = 'Membrane theory meridional MSR';
lgnd{2} = 'Membrane theory circumferential MSR';

% ABAQUS solution
plot(NphiAB, ScoordNorm_AB,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','none','MarkerSize',MS);
plot(NthetaAB, ScoordNorm_AB,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','none','MarkerSize',MS);
lgnd{3} = 'ABAQUS FE meridional MSR';
lgnd{4} = 'ABAQUS FE circumferential MSR';


% AQUINAS solution
plot(LA_out.Stress_Reslts.Membrane.Nphi{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stress_Reslts.Membrane.Ntheta{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stress_Reslts.Membrane.Nphi{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Stress_Reslts.Membrane.Ntheta{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{5} = 'AQUINAS FE solution meridional MSR';
lgnd{6} = 'AQUINAS FE solution circumferential MSR';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Membrane stress resultants (MSR) in the wall [$N/mm$]','interpreter','latex','fontsize',FS);
ylabel('Normalised arc-length distance from the apex [$-$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Bending stress resultant plot
figure('Name','Bending stress resultant plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(4,1);

% ABAQUS solution
plot(MphiAB, ScoordNorm_AB,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','none','MarkerSize',MS); grid on; hold on;
plot(MthetaAB, ScoordNorm_AB,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','none','MarkerSize',MS);
lgnd{1} = 'ABAQUS FE meridional BSR';
lgnd{2} = 'ABAQUS FE circumferential BSR';

% AQUINAS solution
plot(-LA_out.Stress_Reslts.Bending.Mphi{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(-LA_out.Stress_Reslts.Bending.Mtheta{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(-LA_out.Stress_Reslts.Bending.Mphi{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(-LA_out.Stress_Reslts.Bending.Mtheta{1}, LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{3} = 'AQUINAS FE solution meridional BSR';
lgnd{4} = 'AQUINAS FE solution circumferential BSR';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Bending stress resultants (BSR) in the wall [$Nmm/mm$]','interpreter','latex','fontsize',FS);
ylabel('Normalised arc-length distance from the apex [$-$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


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