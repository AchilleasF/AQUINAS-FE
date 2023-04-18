%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example LA 12: LA of a thin toroid shell with uniform internal pressure.
% A BC2f - S3 support is considered at the base of the toroid.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% [1] Abaqus 6.13, Dassault Systèmes Simulia Corp, 2014.
% [2] Adam J. Sadowski, Ludovica Pototschnig, Petrina Constantinou - 'The ‘panel analysis’
%     technique in the computational study of axisymmetric thin-walled shell systems',
%     Finite Elements in Analysis and Design, Volume 152, 2018, Pages 55-68
%     https://doi.org/10.1016/j.finel.2018.07.004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.s

%%%%%%%%%%%%%%%%%
% Preliminaries %
%%%%%%%%%%%%%%%%%
clearvars, clc, close all
addpath('..','../AQUINAS_Analysis','../AQUINAS_Objects','../Maple_Source');

%%%%%%%%%%
% Inputs %
%%%%%%%%%%

% Discretization (only for AQUINAS)
numelemAQ = 100;

% Plot controls
FS = 20; % [-] - Font Size (pts)
FSL = 18; % [-] - Font Size for Legends (pts)
MS = 3; % [-] - Markser size (pts)
LW = 3; % [-] - Line Width (pts)
scale = 10; % [-] - Scale of displacements for deformed shape plot

%%%%%%%%%%%%%%%%%
% Locked Inputs %
%%%%%%%%%%%%%%%%%
% The following input variables SHOULD NOT BE ALTERED, any change to them will affect the problem definition for AQUINAS (and therefore its solution).
% The ABAQUS solution of the toroid shell that will be plotted in this example is numerically hardcoded into this script (by reading a very specific external .csv file).
% Hence, any change to the following inputs will lead to divergence between the two solutions, beyond the purposes of the current comparison.
% The following input variables are only presented here, and then used to define AQUINAS's objects, to give a better understanding of the problem definition.

% Geometry
a = 500.0; % [mm] - The local radius of the circular meridian of the the toroid
h = 2*a; % [mm] - Toroid shell height (projection of meridian in z axis)
rc = 2000; % [mm] - The radial coordinate of the centre of the toroidal shell's axisymmetric cross section
t = 1.0; % [mm] - Toroid shell segment thickness

% Material
MatName = 'Steel'; % Material name
E = 2.0e5; % [N/mm2] - Young's modulus of steel
nu = 0.0; % [-] - Poisson ratio of steel

% Loading
pn = 1; % [N/mm2] - Uniformly distributed normal pressure (positive acting towards the outer surface of the shell)

%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%

M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)));

S1 = AQUINAS_Segment_Object('type','3pointArc','material',MatName,'thickness',t,'rztop',[rc h],'rz3rd',[rc-a h/2],'rzbot',[rc 0.0],'els',numelemAQ,'meshtype','E','g',0.5);
S2 = AQUINAS_Segment_Object('type','3pointArc','material',MatName,'thickness',t,'rztop',[rc h],'rz3rd',[rc+a h/2],'rzbot',[rc 0.0],'els',numelemAQ,'meshtype','E','g',0.5);

C1 = AQUINAS_Constraint_Object('rcoord',rc,'zcoord',0.0,'dofs',{'w'});

P1 = AQUINAS_Distributed_Pressure_Object('segments',{S1},'type','pn','functionHandle', @() -pn);
P2 = AQUINAS_Distributed_Pressure_Object('segments',{S2},'type','pn','functionHandle', @() pn);

SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',6);

A1 = AQUINAS_Analysis_Object('type','LA');

O1 = AQUINAS_Output_Mode_Object('mode','VarArgOut');

LA_out = AQUINAS_Protocol(M1,S1,S2,C1,P1,P2,SOL1,A1,O1);



%%%%%%%%%%%%%%%%%%%
% ABAQUS solution %
%%%%%%%%%%%%%%%%%%%
% Abaqus CAE [1] was chosen as a reliable FEA tool in order to obtain a solution for the problem of a torispherical shell under uniform internal pressure.

numelemAB = 500; % The number of elements used for the discretization, chosen so that a fine enough mesh is assigned to the boundary layer of the toroid
%                (following the 'rule of thumb' that is usually adopted in such scenarios, of at least 10 elements in the boundary layer of the shell, see [2]).
%                The extent of the boundary layer was found through an approximation, using the corresponding formula for a sphere, since the radii of meridional
%                and circumferential curvature of the toroid depend on the meridional coordinates
%                The necessary mesh resolution that was evaluated using this process was then adopted for the entirety of the shell's meridian. In addition instead of the 156 elements
%                that the above 'rule of thumb' would demand to have enough elements in the boundary layer, almost double that number (500 elements) was used, since the scope of these examples is not to test ABAQUS's accuracy.

% ABAQUS data, copied from .csv file into the DataAB array
DataAB = readtable('input_Example_LA12.csv'); DataAB = table2array(DataAB);


ScoordNorm_AB = linspace(0,1,length(DataAB(:,1))+1)';
PHIcoord_AB = 2*pi*ScoordNorm_AB - pi/2;
PHIcoord_AB(PHIcoord_AB<=0) = PHIcoord_AB(PHIcoord_AB<=0)+2*pi;
Rcoord_AB = rc + a*cos(PHIcoord_AB);
Zcoord_AB = a + a*sin(PHIcoord_AB);
Scoord_AB = ScoordNorm_AB*2*pi*a;
uAB = [DataAB(:,2); DataAB(1,2)];
wAB = [DataAB(:,3); DataAB(1,3)];
EpsilonPhi_ABo = [DataAB(:,4); DataAB(1,4)];
EpsilonPhi_ABi = [DataAB(:,5); DataAB(1,5)];
EpsilonPhi_ABm = (EpsilonPhi_ABi+EpsilonPhi_ABo)/2;
EpsilonTheta_ABo = [DataAB(:,6); DataAB(1,6)];
EpsilonTheta_ABi = [DataAB(:,7); DataAB(1,7)];
EpsilonTheta_ABm = (EpsilonTheta_ABi+EpsilonTheta_ABo)/2;
SigmaPhi_ABo = [DataAB(:,8); DataAB(1,8)];
SigmaPhi_ABi = [DataAB(:,9); DataAB(1,9)];
SigmaPhi_ABm = (SigmaPhi_ABi+SigmaPhi_ABo)/2;
SigmaTheta_ABo = [DataAB(:,10); DataAB(1,10)];
SigmaTheta_ABi = [DataAB(:,11); DataAB(1,11)];
SigmaTheta_ABm = (SigmaTheta_ABi+SigmaTheta_ABo)/2;
NphiAB = t*SigmaPhi_ABm;
NthetaAB = t*SigmaTheta_ABm;
MphiAB = (t^2)*(SigmaPhi_ABi-SigmaPhi_ABo)/12;
MthetaAB = (t^2)*(SigmaTheta_ABi-SigmaTheta_ABo)/12;


%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%
%% Deformed shape plot
figure('Name','Deformed shape plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(3,1);

% Undeformed meridional geometry
h1 = plot(LA_out.Shell_Geom.r{1}, LA_out.Shell_Geom.z{1},...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',1.5*LW,'Marker','none');  grid on; hold on; axis square;
plot(LA_out.Shell_Geom.r{2}, LA_out.Shell_Geom.z{2},...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',1.5*LW,'Marker','none');
lgnd{1} = 'Undeformed geometry';

% Shell bending theory solution
h2 = plot(Rcoord_AB+scale*uAB, Zcoord_AB+scale*wAB,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','none');
lgnd{2} = 'ABAQUS FE solution';

% AQUINAS solution
h3 = plot(LA_out.Shell_Geom.r{1}+scale*LA_out.DOFs.u{1}, LA_out.Shell_Geom.z{1}+scale*LA_out.DOFs.w{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',MS);
plot(LA_out.Shell_Geom.r{1}+scale*LA_out.DOFs.u{1}, LA_out.Shell_Geom.z{1}+scale*LA_out.DOFs.w{1},...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','w','MarkerSize',0.8*MS);
plot(LA_out.Shell_Geom.r{2}+scale*LA_out.DOFs.u{2}, LA_out.Shell_Geom.z{2}+scale*LA_out.DOFs.w{2},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',MS);
plot(LA_out.Shell_Geom.r{2}+scale*LA_out.DOFs.u{2}, LA_out.Shell_Geom.z{2}+scale*LA_out.DOFs.w{2},...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','w','MarkerSize',0.8*MS);
lgnd{3} = 'AQUINAS FE solution';

title(strcat("Deformed shape of toroid plot comparison (Displacements are scaled ",num2str(scale)," times)"));
legend([h1 h2 h3],lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Radial coordinate r [$mm$]','interpreter','latex','fontsize',FS);
ylabel('Axial coordinate z [$mm$]','interpreter','latex','fontsize',FS);
set(gca, 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Radial displacement plot
figure('Name','Radial displacement plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(2,1);

% Shell bending theory solution
plot(uAB, ScoordNorm_AB,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
lgnd{1} = 'ABAQUS FE solution';

% AQUINAS solution
% Positive meridional curvature segment (Bottom to top apex)
plot(LA_out.DOFs.u{2}, 0.5*(1-LA_out.Shell_Geom.s{2}/LA_out.Shell_Geom.s{2}(end)),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',MS);
plot(LA_out.DOFs.u{2}, 0.5*(1-LA_out.Shell_Geom.s{2}/LA_out.Shell_Geom.s{2}(end)),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','w','MarkerSize',0.8*MS);
% Negative meridional curvature segment (Top to bottom apex)
plot(LA_out.DOFs.u{1}, 0.5+0.5*LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',MS);
plot(LA_out.DOFs.u{1}, 0.5+0.5*LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','w','MarkerSize',0.8*MS);
lgnd{2} = 'AQUINAS FE solution';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Radial displacement [$mm$]','interpreter','latex','fontsize',FS);
ylabel('Normalised counter-clockwise arc-length distance from bottom apex [$-$]','interpreter','latex','fontsize',FS);
set(gca, 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Axial displacement plot
figure('Name','Axial displacement plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(2,1);

% Shell bending theory solution
plot(wAB, ScoordNorm_AB,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
lgnd{1} = 'ABAQUS FE solution';

% AQUINAS solution
% Positive meridional curvature segment (Bottom to top apex)
plot(LA_out.DOFs.w{2}, 0.5*(1-LA_out.Shell_Geom.s{2}/LA_out.Shell_Geom.s{2}(end)),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',MS);
plot(LA_out.DOFs.w{2}, 0.5*(1-LA_out.Shell_Geom.s{2}/LA_out.Shell_Geom.s{2}(end)),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','w','MarkerSize',0.8*MS);
% Negative meridional curvature segment (Top to bottom apex)
plot(LA_out.DOFs.w{1}, 0.5+0.5*LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',MS);
plot(LA_out.DOFs.w{1}, 0.5+0.5*LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','w','MarkerSize',0.8*MS);
lgnd{2} = 'AQUINAS FE solution';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Axial displacement [$mm$]','interpreter','latex','fontsize',FS);
ylabel('Normalised counter-clockwise arc-length distance from bottom apex [$-$]','interpreter','latex','fontsize',FS);
set(gca, 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Meridional strain plot
figure('Name','Meridional strain plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(6,1);

% ABAQUS solution
h1 = plot(EpsilonPhi_ABi, ScoordNorm_AB,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
h2 = plot(EpsilonPhi_ABm, ScoordNorm_AB,...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
h3 = plot(EpsilonPhi_ABo, ScoordNorm_AB,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{1} = 'ABAQUS FE solution (inner surface)';
lgnd{2} = 'ABAQUS FE solution (midsurface)';
lgnd{3} = 'ABAQUS FE solution (outer surface)';

% AQUINAS solution
% Positive meridional curvature segment (Bottom to top apex)
h4 = plot(LA_out.Strains.Epsilon.Eps_phi.i{2}, 0.5*(1-LA_out.Shell_Geom.s{2}/LA_out.Shell_Geom.s{2}(end)),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
h5 = plot(LA_out.Strains.Epsilon.Eps_phi.m{2}, 0.5*(1-LA_out.Shell_Geom.s{2}/LA_out.Shell_Geom.s{2}(end)),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
h6 = plot(LA_out.Strains.Epsilon.Eps_phi.o{2}, 0.5*(1-LA_out.Shell_Geom.s{2}/LA_out.Shell_Geom.s{2}(end)),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_phi.i{2}, 0.5*(1-LA_out.Shell_Geom.s{2}/LA_out.Shell_Geom.s{2}(end)),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Strains.Epsilon.Eps_phi.m{2}, 0.5*(1-LA_out.Shell_Geom.s{2}/LA_out.Shell_Geom.s{2}(end)),...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(LA_out.Strains.Epsilon.Eps_phi.o{2}, 0.5*(1-LA_out.Shell_Geom.s{2}/LA_out.Shell_Geom.s{2}(end)),...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
% Negative meridional curvature segment (Top to bottom apex).
% Flip inner-outer surfaces to fit ABAQUS' convention and consider as inner the surface closer to the center of primary curvature.
plot(LA_out.Strains.Epsilon.Eps_phi.o{1}, 0.5+0.5*LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_phi.m{1}, 0.5+0.5*LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_phi.i{1}, 0.5+0.5*LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_phi.o{1}, 0.5+0.5*LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Strains.Epsilon.Eps_phi.m{1}, 0.5+0.5*LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(LA_out.Strains.Epsilon.Eps_phi.i{1}, 0.5+0.5*LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{4} = 'AQUINAS FE solution (inner surface)';
lgnd{5} = 'AQUINAS FE solution (midsurface)';
lgnd{6} = 'AQUINAS FE solution (outer surface)';

legend([h1 h2 h3 h4 h5 h6],lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Meridional strains in the toroid [$-$]','interpreter','latex','fontsize',FS);
ylabel('Normalised counter-clockwise arc-length distance from bottom apex [$-$]','interpreter','latex','fontsize',FS);
set(gca, 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Circumferential strain plot
figure('Name','Circumferential strain plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(6,1);

% ABAQUS solution
h1 = plot(EpsilonTheta_ABi, ScoordNorm_AB,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
h2 = plot(EpsilonTheta_ABm, ScoordNorm_AB,...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
h3 = plot(EpsilonTheta_ABo, ScoordNorm_AB,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{1} = 'ABAQUS FE solution (inner surface)';
lgnd{2} = 'ABAQUS FE solution (midsurface)';
lgnd{3} = 'ABAQUS FE solution (outer surface)';

% AQUINAS solution
% Positive meridional curvature segment (Bottom to top apex)
h4 = plot(LA_out.Strains.Epsilon.Eps_theta.i{2}, 0.5*(1-LA_out.Shell_Geom.s{2}/LA_out.Shell_Geom.s{2}(end)),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
h5 = plot(LA_out.Strains.Epsilon.Eps_theta.m{2}, 0.5*(1-LA_out.Shell_Geom.s{2}/LA_out.Shell_Geom.s{2}(end)),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
h6 = plot(LA_out.Strains.Epsilon.Eps_theta.o{2}, 0.5*(1-LA_out.Shell_Geom.s{2}/LA_out.Shell_Geom.s{2}(end)),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_theta.i{2}, 0.5*(1-LA_out.Shell_Geom.s{2}/LA_out.Shell_Geom.s{2}(end)),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Strains.Epsilon.Eps_theta.m{2}, 0.5*(1-LA_out.Shell_Geom.s{2}/LA_out.Shell_Geom.s{2}(end)),...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(LA_out.Strains.Epsilon.Eps_theta.o{2}, 0.5*(1-LA_out.Shell_Geom.s{2}/LA_out.Shell_Geom.s{2}(end)),...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
% Negative meridional curvature segment (Top to bottom apex).
% Flip inner-outer surfaces to fit ABAQUS' convention and consider as inner the surface closer to the center of primary curvature.
plot(LA_out.Strains.Epsilon.Eps_theta.o{1}, 0.5+0.5*LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_theta.m{1}, 0.5+0.5*LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_theta.i{1}, 0.5+0.5*LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_theta.o{1}, 0.5+0.5*LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Strains.Epsilon.Eps_theta.m{1}, 0.5+0.5*LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(LA_out.Strains.Epsilon.Eps_theta.i{1}, 0.5+0.5*LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{4} = 'AQUINAS FE solution (inner surface)';
lgnd{5} = 'AQUINAS FE solution (midsurface)';
lgnd{6} = 'AQUINAS FE solution (outer surface)';

legend([h1 h2 h3 h4 h5 h6],lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Circumferential strains in the toroid [$-$]','interpreter','latex','fontsize',FS);
ylabel('Normalised counter-clockwise arc-length distance from bottom apex [$-$]','interpreter','latex','fontsize',FS);
set(gca, 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Meridional stress plot
figure('Name','Meridional stress plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(6,1);

% ABAQUS solution
h1 = plot(SigmaPhi_ABi, ScoordNorm_AB,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
h2 = plot(SigmaPhi_ABm, ScoordNorm_AB,...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
h3 = plot(SigmaPhi_ABo, ScoordNorm_AB,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{1} = 'ABAQUS FE solution (inner surface)';
lgnd{2} = 'ABAQUS FE solution (midsurface)';
lgnd{3} = 'ABAQUS FE solution (outer surface)';

% AQUINAS solution
% Positive meridional curvature segment (Bottom to top apex)
h4 = plot(LA_out.Stresses.Sig_phi.i{2}, 0.5*(1-LA_out.Shell_Geom.s{2}/LA_out.Shell_Geom.s{2}(end)),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
h5 = plot(LA_out.Stresses.Sig_phi.m{2}, 0.5*(1-LA_out.Shell_Geom.s{2}/LA_out.Shell_Geom.s{2}(end)),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
h6 = plot(LA_out.Stresses.Sig_phi.o{2}, 0.5*(1-LA_out.Shell_Geom.s{2}/LA_out.Shell_Geom.s{2}(end)),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stresses.Sig_phi.i{2}, 0.5*(1-LA_out.Shell_Geom.s{2}/LA_out.Shell_Geom.s{2}(end)),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Stresses.Sig_phi.m{2}, 0.5*(1-LA_out.Shell_Geom.s{2}/LA_out.Shell_Geom.s{2}(end)),...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(LA_out.Stresses.Sig_phi.o{2}, 0.5*(1-LA_out.Shell_Geom.s{2}/LA_out.Shell_Geom.s{2}(end)),...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
% Negative meridional curvature segment (Top to bottom apex).
% Flip inner-outer surfaces to fit ABAQUS' convention and consider as inner the surface closer to the center of primary curvature.
plot(LA_out.Stresses.Sig_phi.o{1}, 0.5+0.5*LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stresses.Sig_phi.m{1}, 0.5+0.5*LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(LA_out.Stresses.Sig_phi.i{1}, 0.5+0.5*LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stresses.Sig_phi.o{1}, 0.5+0.5*LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Stresses.Sig_phi.m{1}, 0.5+0.5*LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(LA_out.Stresses.Sig_phi.i{1}, 0.5+0.5*LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{4} = 'AQUINAS FE solution (inner surface)';
lgnd{5} = 'AQUINAS FE solution (midsurface)';
lgnd{6} = 'AQUINAS FE solution (outer surface)';

legend([h1 h2 h3 h4 h5 h6],lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Meridional stresses in the toroid [$N/mm^2$]','interpreter','latex','fontsize',FS);
ylabel('Normalised counter-clockwise arc-length distance from bottom apex [$-$]','interpreter','latex','fontsize',FS);
set(gca, 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Circumferential stress plot
figure('Name','Circumferential stress plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(6,1);

% ABAQUS solution
h1 = plot(SigmaTheta_ABi, ScoordNorm_AB,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
h2 = plot(SigmaTheta_ABm, ScoordNorm_AB,...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
h3 = plot(SigmaTheta_ABo, ScoordNorm_AB,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{1} = 'ABAQUS FE solution (inner surface)';
lgnd{2} = 'ABAQUS FE solution (midsurface)';
lgnd{3} = 'ABAQUS FE solution (outer surface)';

% AQUINAS solution
% Positive meridional curvature segment (Bottom to top apex)
h4 = plot(LA_out.Stresses.Sig_theta.i{2}, 0.5*(1-LA_out.Shell_Geom.s{2}/LA_out.Shell_Geom.s{2}(end)),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
h5 = plot(LA_out.Stresses.Sig_theta.m{2}, 0.5*(1-LA_out.Shell_Geom.s{2}/LA_out.Shell_Geom.s{2}(end)),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
h6 = plot(LA_out.Stresses.Sig_theta.o{2}, 0.5*(1-LA_out.Shell_Geom.s{2}/LA_out.Shell_Geom.s{2}(end)),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stresses.Sig_theta.i{2}, 0.5*(1-LA_out.Shell_Geom.s{2}/LA_out.Shell_Geom.s{2}(end)),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Stresses.Sig_theta.m{2}, 0.5*(1-LA_out.Shell_Geom.s{2}/LA_out.Shell_Geom.s{2}(end)),...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(LA_out.Stresses.Sig_theta.o{2}, 0.5*(1-LA_out.Shell_Geom.s{2}/LA_out.Shell_Geom.s{2}(end)),...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
% Negative meridional curvature segment (Top to bottom apex).
% Flip inner-outer surfaces to fit ABAQUS' convention and consider as inner the surface closer to the center of primary curvature.
plot(LA_out.Stresses.Sig_theta.o{1}, 0.5+0.5*LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stresses.Sig_theta.m{1}, 0.5+0.5*LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(LA_out.Stresses.Sig_theta.i{1}, 0.5+0.5*LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stresses.Sig_theta.o{1}, 0.5+0.5*LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Stresses.Sig_theta.m{1}, 0.5+0.5*LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(LA_out.Stresses.Sig_theta.i{1}, 0.5+0.5*LA_out.Shell_Geom.s{1}/LA_out.Shell_Geom.s{1}(end),...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{4} = 'AQUINAS FE solution (inner surface)';
lgnd{5} = 'AQUINAS FE solution (midsurface)';
lgnd{6} = 'AQUINAS FE solution (outer surface)';

legend([h1 h2 h3 h4 h5 h6],lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Circumferential stresses in the toroid [$N/mm^2$]','interpreter','latex','fontsize',FS);
ylabel('Normalised counter-clockwise arc-length distance from bottom apex [$-$]','interpreter','latex','fontsize',FS);
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