%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example LA 10: LA of a thin multi-segment shell under normal pressure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% [1] J. Heyman, Equilibrium of Shell Structures, 1st edition, Oxford Engineering Science Series
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
numelemAQ = 100;

% Plot controls
FS = 24; % [-] - Font Size (pts)
FSL = 20; % [-] - Font Size for Legends (pts)
MS = 5; % [-] - Markser size (pts)
LW = 3; % [-] - Line Width (pts)

%%%%%%%%%%%%%%%%%
% Locked Inputs %
%%%%%%%%%%%%%%%%%
% The following input variables SHOULD NOT BE ALTERED, any change to them will affect the problem definition for AQUINAS (and therefore its solution).
% The ABAQUS solution of the elliptical shell cap that will be plotted in this example is numerically hardcoded into this script (by reading a very specific external .csv file).
% Hence, any change to the following inputs will lead to divergence between the two solutions, beyond the purposes of the current comparison.
% The following input variables are only presented here, and then used to define AQUINAS's objects, to give a better understanding of the problem definition.

% Geometry
a = 200.0; % [mm] - Radius of reference alpha, accoring to Fig. 2.4 of [1]
t = 1.0; % [mm] - Elliptical shell segment thickness
alpha = 24.3*pi/180; %[rad] alpha angle, according to Fig. 2.4 of [1]

% Material
MatName = 'Steel'; % Material name
E = 2.0e5; % [N/mm2] - Young's modulus of steel
nu = 0.3; % [-] - Poisson ratio of steel

% Loading
p = 1.0; % [N/mm2] - Distributed normal pressure (positive when acting outwards)


%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%

M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)));

S1 = AQUINAS_Segment_Object('type','Ellipse','material',MatName,'thickness',t,'rzbot',[3*a*sin(alpha) 3*a*cos(alpha)-a],'rztop',[0.0 2*a],'geom',[0 -a 3*a 3*a],'els',numelemAQ);
S2 = AQUINAS_Segment_Object('type','Ellipse','material',MatName,'thickness',t,'rzbot',[3*a*cos(alpha)-a 3*a*sin(alpha)],'rztop',[3*a*sin(alpha) 3*a*cos(alpha)-a],'geom',[2*a*sin(alpha) 2*a*sin(alpha) a a],'els',numelemAQ);
S3 = AQUINAS_Segment_Object('type','Ellipse','material',MatName,'thickness',t,'rzbot',[2*a 0.0],'rztop',[3*a*cos(alpha)-a 3*a*sin(alpha)],'geom',[-a 0 3*a 3*a],'els',numelemAQ);

C1 = AQUINAS_Constraint_Object('rcoord',2*a,'zcoord',0.0,'dofs',{'w'});

P1 = AQUINAS_Distributed_Pressure_Object('segments',{S1,S2,S3},'type','pn','functionHandle',@() p);

SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',6);

A1 = AQUINAS_Analysis_Object('type','LA');

O1 = AQUINAS_Output_Mode_Object('mode','VarArgOut');

LA_out = AQUINAS_Protocol(M1,S1,S2,S3,C1,P1,SOL1,A1,O1);

% Post-processing
Scoord_AQ = [ LA_out.Shell_Geom.s{1}; LA_out.Shell_Geom.s{1}(end) + LA_out.Shell_Geom.s{2}; (LA_out.Shell_Geom.s{1}(end)+LA_out.Shell_Geom.s{2}(end)) + LA_out.Shell_Geom.s{3}];
ScoordNorm_AQ = [ LA_out.Shell_Geom.s{1}; LA_out.Shell_Geom.s{1}(end) + LA_out.Shell_Geom.s{2}; (LA_out.Shell_Geom.s{1}(end)+LA_out.Shell_Geom.s{2}(end)) + LA_out.Shell_Geom.s{3}]/(LA_out.Shell_Geom.s{1}(end)+LA_out.Shell_Geom.s{2}(end)+LA_out.Shell_Geom.s{3}(end));
uAQ = [LA_out.DOFs.u{1}; LA_out.DOFs.u{2}; LA_out.DOFs.u{3}];
wAQ = [LA_out.DOFs.w{1}; LA_out.DOFs.w{2}; LA_out.DOFs.w{3}];
EpsilonPhi_AQo = [LA_out.Strains.Epsilon.Eps_phi.o{1}; LA_out.Strains.Epsilon.Eps_phi.o{2}; LA_out.Strains.Epsilon.Eps_phi.o{3}];
EpsilonPhi_AQm = [LA_out.Strains.Epsilon.Eps_phi.m{1}; LA_out.Strains.Epsilon.Eps_phi.m{2}; LA_out.Strains.Epsilon.Eps_phi.m{3}];
EpsilonPhi_AQi = [LA_out.Strains.Epsilon.Eps_phi.i{1}; LA_out.Strains.Epsilon.Eps_phi.i{2}; LA_out.Strains.Epsilon.Eps_phi.i{3}];
EpsilonTheta_AQo = [LA_out.Strains.Epsilon.Eps_theta.o{1}; LA_out.Strains.Epsilon.Eps_theta.o{2}; LA_out.Strains.Epsilon.Eps_theta.o{3}];
EpsilonTheta_AQm = [LA_out.Strains.Epsilon.Eps_theta.m{1}; LA_out.Strains.Epsilon.Eps_theta.m{2}; LA_out.Strains.Epsilon.Eps_theta.m{3}];
EpsilonTheta_AQi = [LA_out.Strains.Epsilon.Eps_theta.i{1}; LA_out.Strains.Epsilon.Eps_theta.i{2}; LA_out.Strains.Epsilon.Eps_theta.i{3}];
SigmaPhi_AQo = [LA_out.Stresses.Sig_phi.o{1}; LA_out.Stresses.Sig_phi.o{2}; LA_out.Stresses.Sig_phi.o{3}];
SigmaPhi_AQm = [LA_out.Stresses.Sig_phi.m{1}; LA_out.Stresses.Sig_phi.m{2}; LA_out.Stresses.Sig_phi.m{3}];
SigmaPhi_AQi = [LA_out.Stresses.Sig_phi.i{1}; LA_out.Stresses.Sig_phi.i{2}; LA_out.Stresses.Sig_phi.i{3}];
SigmaTheta_AQo = [LA_out.Stresses.Sig_theta.o{1}; LA_out.Stresses.Sig_theta.o{2}; LA_out.Stresses.Sig_theta.o{3}];
SigmaTheta_AQm = [LA_out.Stresses.Sig_theta.m{1}; LA_out.Stresses.Sig_theta.m{2}; LA_out.Stresses.Sig_theta.m{3}];
SigmaTheta_AQi = [LA_out.Stresses.Sig_theta.i{1}; LA_out.Stresses.Sig_theta.i{2}; LA_out.Stresses.Sig_theta.i{3}];

%%%%%%%%%%%%%%%%%%%
% ABAQUS solution %
%%%%%%%%%%%%%%%% [1] Heyman, Jacques%%%%
% Abaqus CAE [2] was chosen as a reliable FEA tool in order to obtain a solution for the problem of an shell wall under uniform pressure.

numelemAB = 300; % The number of elements used for the discretization, chosen so that a fine enough mesh is assigned to the boundary layer of the shell wall
%                (following the 'rule of thumb' that is usually adopted in such scenarios, of at least 10 elements in the boundary layer of the shell, see [3]).
%                The extent of the boundary layer was found through an approximation, using the corresponding formula for a sphere, since the radii of meridional
%                and circumferential curvature of the elliptical shell depend on the meridional coordinate.
%                The necessary mesh resolution that was evaluated using this process was then adopted for the entirety of the shell's meridian. In addition instead of the 156 elements
%                that the above 'rule of thumb' would demand to have enough elements in the boundary layer, almost double that number (300 elements) was used, since the scope of these examples is not to test ABAQUS's accuracy.

% ABAQUS data, copied from .csv file into the DataAB array
DataAB = readtable('input_Example_LA10.csv'); DataAB = table2array(DataAB);

Rcoord_AB = DataAB(:,1);
Scoord_AB = DataAB(:,4);
ScoordNorm_AB = DataAB(:,5);
uAB = DataAB(:,6);
wAB = DataAB(:,7);
EpsilonPhi_ABo = DataAB(:,9);
EpsilonPhi_ABm = (DataAB(:,9)+DataAB(:,8))/2;
EpsilonPhi_ABi = DataAB(:,8);
EpsilonTheta_ABo = DataAB(:,11);
EpsilonTheta_ABm = (DataAB(:,11)+DataAB(:,10))/2;
EpsilonTheta_ABi = DataAB(:,10);
SigmaPhi_ABo = DataAB(:,13);
SigmaPhi_ABm = (DataAB(:,13)+DataAB(:,12))/2;
SigmaPhi_ABi = DataAB(:,12);
SigmaTheta_ABo = DataAB(:,15);
SigmaTheta_ABm = (DataAB(:,15)+DataAB(:,14))/2;
SigmaTheta_ABi = DataAB(:,14);
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
plot(uAQ, ScoordNorm_AQ,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',MS);
plot(uAQ, ScoordNorm_AQ,...
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
plot(wAQ, ScoordNorm_AQ,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',MS);
plot(wAQ, ScoordNorm_AQ,...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','w','MarkerSize',0.8*MS);
lgnd{2} = 'AQUINAS FE solution';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Axial displacement [$mm$]','interpreter','latex','fontsize',FS);
ylabel('Normalised arc-length distance from the apex [$-$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Meridional strain plot
figure('Name','Meridional strain plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(6,1);

% ABAQUS solution
plot(EpsilonPhi_ABi, ScoordNorm_AB,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
plot(EpsilonPhi_ABm, ScoordNorm_AB,...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(EpsilonPhi_ABo, ScoordNorm_AB,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{1} = 'ABAQUS FE solution (inner surface)';
lgnd{2} = 'ABAQUS FE solution (midsurface)';
lgnd{3} = 'ABAQUS FE solution (outer surface)';

% AQUINAS solution
plot(EpsilonPhi_AQi, ScoordNorm_AQ,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(EpsilonPhi_AQm, ScoordNorm_AQ,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(EpsilonPhi_AQo, ScoordNorm_AQ,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(EpsilonPhi_AQi, ScoordNorm_AQ,...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(EpsilonPhi_AQm, ScoordNorm_AQ,...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(EpsilonPhi_AQo, ScoordNorm_AQ,...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{4} = 'AQUINAS FE solution (inner surface)';
lgnd{5} = 'AQUINAS FE solution (midsurface)';
lgnd{6} = 'AQUINAS FE solution (outer surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Meridional strains in the shell wall [$-$]','interpreter','latex','fontsize',FS);
ylabel('Normalised arc-length distance from the apex [$-$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Circumferential strain plot
figure('Name','Circumferential strain plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(6,1);

% ABAQUS solution
plot(EpsilonTheta_ABi, ScoordNorm_AB,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
plot(EpsilonTheta_ABm, ScoordNorm_AB,...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(EpsilonTheta_ABo, ScoordNorm_AB,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{1} = 'ABAQUS FE solution (inner surface)';
lgnd{2} = 'ABAQUS FE solution (midsurface)';
lgnd{3} = 'ABAQUS FE solution (outer surface)';

% AQUINAS solution
plot(EpsilonTheta_AQi, ScoordNorm_AQ,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(EpsilonTheta_AQm, ScoordNorm_AQ,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(EpsilonTheta_AQo, ScoordNorm_AQ,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(EpsilonTheta_AQi, ScoordNorm_AQ,...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(EpsilonTheta_AQm, ScoordNorm_AQ,...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(EpsilonTheta_AQo, ScoordNorm_AQ,...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{4} = 'AQUINAS FE solution (inner surface)';
lgnd{5} = 'AQUINAS FE solution (midsurface)';
lgnd{6} = 'AQUINAS FE solution (outer surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Circumferential strains in the shell wall [$-$]','interpreter','latex','fontsize',FS);
ylabel('Normalised arc-length distance from the apex [$-$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Meridional stress plot
figure('Name','Meridional stress plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(6,1);

% ABAQUS solution
plot(SigmaPhi_ABi, ScoordNorm_AB,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
plot(SigmaPhi_ABm, ScoordNorm_AB,...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(SigmaPhi_ABo, ScoordNorm_AB,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{1} = 'ABAQUS FE solution (inner surface)';
lgnd{2} = 'ABAQUS FE solution (midsurface)';
lgnd{3} = 'ABAQUS FE solution (outer surface)';

% AQUINAS solution
plot(SigmaPhi_AQi, ScoordNorm_AQ,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(SigmaPhi_AQm, ScoordNorm_AQ,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(SigmaPhi_AQo, ScoordNorm_AQ,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(SigmaPhi_AQi, ScoordNorm_AQ,...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(SigmaPhi_AQm, ScoordNorm_AQ,...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(SigmaPhi_AQo, ScoordNorm_AQ,...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{4} = 'AQUINAS FE solution (inner surface)';
lgnd{5} = 'AQUINAS FE solution (midsurface)';
lgnd{6} = 'AQUINAS FE solution (outer surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Meridional stresses in the shell wall [$N/mm^2$]','interpreter','latex','fontsize',FS);
ylabel('Normalised arc-length distance from the apex [$-$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Circumferential stress plot
figure('Name','Circumferential stress plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(6,1);

% ABAQUS solution
plot(SigmaTheta_ABi, ScoordNorm_AB,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
plot(SigmaTheta_ABm, ScoordNorm_AB,...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(SigmaTheta_ABo, ScoordNorm_AB,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{1} = 'ABAQUS FE solution (inner surface)';
lgnd{2} = 'ABAQUS FE solution (midsurface)';
lgnd{3} = 'ABAQUS FE solution (outer surface)';

plot(SigmaTheta_AQi, ScoordNorm_AQ,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(SigmaTheta_AQm, ScoordNorm_AQ,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(SigmaTheta_AQo, ScoordNorm_AQ,...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(SigmaTheta_AQi, ScoordNorm_AQ,...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(SigmaTheta_AQm, ScoordNorm_AQ,...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(SigmaTheta_AQo, ScoordNorm_AQ,...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{4} = 'AQUINAS FE solution (inner surface)';
lgnd{5} = 'AQUINAS FE solution (midsurface)';
lgnd{6} = 'AQUINAS FE solution (outer surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Circumferential stresses in the shell wall [$N/mm^2$]','interpreter','latex','fontsize',FS);
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