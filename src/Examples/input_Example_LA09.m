%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example LA 09: LA of a thin hyperboloid shell under an axial top edge load.
% The shell segment is clamped at its base.
% It should be noted that for both AQUINAS and ABAQUS the hyperbolic meridional
% geometry of the shell is approximated through a circular arc defined by 3 points.
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

%%%%%%%%%%%%%%%%%
% Locked Inputs %
%%%%%%%%%%%%%%%%%
% The following input variables SHOULD NOT BE ALTERED, any change to them will affect the problem definition for AQUINAS (and therefore its solution).
% The ABAQUS solution of the hyperboloid shell segment that will be plotted in this example is numerically hardcoded into this script (by reading a very specific external .csv file).
% Hence, any change to the following inputs will lead to divergence between the two solutions, beyond the purposes of the current comparison.
% The following input variables are only presented here, and then used to define AQUINAS's objects, to give a better understanding of the problem definition.


% Geometry
h = 1000.0; % [mm] - Hyperboloid shell height (projection of meridian in z axis)
r_in = 750.0; % [mm] - The radial coordinate (distance from the axis of revolution, R=0) to the meridian (point of maximum inwards bulge)
r_edge = 1000; % [mm] - The radial coordinate of the top and bottom edges of the hyperboloid shell segment
t = 1.0; % [mm] - Hyperboloid shell segment thickness

% Material
MatName = 'Steel'; % Material name
E = 2.0e5; % [N/mm2] - Young's modulus of steel
nu = 0.3; % [-] - Poisson ratio of steel

% Loading
N = 1000.0; % [N/mm] - Distributed axisymmetric axial load along top edge (positive acting axially upwards)


%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%

M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)));

S1 = AQUINAS_Segment_Object('type','3pointArc','material',MatName,'thickness',t,'rztop',[r_edge h],'rz3rd',[r_in h/2],'rzbot',[r_edge 0.0],'els',numelemAQ,'meshtype','E','g',0.5);

C1 = AQUINAS_Constraint_Object('rcoord',r_edge,'zcoord',0.0,'dofs',{'u'});
C2 = AQUINAS_Constraint_Object('rcoord',r_edge,'zcoord',0.0,'dofs',{'w'});
C3 = AQUINAS_Constraint_Object('rcoord',r_edge,'zcoord',0.0,'dofs',{'b'});

F1 = AQUINAS_Nodal_Force_Object('type','Fw','rcoord',r_edge,'zcoord',h,'magnitude',N);

SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',6);

A1 = AQUINAS_Analysis_Object('type','LA');

O1 = AQUINAS_Output_Mode_Object('mode','VarArgOut');

LA_out = AQUINAS_Protocol(M1,S1,C1,C2,C3,F1,SOL1,A1,O1);

%%%%%%%%%%%%%%%%%%%
% ABAQUS solution %
%%%%%%%%%%%%%%%%%%%
% Abaqus CAE [1] was chosen as a reliable FEA tool in order to obtain a solution for the problem of a hyperbolic under the effect of an axial top edge line load.

numelemAB = 300; % The number of elements used for the discretization, chosen so that a fine enough mesh is assigned to the boundary layer of the hyperboloid
%                (following the 'rule of thumb' that is usually adopted in such scenarios, of at least 10 elements in the boundary layer of the shell, see [2]).
%                The extent of the boundary layer was found through an approximation, using the corresponding formula for a sphere, since the radii of meridional
%                and circumferential curvature of the hyperboloid depend on the meridional coordinate.
%                The necessary mesh resolution that was evaluated using this process was then adopted for the entirety of the shell's meridian. In addition instead of the 156 elements
%                that the above 'rule of thumb' would demand to have enough elements in the boundary layer, almost double that number (300 elements) was used, since the scope of these examples is not to test ABAQUS's accuracy.

% ABAQUS data, copied from .csv file into the DataAB array
DataAB = readtable('input_Example_LA09.csv'); DataAB = table2array(DataAB);


Rcoord_AB = DataAB(:,1);
Zcoord_AB = DataAB(:,3);
Scoord_AB = DataAB(:,4);
ScoordNorm_AB = DataAB(:,5);
uAB = DataAB(:,6);
wAB = DataAB(:,7);
EpsilonPhi_ABo = DataAB(:,8);
EpsilonPhi_ABm = (DataAB(:,8)+DataAB(:,9))/2;
EpsilonPhi_ABi = DataAB(:,9);
EpsilonTheta_ABo = DataAB(:,10);
EpsilonTheta_ABm = (DataAB(:,10)+DataAB(:,11))/2;
EpsilonTheta_ABi = DataAB(:,11);
SigmaPhi_ABo = DataAB(:,12);
SigmaPhi_ABm = (DataAB(:,12)+DataAB(:,13))/2;
SigmaPhi_ABi = DataAB(:,13);
SigmaTheta_ABo = DataAB(:,14);
SigmaTheta_ABm = (DataAB(:,14)+DataAB(:,15))/2;
SigmaTheta_ABi = DataAB(:,15);
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

% Shell bending theory solution
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
ylabel('Normalised arc-length distance from the top edge [$-$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Axial displacement plot
figure('Name','Axial displacement plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(2,1);

% Shell bending theory solution
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
ylabel('Normalised arc-length distance from the top edge [$-$]','interpreter','latex','fontsize',FS);
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
lgnd{4} = 'AQUINAS FE solution (inner surface)';
lgnd{5} = 'AQUINAS FE solution (midsurface)';
lgnd{6} = 'AQUINAS FE solution (outer surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Meridional strains in the hyperboloid [$-$]','interpreter','latex','fontsize',FS);
ylabel('Normalised arc-length distance from the top edge [$-$]','interpreter','latex','fontsize',FS);
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
lgnd{4} = 'AQUINAS FE solution (inner surface)';
lgnd{5} = 'AQUINAS FE solution (midsurface)';
lgnd{6} = 'AQUINAS FE solution (outer surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Circumferential strains in the hyperboloid [$-$]','interpreter','latex','fontsize',FS);
ylabel('Normalised arc-length distance from the top edge [$-$]','interpreter','latex','fontsize',FS);
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
lgnd{4} = 'AQUINAS FE solution (inner surface)';
lgnd{5} = 'AQUINAS FE solution (midsurface)';
lgnd{6} = 'AQUINAS FE solution (outer surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Meridional stresses in the hyperboloid [$N/mm^2$]','interpreter','latex','fontsize',FS);
ylabel('Normalised arc-length distance from the top edge [$-$]','interpreter','latex','fontsize',FS);
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
lgnd{4} = 'AQUINAS FE solution (inner surface)';
lgnd{5} = 'AQUINAS FE solution (midsurface)';
lgnd{6} = 'AQUINAS FE solution (outer surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Circumferential stresses in the hyperboloid [$N/mm^2$]','interpreter','latex','fontsize',FS);
ylabel('Normalised arc-length distance from the top edge [$-$]','interpreter','latex','fontsize',FS);
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