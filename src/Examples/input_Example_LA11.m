%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example LA 11: LA of a thin circular plate with clamped edges under uniform pressure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% [1] Timoshenko, Stephen, and S. Woinowsky-Krieger - Theory of Plates and Shells,
% 2nd ed, New York, McGraw-Hill, 1959.
% The analytical solution to this problem may also be found in
% several relevant classical texts. The original solution to this problem has been
% provided by Poisson - 'Memoirs of the Academy', vol. 8, Paris 1829.
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
r = 50.0; % [mm] - Circular plate radius
t = 1.0; % [mm] - Circular plate thickness

% Material
MatName = 'Steel'; % Material name
E = 2.0e5; % [N/mm2] - Young's modulus of steel
nu = 0.3; % [-] - Poisson ratio of steel

% Loading
pn = -1.0; % [N/mm2] - Distributed normal pressure (positive when acting axially upwards)

% Plot controls
FS = 18; % [-] - Font Size (pts)
FSL = 18; % [-] - Font Size for Legends (pts)
MS = 3; % [-] - Markser size (pts)
LW = 3; % [-] - Line Width (pts)


%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%
M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)));

S1 = AQUINAS_Segment_Object('type','Plate','thickness',t,'rzbot',[r 0.0],'rztop',[0.0 0.0],'material',MatName,'els',100);

C1 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'u'});
C2 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'v'});
C3 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'w'});
C4 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'b'});

P1 = AQUINAS_Distributed_Pressure_Object('segments',{S1},'type','pn','functionHandle',@() pn,'withRespectTo','-');

SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',6);

A1 = AQUINAS_Analysis_Object('type','LA');

O1 = AQUINAS_Output_Mode_Object('mode','VarArgOut');

LA_out = AQUINAS_Protocol(M1,S1,C1,C2,C3,C4,P1,SOL1,A1,O1);


%%%%%%%%%%%%%%%%%%%%%%%
% Analytical solution %
%%%%%%%%%%%%%%%%%%%%%%%

% Membrane theory solution
% There is no membrane theory solution for this problem,
% since the membrane stress resultants of a circular plate under uniform normal pressure are zero.

% Bending theory solution
D = E*t*t*t/(12.0*(1.0 - nu*nu)); % Bending stiffness
wb = @(R) pn/(64*D)*((r^2)-(R.^2)).^2; % Complete axial displacement (deflection), according to eq. 62 of [1]
KappaR = @(R) (pn/(16*D))*((r^2)-3*(R.^2)); % Radial curvature, obtained through differentiation of eq. 61 of [1]
KappaTheta = @(R) pn*((r^2)-(R.^2))/(16*D); % Circumferential curvature, equal to the slope (eq. 61 of [1]) divided by the radial coordinate R (obtained through observation of eq. 53 of [1])
MR_b = @(R) (pn/16)*((r^2)*(1+nu)-(R.^2)*(3+nu)); % Radial bending moment stress resultant, according to eq. 63 of [1]
MTh_b = @(R) (pn/16)*((r^2)*(1+nu)-(R.^2)*(1+3*nu)); % Circumferential bending moment stress resultant, according to eq. 64 of [1]
% The following stresses-strains are obtained from general shell theory relations, assuming membrane stress resultants equal to zero
EpsR_bt = @(R) 0.0 + 0.5*t*KappaR(R); % Complete radial strain - top surface
EpsR_bm = @(R) zeros(size(R)); % Complete radial strain - midsurface
EpsR_bb = @(R) 0.0 - 0.5*t*KappaR(R); % Complete radial strain - bottom surface
EpsTh_bt = @(R) 0.0 + 0.5*t*KappaTheta(R); % Complete circumferential strain - top surface
EpsTh_bm = @(R) zeros(size(R)); % Complete circumferential strain - midsurface
EpsTh_bb = @(R) 0.0 - 0.5*t*KappaTheta(R); % Complete circumferential strain - bottom surface
SigmaR_bt = @(R) 0.0 + 6.0*MR_b(R)/(t*t); % Complete radial membrane stress - top surface
SigmaR_bm = @(R) zeros(size(R)); % Complete radial membrane stress - midfsurface
SigmaR_bb = @(R) 0.0 - 6.0*MR_b(R)/(t*t); % Complete radial membrane stress - bottom surface
SigmaTh_bt = @(R) 0.0 + 6.0*MTh_b(R)/(t*t); % Complete circumferential membrane stress - top surface
SigmaTh_bm = @(R) zeros(size(R)); % Complete circumferential membrane stress - midsurface
SigmaTh_bb = @(R) 0.0 - 6.0*MTh_b(R)/(t*t); % Complete circumferential membrane stress - bottom surface


%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%

%% Axial displacement plot
figure('Name','Axial displacement plot','NumberTitle','off','units','normalized','outerposition',[0 0 1 1],'WindowState','Maximized','color','w');
lgnd = cell(2,1);

% Shell bending theory solution
plot(0:r/1000:r, wb(0:r/1000:r),...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
lgnd{1} = 'Bending theory displacement';

% AQUINAS solution
plot(LA_out.Shell_Geom.r{1}, LA_out.DOFs.w{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',MS);
plot(LA_out.Shell_Geom.r{1}, LA_out.DOFs.w{1},...
    'LineStyle','none','Marker','o','MarkerFaceColor','w','MarkerEdgeColor','none','MarkerSize',0.8*MS);
lgnd{2} = 'AQUINAS FE solution';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Circular plate radial coordinate [$mm$]','interpreter','latex','fontsize',FS);
ylabel('Axial deflection of circular plate [$mm$]','interpreter','latex','fontsize',FS);
set(gca,'ticklabelinterpreter','latex','fontsize',FS,'tickdir','out');


%% Radial strain plot
figure('Name','Radial strain plot','NumberTitle','off','units','normalized','outerposition',[0 0 1 1],'WindowState','Maximized','color','w');
lgnd = cell(6,1);

% Shell bending theory solution
plot(0:r/1000:r, EpsR_bt(0:r/1000:r),...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerEdgeColor','none','MarkerSize',MS); grid on; hold on;
plot(0:r/1000:r, EpsR_bm(0:r/1000:r),...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(0:r/1000:r, EpsR_bb(0:r/1000:r),...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{1} = 'Bending theory strain (top surface)';
lgnd{2} = 'Bending theory strain (midsurface)';
lgnd{3} = 'Bending theory strain (bottom surface)';

% AQUINAS solution
plot(LA_out.Shell_Geom.r{1}, LA_out.Strains.Epsilon.Eps_phi.o{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Shell_Geom.r{1}, LA_out.Strains.Epsilon.Eps_phi.m{1},...
'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(LA_out.Shell_Geom.r{1}, LA_out.Strains.Epsilon.Eps_phi.i{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Shell_Geom.r{1}, LA_out.Strains.Epsilon.Eps_phi.o{1},...
    'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerEdgeColor','none','MarkerSize',0.8*MS);
plot(LA_out.Shell_Geom.r{1}, LA_out.Strains.Epsilon.Eps_phi.m{1},...
    'LineStyle','none','Marker','d','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor','none','MarkerSize',0.8*MS);
plot(LA_out.Shell_Geom.r{1}, LA_out.Strains.Epsilon.Eps_phi.i{1},...
    'LineStyle','none','Marker','s','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerSize',0.8*MS);
lgnd{4} = 'AQUINAS FE solution (top surface)';
lgnd{5} = 'AQUINAS FE solution (midsurface)';
lgnd{6} = 'AQUINAS FE solution (bottom surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Circular plate radial coordinate [$mm$]','interpreter','latex','fontsize',FS);
ylabel('Radial strains in the wall [$-$]','interpreter','latex','fontsize',FS);
set(gca,'ticklabelinterpreter','latex','fontsize',FS,'tickdir','out');


%% Circumferential strain plot
figure('Name','Circumferential strain plot','NumberTitle','off','units','normalized','outerposition',[0 0 1 1],'WindowState','Maximized','color','w');
lgnd = cell(6,1);

% Shell bending theory solution
plot(0:r/1000:r, EpsTh_bt(0:r/1000:r),...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
plot(0:r/1000:r, EpsTh_bm(0:r/1000:r),...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(0:r/1000:r, EpsTh_bb(0:r/1000:r),...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{1} = 'Bending theory strain (top surface)';
lgnd{2} = 'Bending theory strain (midsurface)';
lgnd{3} = 'Bending theory strain (bottom surface)';

% AQUINAS solution
plot(LA_out.Shell_Geom.r{1}, LA_out.Strains.Epsilon.Eps_theta.o{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Shell_Geom.r{1}, LA_out.Strains.Epsilon.Eps_theta.m{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(LA_out.Shell_Geom.r{1}, LA_out.Strains.Epsilon.Eps_theta.i{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Shell_Geom.r{1}, LA_out.Strains.Epsilon.Eps_theta.o{1},...
    'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerEdgeColor','none','MarkerSize',0.8*MS);
plot(LA_out.Shell_Geom.r{1}, LA_out.Strains.Epsilon.Eps_theta.m{1},...
    'LineStyle','none','Marker','d','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor','none','MarkerSize',0.8*MS);
plot(LA_out.Shell_Geom.r{1}, LA_out.Strains.Epsilon.Eps_theta.i{1},...
    'LineStyle','none','Marker','s','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerSize',0.8*MS);
lgnd{4} = 'AQUINAS FE solution (bottom surface)';
lgnd{5} = 'AQUINAS FE solution (midsurface)';
lgnd{6} = 'AQUINAS FE solution (top surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Circular plate radial coordinate [$mm$]','interpreter','latex','fontsize',FS);
ylabel('Circumferential strains in the plate [$-$]','interpreter','latex','fontsize',FS);
set(gca,'ticklabelinterpreter','latex','fontsize',FS,'tickdir','out');


%% Radial stress plot
figure('Name','Radial stress plot','NumberTitle','off','units','normalized','outerposition',[0 0 1 1],'WindowState','Maximized','color','w');
lgnd = cell(6,1);

% Shell bending theory solution
plot(0:r/1000:r, SigmaR_bt(0:r/1000:r),...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
plot(0:r/1000:r, SigmaR_bm(0:r/1000:r),...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(0:r/1000:r, SigmaR_bb(0:r/1000:r),...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{1} = 'Bending theory stress (top surface)';
lgnd{2} = 'Bending theory stress (midsurface)';
lgnd{3} = 'Bending theory stress (bottom surface)';

% AQUINAS solution
plot(LA_out.Shell_Geom.r{1}, LA_out.Stresses.Sig_phi.o{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Shell_Geom.r{1}, LA_out.Stresses.Sig_phi.m{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(LA_out.Shell_Geom.r{1}, LA_out.Stresses.Sig_phi.i{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Shell_Geom.r{1}, LA_out.Stresses.Sig_phi.o{1},...
    'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerEdgeColor','none','MarkerSize',0.8*MS);
plot(LA_out.Shell_Geom.r{1}, LA_out.Stresses.Sig_phi.m{1},...
    'LineStyle','none','Marker','d','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor','none','MarkerSize',0.8*MS);
plot(LA_out.Shell_Geom.r{1}, LA_out.Stresses.Sig_phi.i{1},...
    'LineStyle','none','Marker','s','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerSize',0.8*MS);
lgnd{4} = 'AQUINAS FE solution (top surface)';
lgnd{5} = 'AQUINAS FE solution (midsurface)';
lgnd{6} = 'AQUINAS FE solution (bottom surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Circular plate radial coordinate [$mm$]','interpreter','latex','fontsize',FS);
ylabel('Radial stresses in the plate [$N/mm^2$]','interpreter','latex','fontsize',FS);
set(gca,'ticklabelinterpreter','latex','fontsize',FS,'tickdir','out');


%% Circumferential stress plot
figure('Name','Circumferential stress plot','NumberTitle','off','units','normalized','outerposition',[0 0 1 1],'WindowState','Maximized','color','w');
lgnd = cell(6,1);

% Shell bending theory solution
plot(0:r/1000:r, SigmaTh_bt(0:r/1000:r),...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
plot(0:r/1000:r, SigmaTh_bm(0:r/1000:r),...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(0:r/1000:r, SigmaTh_bb(0:r/1000:r),...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{1} = 'Bending theory stress (top surface)';
lgnd{2} = 'Bending theory stress (midsurface)';
lgnd{3} = 'Bending theory stress (bottom surface)';

% AQUINAS solution
plot(LA_out.Shell_Geom.r{1}, LA_out.Stresses.Sig_theta.o{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Shell_Geom.r{1}, LA_out.Stresses.Sig_theta.m{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(LA_out.Shell_Geom.r{1}, LA_out.Stresses.Sig_theta.i{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Shell_Geom.r{1}, LA_out.Stresses.Sig_theta.o{1},...
    'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerEdgeColor','none','MarkerSize',0.8*MS);
plot(LA_out.Shell_Geom.r{1}, LA_out.Stresses.Sig_theta.m{1},...
    'LineStyle','none','Marker','d','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor','none','MarkerSize',0.8*MS);
plot(LA_out.Shell_Geom.r{1}, LA_out.Stresses.Sig_theta.i{1},...
    'LineStyle','none','Marker','s','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerSize',0.8*MS);
lgnd{4} = 'AQUINAS FE solution (top surface)';
lgnd{5} = 'AQUINAS FE solution (midsurface)';
lgnd{6} = 'AQUINAS FE solution (bottom surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Circular plate radial coordinate [$mm$]','interpreter','latex','fontsize',FS);
ylabel('Circumferential stresses in the plate [$N/mm^2$]','interpreter','latex','fontsize',FS);
set(gca,'ticklabelinterpreter','latex','fontsize',FS,'tickdir','out');


%% Bending stress resultants
figure('Name','Bending stress resultant plot','NumberTitle','off','units','normalized','outerposition',[0 0 1 1],'WindowState','Maximized','color','w');
lgnd = cell(4,1);

% Shell bending theory solution
plot(0:r/1000:r, MR_b(0:r/1000:r),...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
plot(0:r/1000:r, MTh_b(0:r/1000:r),...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{1} = 'Bending theory radial BSR';
lgnd{2} = 'Bending theory circumferential BSR';

% AQUINAS solution
plot(LA_out.Shell_Geom.r{1}, LA_out.Stress_Reslts.Bending.Mphi{1},...
'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Shell_Geom.r{1}, LA_out.Stress_Reslts.Bending.Mtheta{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Shell_Geom.r{1}, LA_out.Stress_Reslts.Bending.Mphi{1},...
    'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerEdgeColor','none','MarkerSize',0.8*MS);
plot(LA_out.Shell_Geom.r{1}, LA_out.Stress_Reslts.Bending.Mtheta{1},...
    'LineStyle','none','Marker','d','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerSize',0.8*MS);
lgnd{3} = 'AQUINAS FE solution radial BSR';
lgnd{4} = 'AQUINAS FE solution circumferential BSR';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Circular plate radial coordinate [$mm$]','interpreter','latex','fontsize',FS);
ylabel('Bending stress resultants (BSR) in the plate [$Nmm/mm$]','interpreter','latex','fontsize',FS);
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