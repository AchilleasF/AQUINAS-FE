%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example GNA 04: GNA of a thin spherical cap under a ring or an apex load.
% The edges of the spherical cap are clamped. The equilibrium path
% of the apex deflection is explored.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% [1] Teng, J. G., & Rotter, J. M. (1989) "Elastic-plastic large deflection
% analysis of axisymmetric shells", Computers & structures, 31(2), 211-233
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.

%%%%%%%%%%%%%%%%%
% Preliminaries %
%%%%%%%%%%%%%%%%%
clearvars, clc, clear, close all
addpath('..','../AQUINAS_Analysis','../AQUINAS_Objects','../Maple_Source');

%%%%%%%%%%
% Inputs %
%%%%%%%%%%
% Geometry
r = 4.758; % [in] - Spherical cap radius
t = 0.01576; % [in] - Spherical cap thickness
a = 0.9; % [in] - Spherical cap half-projection to the radial axis
r_a_ring = 0.42; % [-] - r/a ratio, where r is the distance from the apex where the ring load is applied

r_ring = a*r_a_ring; % [in] - Radial distance from the apex where the ring load is apllied
phi0 = asin(a/r); % [rad] - Spherical cap angle span
h = r - r*cos(phi0); % [in] - Spherical cap projection to the axial axis

% Material
MatName = 'Steel'; % Material name
E = 10e6; % [lb/in2] - Young's modulus of steel
nu = 0.3; % [-] - Poisson ratio of steel

% Loading
P_ring = 1/2/pi/r_ring; % [N/mm] - Ring axial load (acting vertically downwards), applied at a distance r_ring from the apex
P_apex = 1/2/pi; % [N/mm] - Apex axial load (acting vertically downwards)
Pmax = 80; % [N/mm] - Maximum ring load to be reached (for a unit initially applied load this is also the maximum LPF to be considered)

% Plot controls
FS = 24; % [-] - Font Size (pts)
FSL = 20; % [-] - Font Size for Legend (pts)
MS = 4; % [-] - Markser size (pts)
LW = 2; % [-] - Line Width (pts)

%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%
M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)));

% C1 - clamped boundary condition at base of spherical cap
C1 = AQUINAS_Constraint_Object('rcoord',a,'zcoord',0,'dofs',{'u'});
C2 = AQUINAS_Constraint_Object('rcoord',a,'zcoord',0,'dofs',{'v'});
C3 = AQUINAS_Constraint_Object('rcoord',a,'zcoord',0,'dofs',{'w'});
C4 = AQUINAS_Constraint_Object('rcoord',a,'zcoord',0,'dofs',{'b'});

SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',6,'constrLagrMultF',0.0001,'dLPF',1,'maxAttempts',4,...
    'terminationConditions','B','LPFmax',Pmax,'Jd',4,'noMaxSteps',1000,'Jmax',20);

A1 = AQUINAS_Analysis_Object('type','GNA');

O1 = AQUINAS_Output_Mode_Object('mode','VarArgOut');

%% Ring load

S1 = AQUINAS_Segment_Object('type','Ellipse','material',MatName,'thickness',t,'rzbot',[r_ring h-(r-r*cos(asin(r_ring/r)))],'rztop',[0 h],'geom',[0 -r*cos(phi0) r r],'els',50);
S2 = AQUINAS_Segment_Object('type','Ellipse','material',MatName,'thickness',t,'rzbot',[a 0.0],'rztop',[r_ring h-(r-r*cos(asin(r_ring/r)))],'geom',[0 -r*cos(phi0) r r],'els',50);

F1 = AQUINAS_Nodal_Force_Object('type','Fw','rcoord',r_ring,'zcoord',h-(r-r*cos(asin(r_ring/r))),'magnitude',-P_ring);

DOF1 = AQUINAS_DOF_Tracker_Object('name','wc','segment',S1,'activeEnd','top','typeOfDOF','w');

GNA_out = AQUINAS_Protocol(M1,S1,S2,C1,C2,C3,C4,F1,DOF1,SOL1,A1,O1);

% Post processing
wc_AQ_ring = zeros(length(GNA_out.Step)+1,1);
P_AQ_ring = zeros(length(GNA_out.Step)+1,1);
for I = 1:length(GNA_out.Step)
    wc_AQ_ring(I+1) = -GNA_out.Step{I}.DOF_Trackers('wc');
    P_AQ_ring(I+1) = GNA_out.Step{I}.LPF;
end

%% Apex load

S1 = AQUINAS_Segment_Object('type','Ellipse','material',MatName,'thickness',t,'rzbot',[a 0.0],'rztop',[0 h],'geom',[0 -r*cos(phi0) r r],'els',80);

F1 = AQUINAS_Nodal_Force_Object('type','Fw','rcoord',0,'zcoord',h,'magnitude',-P_apex);

DOF1 = AQUINAS_DOF_Tracker_Object('name','wc','segment',S1,'activeEnd','top','typeOfDOF','w');

GNA_out = AQUINAS_Protocol(M1,S1,C1,C2,C3,C4,F1,DOF1,SOL1,A1,O1);

% Post processing
wc_AQ_apex = zeros(length(GNA_out.Step)+1,1);
P_AQ_apex = zeros(length(GNA_out.Step)+1,1);
for I = 1:length(GNA_out.Step)
    wc_AQ_apex(I+1) = -GNA_out.Step{I}.DOF_Trackers('wc');
    P_AQ_apex(I+1) = GNA_out.Step{I}.LPF;
end

%% Load digitised data
% Load data to recreate Fig 12 of Teng and Rotter (1989a) [1]
Data_RT_and_WZ = readtable('input_Example_GNA04.csv'); Data_RT_and_WZ = table2array(Data_RT_and_WZ);
Data_RT_ring = Data_RT_and_WZ(:,[1 2]);
Data_WZ_ring = Data_RT_and_WZ(:,[3 4]);
Data_RT_apex = Data_RT_and_WZ(:,[5 6]);
Data_WZ_apex = Data_RT_and_WZ(:,[7 8]);


%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%
%% Apex deflection plot
figure('Name','Apex deflection plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(6,1);

% Teng and Rotter (1989a)
ph1 = plot(Data_RT_ring(:,1), Data_RT_ring(:,2),...
    'LineStyle','none','Marker','^','MarkerSize',MS,'MarkerEdgeColor','k','MarkerFaceColor','r'); grid on; hold on;
lgnd{1} = 'Teng and Rotter (1989a) - Ring Load ($r/a = 0.42$) [Digitised]';

ph2 = plot(Data_RT_apex(:,1), Data_RT_apex(:,2),...
    'LineStyle','none','Marker','^','MarkerSize',MS,'MarkerEdgeColor','k','MarkerFaceColor','b');
lgnd{2} = 'Teng and Rotter (1989a) - Apex Load ($r/a = 0$) [Digitised]';

% Wood and Zienkiewicz (1977)
plot(Data_WZ_ring(:,1), Data_WZ_ring(:,2),...
    'LineStyle','none','LineWidth',LW,'Marker','o','MarkerSize',MS,'MarkerEdgeColor','k','MarkerFaceColor','r');
ph3 = plot(Data_WZ_ring(:,1), Data_WZ_ring(:,2),...
    'LineStyle','none','LineWidth',LW,'Marker','o','MarkerSize',0.8*MS,'MarkerEdgeColor','k','MarkerFaceColor','r');
lgnd{3} = 'Wood and Zienkiewicz (1977) - Ring Load ($r/a = 0.42$) [Digitised]';

plot(Data_WZ_apex(:,1), Data_WZ_apex(:,2),...
    'LineStyle','none','LineWidth',LW,'Marker','o','MarkerSize',MS,'MarkerEdgeColor','k','MarkerFaceColor','b');
ph4 = plot(Data_WZ_apex(:,1), Data_WZ_apex(:,2),...
    'LineStyle','none','LineWidth',LW,'Marker','o','MarkerSize',0.8*MS,'MarkerEdgeColor','k','MarkerFaceColor','b');
lgnd{4} = 'Wood and Zienkiewicz (1977) - Apex Load ($r/a = 0$) [Digitised]';

% AQUINAS solution
ph5 = plot(wc_AQ_ring, P_AQ_ring,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','d','MarkerSize',MS,'MarkerEdgeColor','k','MarkerFaceColor','r');
lgnd{5} = 'AQUINAS';
ph6 = plot(wc_AQ_apex, P_AQ_apex,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','d','MarkerSize',MS,'MarkerEdgeColor','k','MarkerFaceColor','b');
lgnd{6} = 'AQUINAS';

legend([ph1,ph2,ph3,ph4,ph5,ph6],lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Apex Deflection $w_{c}$  [in]','interpreter','latex','fontsize',FS);
ylabel('Total Load $P$  [lb]','interpreter','latex','fontsize',FS);
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