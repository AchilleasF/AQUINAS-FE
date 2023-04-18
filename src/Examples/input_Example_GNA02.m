%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example GNA 02: GNA of the Belleville spring (conical axisymmetric shell).
% The equilibrium path of the top - inner edge deflection is explored.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% [1] Teng, J. G., & Rotter, J. M. (1989) "Elastic-plastic large deflection analysis of axisymmetric shells"
% Computers & structures, 31(2), 211-233
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
r1 = 7.62; % [cm] - Belleville spring bottom edge radius
r2 = 2.54; % [cm] - Belleville spring top edge radius
t = 0.508; % [cm] - Thickness of the conical shell
h = 1.27; % [cm] - Height of the Belleville spring (projection of the cone on the axial axis)

% Material
MatName = 'Steel'; % Material name
E = 2.1e3; % [kg/cm2] - Young's modulus of steel
nu = 0.3; % [-] - Poisson ratio of steel

% Loading
N = 1; % [kg/cm] - Top edge line load (positive acting vertically downwards)

% Plot controls
FS = 24; % [-] - Font Size (pts)
MS = 5; % [-] - Markser size (pts)
LW = 2; % [-] - Line Width (pts)

%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%
M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)));

S1 = AQUINAS_Segment_Object('type','Cone','material',MatName,'thickness',t,'rzbot',[r1 0],'rztop',[r2 h],'els',100);

% C1 - clamped boundary condition at base of Belleville spring
C1 = AQUINAS_Constraint_Object('rcoord',r1,'zcoord',0,'dofs',{'u'});
C2 = AQUINAS_Constraint_Object('rcoord',r1,'zcoord',0,'dofs',{'v'});
C3 = AQUINAS_Constraint_Object('rcoord',r1,'zcoord',0,'dofs',{'w'});
C4 = AQUINAS_Constraint_Object('rcoord',r1,'zcoord',0,'dofs',{'b'});

F1 = AQUINAS_Nodal_Force_Object('type','Fw','rcoord',r2,'zcoord',h,'magnitude',-N);

DOF1 = AQUINAS_DOF_Tracker_Object('name','wc','segment',S1,'activeEnd','top','typeOfDOF','w');

SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',6,'constrLagrMultF',0.0001,'dLPF',0.1,'maxAttempts',4,...
    'terminationConditions','B','LPFmax',1,'Jd',4,'bifurcationAfterLPF',Inf,'noMaxSteps',1000,'Jmax',20);

A1 = AQUINAS_Analysis_Object('type','GNA');

O1 = AQUINAS_Output_Mode_Object('mode','VarArgOut');

GNA_out = AQUINAS_Protocol(M1,S1,C3,F1,DOF1,SOL1,A1,O1);

%% Post processing
wc_AQ = zeros(length(GNA_out.Step)+1,1);
p_AQ = zeros(length(GNA_out.Step)+1,1);
for I = 1:length(GNA_out.Step)
    wc_AQ(I+1) = -GNA_out.Step{I}.DOF_Trackers('wc');
    p_AQ(I+1) = GNA_out.Step{I}.LPF;
end


%% Load digitised data
% Load data to recreate Fig 11 of Teng and Rotter (1989a)
Data_RT_and_Su = readtable('input_Example_GNA02.csv'); Data_RT_and_Su = table2array(Data_RT_and_Su);
Data_RT = Data_RT_and_Su(:,[1 2]);
Data_Su = Data_RT_and_Su(:,[3 4]);


%% Top - inner edge deflection plot
figure('Name','Top - inner edge deflection plot of Belleville spring','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(3,1);

% Teng and Rotter (1989a)
ph1 = plot(Data_RT(:,1), Data_RT(:,2),...
    'LineStyle','none','Marker','^','MarkerSize',MS,'MarkerEdgeColor','k','MarkerFaceColor','k'); grid on; hold on;
lgnd{1} = 'Teng and Rotter (1989a) [Digitised]';

% Surana (1982)
ph2 = plot(Data_Su(:,1), Data_Su(:,2),...
    'LineStyle','none','Marker','o','MarkerSize',MS,'MarkerEdgeColor','k','MarkerFaceColor','w');
plot(Data_Su(:,1), Data_Su(:,2),...
    'LineStyle','none','Marker','o','MarkerSize',0.8*MS,'MarkerEdgeColor','w','MarkerFaceColor','w');
lgnd{2} = 'Surana (1982) [Digitised]';

% AQUINAS solution
ph3 = plot(wc_AQ, p_AQ,...
    'Color','k','LineStyle','-','LineWidth',LW,'Marker','d','MarkerSize',MS);
lgnd{3} = 'AQUINAS';

legend([ph1,ph2,ph3],lgnd,'interpreter','latex','fontsize',FS,'location','best');
xlabel('Inner Edge Deflection $\delta$  [cm]','interpreter','latex','fontsize',FS);
ylabel('Load $N$  [kg/cm]','interpreter','latex','fontsize',FS);
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