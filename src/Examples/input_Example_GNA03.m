%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example GNA 03: GNA of a thin spherical cap under uniform normal
% pressure q. The equilibrium path of the apex deflection is explored.
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
r = 2540.0; % [mm] - Spherical cap radius (R = 100 [in])
t = 12.7; % [mm] - Spherical cap thickness (t = 0.5 [in])
phi0 = 7.1*pi/180; % [rad] - Spherical cap angle span

a = r*sin(phi0);

% Material
MatName = 'Steel'; % Material name
E = 2.07e5; % [N/mm2] - Young's modulus of steel
nu = 0.3; % [-] - Poisson ratio of steel

% Loading
q = E/((a/t)^4); % [N/mm] - Distributed axisymmetric axial load along top edge (acting vertically downwards)

% Plot controls
FS = 20; % [-] - Font Size (pts)
MS = 3; % [-] - Markser size (pts)
LW = 3; % [-] - Line Width (pts)

%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%
M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)));

S1 = AQUINAS_Segment_Object('type','Ellipse','material',MatName,'thickness',t,'rzbot',[a 0.0],'rztop',[0 r*(1-cos(phi0))],'geom',[0 -r*cos(phi0) r r],'els',100);

% C1 - clamped boundary condition at base of spherical cap
C1 = AQUINAS_Constraint_Object('rcoord',a,'zcoord',0,'dofs',{'u'});
C2 = AQUINAS_Constraint_Object('rcoord',a,'zcoord',0,'dofs',{'v'});
C3 = AQUINAS_Constraint_Object('rcoord',a,'zcoord',0,'dofs',{'w'});
C4 = AQUINAS_Constraint_Object('rcoord',a,'zcoord',0,'dofs',{'b'});

P1 = AQUINAS_Distributed_Pressure_Object('segments',{S1},'type','pn','functionHandle',@() -q);

SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',6,'constrLagrMultF',0.0001,'dLPF',0.05,'maxAttempts',4,...
    'terminationConditions','B','LPFmax',15,'Jd',4,'bifurcationAfterLPF',Inf,'noMaxSteps',1000,'Jmax',20);

A1 = AQUINAS_Analysis_Object('type','GNA');

DOF1 = AQUINAS_DOF_Tracker_Object('name','wc','segment',S1,'activeEnd','top','typeOfDOF','w');

O1 = AQUINAS_Output_Mode_Object('mode','VarArgOut');

GNA_out = AQUINAS_Protocol(M1,S1,C1,C2,C3,C4,P1,DOF1,SOL1,A1,O1);

% Post processing
wc_AQ = zeros(length(GNA_out.Step),1);
q_AQ = zeros(length(GNA_out.Step),1);
for I = 1:length(GNA_out.Step)
    wc_AQ(I) = -GNA_out.Step{I}.DOF_Trackers('wc');
    q_AQ(I) = GNA_out.Step{I}.LPF*q;
end
wc_AQ_normalised = wc_AQ/t;
q_AQ_normalised = ((a/t)^4)*q_AQ/E;

% Load data to recreate Fig 13 of Teng and Rotter (1989a) [1]
Data_RT_and_Bathe = readtable('input_Example_GNA03.csv'); Data_RT_and_Bathe = table2array(Data_RT_and_Bathe);
Data_RT = Data_RT_and_Bathe(:,[1 2]);
Data_Bathe = Data_RT_and_Bathe(:,[3 4]);


%% Apex deflection plot
figure('Name','Apex deflection plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(3,1);

% Teng and Rotter (1989a)
ph1 = plot(Data_RT(:,1), Data_RT(:,2),...
    'LineStyle','none','Marker','^','MarkerSize',MS,'MarkerEdgeColor','k','MarkerFaceColor','k'); grid on; hold on;
lgnd{1} = 'Teng and Rotter (1989b) [Digitised]';

% Bathe et al (1975)
ph2 = plot(Data_Bathe(:,1), Data_Bathe(:,2),...
    'LineStyle','none','LineWidth',LW,'Marker','o','MarkerSize',MS,'MarkerEdgeColor','k','MarkerFaceColor','w');
plot(Data_Bathe(:,1), Data_Bathe(:,2),...
    'LineStyle','none','LineWidth',LW,'Marker','o','MarkerSize',0.8*MS,'MarkerEdgeColor','w','MarkerFaceColor','w');
lgnd{2} = 'Bathe et al (1975) [Digitised]';

% AQUINAS solution
ph3 = plot(wc_AQ_normalised, q_AQ_normalised,...
    'Color','k','LineStyle','-','LineWidth',LW,'Marker','d','MarkerSize',MS);
lgnd{3} = 'AQUINAS';

legend([ph1,ph2,ph3],lgnd,'interpreter','latex','fontsize',FS,'location','best');
xlabel('Dimensionless Deflection $w_{c}/t$  [-]','interpreter','latex','fontsize',FS);
ylabel('Dimensionless Load $(a/t)^{4}q/E$  [-]','interpreter','latex','fontsize',FS);
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