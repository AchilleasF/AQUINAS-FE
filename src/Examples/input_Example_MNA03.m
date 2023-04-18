%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example MNA 02: MNAs of spherical caps under normal pressure.
% A BC1r boundary conditions is considered at the edge of the spherical caps.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% [1] Teng, J. G., & Rotter, J. M. (1989) "Elastic-plastic large deflection analysis of axisymmetric shells"
% Computers & structures, 31(2), 211-233
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
phi0s = [15 20 30 40 50 60]*pi/180; % [rad] - Spherical cap angle span
R = 12.5; % [mm] - Spherical cap radius
t = 1.0; % [mm] - Spherical cap thickness

% Material
MatName = 'Steel'; % Material name
E = 2e5; % [psi] - Young's modulus of steel
Et = 0; % [psi] - post yielding modulus of steel
nu = 0.3; % [-] - Poisson ratio of steel
fy = 250; % [N/mm2] - Yield stress

% Discretization - Integration
nels = [20 20 30 40 50 60]; % Number of elements along the meridian of the spherical cap
ngps = 3; % Number of Gauss integration points
nsps = 9; % Number of Simpson integration points

% Plot controls
FS = 20; % [-] - Font Size (pts)
MS = 3; % [-] - Markser size (pts)
LW = 3; % [-] - Line Width (pts)

%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%
LimitPressures = nan(length(phi0s),1);
for I = 1:length(phi0s)
    disp(['Computing spherical cap for a phi0 angle of ',num2str(phi0s(I)*180/pi),' degrees.']);
    phi0 = phi0s(I);

    M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)),'sy',fy,'ep',0,'curveDef','True');

    S1 = AQUINAS_Segment_Object('type','Ellipse','thickness',t,'rzbot',[R*sin(phi0) R*cos(phi0)],'rztop',[0.0 R],'material',MatName,'geom',[0 0 R R],'els',nels(I),'g',0.0,'meshtype','S');

    C1 = AQUINAS_Constraint_Object('rcoord',R*sin(phi0),'zcoord',R*cos(phi0),'dofs',{'u'});
    C2 = AQUINAS_Constraint_Object('rcoord',R*sin(phi0),'zcoord',R*cos(phi0),'dofs',{'w'});
    C3 = AQUINAS_Constraint_Object('rcoord',R*sin(phi0),'zcoord',R*cos(phi0),'dofs',{'beta'});

    P1 = AQUINAS_Distributed_Pressure_Object('segments',{S1},'type','pn','functionHandle',@() -1);

    SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',6,'constrLagrMultF',0.00001,'noGaussStations',3,'noSimpsonStations',9,...
        'dLPF',0.5,'dLPFmax',0.5,'maxAttempts',5,'Jd',4,'terminationConditions','AD','noMaxSteps',500,'Jmax',50);

    DOF1 = AQUINAS_DOF_Tracker_Object('name','wc','segment',S1,'activeEnd','top','typeOfDOF','w');

    A1 = AQUINAS_Analysis_Object('type','MNA');

    O1 = AQUINAS_Output_Mode_Object('mode','VarArgOut');

    MNA_out = AQUINAS_Protocol(M1,S1,C1,C2,C3,DOF1,P1,SOL1,A1,O1);
    LimitPressures(I) = MNA_out.Step{end}.LPF;
end

%%%%%%%%%
% Table %
%%%%%%%%%
%% Load data to compare with digitised limit pressures of Teng and Rotter (1989a) [1]
Data_RT = readtable('input_Example_MNA03.csv'); Data_RT = table2array(Data_RT);

%% Tabular comparison
T = table(Data_RT(:,1),Data_RT(:,2),LimitPressures,...
          'VariableNames',{'phi0 angles','Teng and Rotter (1989a) limit load q','AQUINAS limit load q'})

%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%
%% Limit loads of spherical caps of different angle spans phi0
figure('Name','Limit loads of spherical caps of different angle spans phi0','NumberTitle','off','Color','w','WindowState','maximized');
lgnd = cell(2,1);

plot(phi0s*180/pi,LimitPressures,'Color','k','LineStyle','-','LineWidth',LW,'Marker','o','MarkerSize',MS); grid on; hold on;
lgnd{1} = 'AQUINAS';
plot(Data_RT(:,1),Data_RT(:,2),'Color','r','LineStyle','none','LineWidth',LW,'Marker','d','MarkerSize',MS);
lgnd{2} = 'Teng and Rotter (1989a)';

xlabel('Spherical cap angle span $\phi_0$ [$^o$]','FontSize',FS,'Interpreter','latex');
ylabel('Limit pressure $q$ [psi]','FontSize',FS,'Interpreter','latex');
legend(lgnd,'FontSize',FS,'Interpreter','latex','Location','best');


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