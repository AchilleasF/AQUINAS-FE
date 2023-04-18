%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example MNA 01: MNA of unpressurised thin cylindrical shells under
% uniform axial compression. A BC1f / S1 condition is
% assumed at the base, and a BC2f / S3 condition at the top.
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
r = 1000.0; % [mm] - Cylinder radius
t = 1.0; % [mm] - Cylinder thickness
h = 1000.0; % [mm] - Cylinder height

% Material
MatName = 'Steel'; % Material name
E = 2.0e5; % [N/mm2] - Young's modulus of steel
nu = 0.3; % [-] - Poisson ratio of steel
fy = 250; % [N/mm2] - Yield stress

% Loading
Ny = fy*t; % [N/mm] - Distributed axisymmetric axial load along top edge (acting vertically downwards)

% Plot controls
FS = 20; % [-] - Font Size (pts)
MS = 6; % [-] - Markser size (pts)
LW = 2; % [-] - Line Width (pts)

%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%
%% AQUINAS FE model
M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)),'sy',fy,'ep',0,'curveDef','True');

S1 = AQUINAS_Segment_Object('type','Cone','material',MatName,'thickness',t,'rzbot',[r 0.0],'rztop',[r h],'els',100);

F1 = AQUINAS_Nodal_Force_Object('type','Fw','rcoord',r,'zcoord',h,'magnitude', -Ny);

% Boundary conditions at base of cylinder
C1 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'u'});
C2 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'w'});

% Boundary conditions at top of cylinder
C3 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',h,'dofs',{'u'});

DOF1 = AQUINAS_DOF_Tracker_Object('name','wtop','segment',S1,'activeEnd','top','typeOfDOF','w');

A1 = AQUINAS_Analysis_Object('type','MNA');

% Solver object for complete tracing of the nonlinear equilibrium path of the cylinder
SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',6,'constrLagrMultF',0.00001,'noGaussStations',3,'noSimpsonStations',7,'dLPF',0.01,'dLPFmax',0.01,'maxAttempts',10,...
    'Jd',4,'terminationConditions','A','bifurcationAfterLPF',Inf,'noMaxSteps',300,'Jmax',100,'consoleOutput',true);

O1 = AQUINAS_Output_Mode_Object('mode','VarArgOut','level','basic');

MNA_out = AQUINAS_Protocol(M1,S1,C1,C2,C3,DOF1,F1,SOL1,A1,O1);

%% Post Processing
eqpath = zeros(length(MNA_out.Step)+1,2);
for i = 1:length(MNA_out.Step)
    eqpath(i+1,:) = [MNA_out.Step{i}.DOF_Trackers('wtop'), MNA_out.Step{i}.LPF];
end


%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%
%% Top edge deflection equilibrium path plot
figure('Name','Top edge deflection equilibrium path plot of cylinder under axial compression','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(1,1);

plot([-eqpath(1,1), -eqpath(end,1)],[1.0, 1.0],'Color','k','LineStyle','--','LineWidth',LW,'Marker','none'); grid on; hold on;

% AQUINAS solution
ph1 = plot(-eqpath(:,1), eqpath(:,2),'Color','k','LineStyle','-','LineWidth',LW,'Marker','.','MarkerSize',MS,'MarkerEdgeColor','k','MarkerFaceColor','k');
lgnd{1} = 'AQUINAS FE solution';

xlim([-eqpath(1,1),-eqpath(end,1)])
legend([ph1],lgnd,'interpreter','latex','fontsize',FS,'location','best');
xlabel('Cylinder''s top edge deflection $\delta_{top}$  [mm]','interpreter','latex','fontsize',FS);
ylabel('Normalised edge line load $\frac{N}{\sigma_y t}$ [N/mm]','interpreter','latex','fontsize',FS);
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