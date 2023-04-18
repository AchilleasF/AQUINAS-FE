%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example GNA 06: GNA of thin spherical caps under uniform normal
% pressure p0. The buckling modes of the spherical caps are explored.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% [1] Teng, J. G., & Rotter, J. M. (1989) "Non-symmetric bifurcation of geometrically
% nonlinear elastic-plastic axisymmetric shells under combined loads including torsion"
% Computers & structures, 32(2), 453-475
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
r = 100.0; % [mm] - Spherical cap radius
t = 1.0; % [mm] - Spherical cap thickness
lambdas = [6 7 8 9 10 12 14 16]; % [-] - Rise parameter lambda

% Material
MatName = 'Steel'; % Material name
E = 2.0e5; % [N/mm2] - Young's modulus of steel
nu = 1/3; % [-] - Poisson ratio of steel

% Loading
p0 = 2*E*t^2/r^2/sqrt(3*(1-nu^2)); % [N/mm] - Distributed normal pressure (positive pointing inwards)

% Dependent geometric parameters
Hs = (t*lambdas.^2)/4/sqrt(3*(1-nu^2));
phi0s = acos((r-Hs)/r);
as = r*sin(phi0s);

% Plot controls
FS = 20; % [-] - Font Size (pts)
MS = 6; % [-] - Markser size (pts)
LW = 2; % [-] - Line Width (pts)

%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%
crit_LPFs = nan(length(lambdas),1); crit_circ_modes = nan(length(lambdas),1);
for I = 1:length(lambdas)
    disp(['Executing GNA for rise parameter ',num2str(lambdas(I))]);
    M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)));

    S1 = AQUINAS_Segment_Object('type','Ellipse','material',MatName,'thickness',t,'rzbot',[as(I) 0.0],'rztop',[0.0 Hs(I)],'geom',[0 (Hs(I)-r) r r],'els',10);

    % C1 - clamped boundary condition at base of spherical cap
    C1 = AQUINAS_Constraint_Object('rcoord',as(I),'zcoord',0,'dofs',{'u'});
    C2 = AQUINAS_Constraint_Object('rcoord',as(I),'zcoord',0,'dofs',{'v'});
    C3 = AQUINAS_Constraint_Object('rcoord',as(I),'zcoord',0,'dofs',{'w'});
    C4 = AQUINAS_Constraint_Object('rcoord',as(I),'zcoord',0,'dofs',{'b'});

    P1 = AQUINAS_Distributed_Pressure_Object('segments',{S1},'type','pn','functionHandle',@() -p0);

    SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',6,'constrLagrMultF',0.0001,'dLPF',0.05,...
        'maxAttempts',10,'ksi',5e-4,'terminationConditions','AC','noMaxSteps',300,'Jmax',30);

    A1 = AQUINAS_Analysis_Object('type','GNA','circumferentialModes',0:20);

    O1 = AQUINAS_Output_Mode_Object('mode','VarArgOut');

    GNA_out = AQUINAS_Protocol(M1,S1,C1,C2,C3,C4,P1,SOL1,A1,O1);

    % Post processing
    if isempty(GNA_out.Step{end}.BuckledIntoMode)
        warning(['AQUINAS Warning: The spherical cap for a rise parameter ',num2str(lambdas(I)),' has not be found to buckle']);
    else
        crit_LPFs(I) = GNA_out.Step{end}.LPF;
        crit_circ_modes(I) = GNA_out.Step{end}.BuckledIntoMode;
    end
end

%%%%%%%%%
% Table %
%%%%%%%%%
% Load data to recreate Table 3 of Teng and Rotter (1989b) [1]
Data_RT_and_Huang = readtable('input_Example_GNA06.csv'); Data_RT_and_Huang = table2array(Data_RT_and_Huang);
Data_Huang = Data_RT_and_Huang(:,[2 3]);
Data_RT = Data_RT_and_Huang(:,[4 5]);

%% Tabular comparison
T = table(lambdas',Data_Huang(:,1),Data_Huang(:,2),Data_RT(:,1),Data_RT(:,2),crit_LPFs,crit_circ_modes,...
          'VariableNames',{'Rise parameter','Huang(1964) pcr/p0','Huang(1964) ncr','Teng&Rotter(1989b) pcr/p0','Teng&Rotter(1989b) ncr','AQUINAS pcr/p0','AQUINAS ncr'})


%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%
%% Bifurcaiton resistance of spherical caps of different rise parameters lambda under external pressure plot
figure('Name','Bifurcation resistance of spherical caps of different rise parameters lambda under external pressure','NumberTitle','off','Color','w','WindowState','maximized');
lgnd = cell(3,1);

plot(lambdas,crit_LPFs,'Color','k','LineStyle','-','LineWidth',LW,'Marker','o','MarkerSize',MS); grid on; hold on;
lgnd{1} = 'AQUINAS';
plot(Data_RT_and_Huang(:,1),Data_Huang(:,1),'Color','b','LineStyle','none','LineWidth',LW,'Marker','s','MarkerSize',MS);
lgnd{2} = 'Huang (1964)';
plot(Data_RT_and_Huang(:,1),Data_RT(:,1),'Color','r','LineStyle','none','LineWidth',LW,'Marker','d','MarkerSize',MS);
lgnd{3} = 'Teng and Rotter (1989b)';

xlabel('Rise parameter $\lambda$','FontSize',FS,'Interpreter','latex');
ylabel('Normalised critical bifurcation pressure $\frac{p_{cr}}{2E}\frac{r^2\sqrt{3(1-nu^2)}}{t^2}$ [-]','FontSize',FS,'Interpreter','latex');
legend(lgnd,'FontSize',FS,'Interpreter','latex','Location','best');

%% Critical circumferential mode of spherical caps of different rise parameters lambda under external pressure plot
figure('Name','Bifurcation resistance of spherical caps of different rise parameters lambda under external pressure','NumberTitle','off','Color','w','WindowState','maximized');
lgnd = cell(3,1);

plot(lambdas,crit_circ_modes,'Color','k','LineStyle','-','LineWidth',LW,'Marker','o','MarkerSize',MS); grid on; hold on;
lgnd{1} = 'AQUINAS';
plot(Data_RT_and_Huang(:,1),Data_Huang(:,2),'Color','b','LineStyle','none','LineWidth',LW,'Marker','s','MarkerSize',MS);
lgnd{2} = 'Huang (1964)';
plot(Data_RT_and_Huang(:,1),Data_RT(:,2),'Color','r','LineStyle','none','LineWidth',LW,'Marker','d','MarkerSize',MS);
lgnd{3} = 'Teng and Rotter (1989b)';

xlabel('Rise parameter $\lambda$','FontSize',FS,'Interpreter','latex');
ylabel('Critical circumferential mode $n$','FontSize',FS,'Interpreter','latex');
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