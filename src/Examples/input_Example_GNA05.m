%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example GNA 05: GNA of thin spherical caps under uniform normal
% pressure p0. The axisymmetric snap-through loads are explored.
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
r = 100.0; % [mm] - Spherical cap radius
t = 1.0; % [mm] - Spherical cap thickness
lambdas = [4 5 5.5 6 7 8]; % [-] - Rise parameter lambda

% Material
MatName = 'Steel'; % Material name
E = 2.0e5; % [N/mm2] - Young's modulus of steel
nu = 0.3; % [-] - Poisson ratio of steel

% Loading
p0 = 2*E*t^2/r^2/sqrt(3*(1-nu^2)); % [N/mm] - Distributed normal pressure (positive pointing inwards)

% Dependent geometric parameters
hs = (t*lambdas.^2)/4/sqrt(3*(1-nu^2));
phi0s = acos((r-hs)/r);
as = r*sin(phi0s);

% Plot controls
FS = 20; % [-] - Font Size (pts)
MS = 6; % [-] - Markser size (pts)
LW = 3; % [-] - Line Width (pts)

%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%
limit_points = nan(length(lambdas),1);
for I = 1:length(lambdas)
    disp(['Executing GNA for rise parameter ',num2str(lambdas(I))]);
    M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)));

    S1 = AQUINAS_Segment_Object('type','Ellipse','material',MatName,'thickness',t,'rzbot',[as(I) 0.0],'rztop',[0.0 hs(I)],'geom',[0 (hs(I)-r) r r],'els',100);

    % C1 - clamped boundary condition at base of spherical cap
    C1 = AQUINAS_Constraint_Object('rcoord',as(I),'zcoord',0,'dofs',{'u'});
    C2 = AQUINAS_Constraint_Object('rcoord',as(I),'zcoord',0,'dofs',{'v'});
    C3 = AQUINAS_Constraint_Object('rcoord',as(I),'zcoord',0,'dofs',{'w'});
    C4 = AQUINAS_Constraint_Object('rcoord',as(I),'zcoord',0,'dofs',{'b'});

    P1 = AQUINAS_Distributed_Pressure_Object('segments',{S1},'type','pn','functionHandle',@() -p0);

    DOF1 = AQUINAS_DOF_Tracker_Object('name','wc','segment',S1,'activeEnd','top','typeOfDOF','w');

    SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',6,'constrLagrMultF',0.0001,'dLPF',0.01,...
        'maxAttempts',4,'bifurcationAfterLPF',Inf,'terminationConditions','D','noMaxSteps',100,'Jmax',30);

    A1 = AQUINAS_Analysis_Object('type','GNA');

    O1 = AQUINAS_Output_Mode_Object('mode','VarArgOut');

    GNA_out = AQUINAS_Protocol(M1,S1,C1,C2,C3,C4,DOF1,P1,SOL1,A1,O1);

    % Post processing
    if ~strcmp(GNA_out.Step{end}.termination_cause,"Termination condition met. Limit point detected.")
        warning(strcat("AQUINAS Warning: No limit point has been found for a spherical cap with a rise parameter ",num2str(lambdas(I))));
    else
        limit_points(I) = max(GNA_out.Step{end-1}.LPF,GNA_out.Step{end}.LPF);
    end
end



%%%%%%%%%
% Table %
%%%%%%%%%
%% Tabular comparison

% Load data to recreate Table 1 of Teng and Rotter (1989a) [1]
Data_RT_Budiansky_DN = readtable('input_Example_GNA05.csv'); Data_RT_Budiansky_DN = table2array(Data_RT_Budiansky_DN);
Data_Budiansky = Data_RT_Budiansky_DN(:,2);
Data_DN = Data_RT_Budiansky_DN(:,3);
Data_RT = Data_RT_Budiansky_DN(:,4);

% The reader is cautioned that the 'p' is used for the normal pressure notation instead of the 'q'
T = table(lambdas',Data_Budiansky,Data_DN,Data_RT,limit_points,...
          'VariableNames',{'Rise parameter','Budiansky(1959) pcr/p0','Dumir&Nath(1984) pcr/p0','Teng&Rotter(1989a) pcr/p0','AQUINAS pcr/p0'})


%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%
%% Axisymmetric snap through loads for spherical caps of different rise parameters lambda under external pressure plot
figure('Name','Axisymmetric snap through loads for spherical caps of different rise parameters lambda under external pressure','NumberTitle','off','Color','w','WindowState','maximized');
lgnd = cell(4,1);

plot(lambdas,limit_points,'Color','k','LineStyle','-','LineWidth',LW,'Marker','o','MarkerSize',MS); grid on; hold on;
lgnd{1} = 'AQUINAS';
plot(Data_RT_Budiansky_DN(:,1),Data_RT,'Color','r','LineStyle','none','LineWidth',LW,'Marker','d','MarkerSize',MS);
lgnd{2} = 'Teng and Rotter (1989a)';
plot(Data_RT_Budiansky_DN(:,1),Data_Budiansky,'Color','g','LineStyle','none','LineWidth',LW,'Marker','^','MarkerSize',MS);
lgnd{3} = 'Budiansky(1959)';
plot(Data_RT_Budiansky_DN(:,1),Data_DN,'Color','b','LineStyle','none','LineWidth',LW,'Marker','s','MarkerSize',MS);
lgnd{4} = 'Dumir and Nath (1984)';

xlabel('Rise parameter $\lambda$','FontSize',FS,'Interpreter','latex');
ylabel('Normalised snap-through pressure $\frac{p_{cr}}{2E}\frac{r^2\sqrt{3(1-nu^2)}}{t^2}$ [-]','FontSize',FS,'Interpreter','latex');
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