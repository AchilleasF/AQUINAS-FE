%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example GMNA 02: GMNA of pressurised thin cylindrical shell under
% uniform axial compression. A few p/Nx pressure ratios are considered.
% A BC1f-S1 boundary is considered at the base of the cylinder and a BC3-S4 at the top.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% [1] - J.M. Rotter, Local collapse of axially compressed pressurized thin steel cylinders.
% Journal of Structural Engineering 116.7 (1990): 1955-1970.
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
Z = 142; % [-] - Dimensionless cylinder length
t = 1.0; % [mm] - Cylinder thickness
rts = [250 500 1000 2000]; % [-] - Cylinder r/t ratios to be considered (same as the ones considered in [1])
pr_fyt = [0.0 0.25 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 6.0 7.0 8.0 9.0 10.0 12.0 14.0 16.0 18.0 20.0]; % [N/mm2] p*r/fy/t ratios of internal pressure to be considered

% Material
MatName = 'Steel'; % Material name
E = 2.0e5; % [N/mm2] - Young's modulus of steel
nu = 0.3; % [-] - Poisson ratio of steel
fy = 250; % [N/mm2] - Yield stress


% Loading
N = fy*t; % [N/mm] - Distributed axisymmetric axial load along top edge (acting vertically upwards)

% Plot controls
FS = 24; % [-] - Font Size (pts)
MS = 6; % [-] - Markser size (pts)
LW = 3; % [-] - Line Width (pts)

%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%
crit_LPFs = nan(length(rts),length(pr_fyt)); crit_prts = nan(length(rts),length(pr_fyt));
for I = 1:length(rts)
    r = rts(I)*t; % Cylinder's radius for current r/t ratio
    L = sqrt(Z)*sqrt(r*t)/sqrt(sqrt(1-nu^2)); % Cylinder's length for Batdorf parameter Z
    lambda = pi*sqrt(r*t)/((3*(1-nu^2))^0.25); % Bending layer lambda

    for J = 1:length(pr_fyt)
        disp(['Computing cylinder with an r/t ratio of ',num2str(rts(I)),' for pressure ratio p*r/fy/t = ',num2str(pr_fyt(J)),' ratio']);
        p = pr_fyt(J)*fy*t/r;

        M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)),'sy',fy,'ep',0,'curveDef','True');

        S1 = AQUINAS_Segment_Object('type','Cone','material',MatName,'thickness',t,'rzbot',[r 0],'rztop',[r L],'els',ceil(6*L/lambda));

        F1 = AQUINAS_Nodal_Force_Object('type','Fw','rcoord',r,'zcoord',L,'magnitude', -N);

        P1 = AQUINAS_Distributed_Pressure_Object('segments',{S1},'type','pn','functionHandle',@() p);

        A1 = AQUINAS_Analysis_Object('type','GMNA','circumferentialModes',0);

        DOF1 = AQUINAS_DOF_Tracker_Object('name','wtop','segment',S1,'activeEnd','top','typeOfDOF','w');

        SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',8,'constrLagrMultF',0.000001,'noGaussStations',4,'noSimpsonStations',9,'dLPF',0.005,'dLPFmax',0.005,'maxAttempts',10,'Jd',4,...
            'terminationConditions','AD','noMaxSteps',1000,'Jmax',50,'consoleOutput',false,'bifurcationAfterLPF',Inf,'visualiseNonlinear',false,'DOFTrackerForCIP','wtop');

        O1 = AQUINAS_Output_Mode_Object('mode','VarArgOut');

        % Boundary conditions at base of cylinder
        C1 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'u'});
        C2 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'w'});
        % Boundary conditions at top edge of cylinder
        C3 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',L,'dofs',{'b'});

        % Submit for analysis
        GMNA_out = AQUINAS_Protocol(M1,S1,C1,C2,C3,F1,P1,SOL1,DOF1,A1,O1);
        crit_LPFs(I,J) = GMNA_out.Step{end}.LPF; crit_prts(I,J) = crit_LPFs(I,J)*p*r/t;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotter (1990) solution %
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load digitsed results from [1]
DataRotter = readtable('input_Example_GMNA02.csv'); DataRotter = table2array(DataRotter);

%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%
%% Elasto-plastic collapse predictions for different r/t ratios
figure('Name','Elasto-plastic collapse predictions for different r/t ratios','NumberTitle','off','WindowState','Maximized','color','w','PaperType','A3','PaperOrientation','landscape');
lgnd = cell(8,1);
hold on; grid on;

ph1 = plot(DataRotter(:,1),DataRotter(:,2),'LineStyle','--','Color',[1.0 0.7 0.7],'LineWidth',LW);
lgnd{1} = 'Rotter (1990) r/t = 250';
ph2 = plot(DataRotter(:,3),DataRotter(:,4),'LineStyle','--','Color',[0.7 1.0 0.7],'LineWidth',LW);
lgnd{2} = 'Rotter (1990) r/t = 500';
ph3 = plot(DataRotter(:,5),DataRotter(:,6),'LineStyle','--','Color',[0.7 0.7 1.0],'LineWidth',LW);
lgnd{3} = 'Rotter (1990) r/t = 1000';
ph4 = plot(DataRotter(:,7),DataRotter(:,8),'LineStyle','--','Color',[0.7 1.0 1.0],'LineWidth',LW);
lgnd{4} = 'Rotter (1990) r/t = 2000';
ph5 = plot(crit_prts(1,:),fy*crit_LPFs(1,:),'LineStyle','-','Color','r','LineWidth',0.5*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r');
lgnd{5} = 'AQUINAS r/t = 250';
ph6 = plot(crit_prts(2,:),fy*crit_LPFs(2,:),'LineStyle','-','Color','g','LineWidth',0.5*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','g');
lgnd{6} = 'AQUINAS r/t = 500';
ph7 = plot(crit_prts(3,:),fy*crit_LPFs(3,:),'LineStyle','-','Color','b','LineWidth',0.5*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','b');
lgnd{7} = 'AQUINAS r/t = 1000';
ph8 = plot(crit_prts(4,:),fy*crit_LPFs(4,:),'LineStyle','-','Color','c','LineWidth',0.5*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','c');
lgnd{8} = 'AQUINAS r/t = 2000';
legend([ph1 ph2 ph3 ph4 ph5 ph6 ph7 ph8],lgnd,'Interpreter','latex','FontSize',FS);
xlabel('Shell Body Circumferential Stress $\frac{pr}{t}$ [MPa]','Interpreter','latex','FontSize',FS);
ylabel('Axial Compressive Stress $\sigma_x$ [MPa]','Interpreter','latex','FontSize',FS);
set(gca,'FontSize',FS);


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