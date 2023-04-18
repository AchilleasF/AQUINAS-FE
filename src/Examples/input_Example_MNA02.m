%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example MNA 02: MNA of a circular plate under an axial load applied at its
% axis of revolution. A BC1f - S1 boundary condition is considered at the
% edge of the circular plate.
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
r = 100.0; % [mm] - Plate radius
t = 10.0; % [mm] - Plate thickness

% Material
MatName = 'Steel'; % Material name
E = 2.0e5; % [N/mm2] - Young's modulus of steel
nu = 0.3; % [-] - Poisson ratio of steel
fy = 250; % [N/mm2] - Yield stress

% Plot controls
FS = 22; % [-] - Font Size (pts)
MS = 6; % [-] - Markser size (pts)
LW = 2; % [-] - Line Width (pts)

% Discretisation and integration
nel = 10; % Number of elements used
ngps = 3; % Number of Gauss points used
nsps = 7; % Number of Simpson (through thickness) points used

% LPFs for the bending stress resultants to be plotted
plotLPFs = [1.56 2.0 3.0 4.0 5.0 6.0];
blackcales = [0.8 0.7 0.6 0.45 0.3 0.0];
markers = ['*','s','+','d','x','o'];
% Plastic moment for the circular plate
Mp = fy*t^2/25;

%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%
M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)),'sy',fy,'ep',0,'curveDef','True');

S1 = AQUINAS_Segment_Object('type','Plate','thickness',t,'rzbot',[r 0.0],'rztop',[0.0 0.0],'material',MatName,'els',nel);

C1 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'u'});
C2 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'w'});

F1 = AQUINAS_Nodal_Force_Object('type','Fw','rcoord',0,'zcoord',0,'magnitude',-Mp);

SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',6,'constrLagrMultF',0.00001,'noGaussStations',ngps,'noSimpsonStations',nsps,'noMaxSteps',500,'dLPF',0.2,'dLPFmax',0.2,'maxAttempts',0,'Jd',3,...
    'terminationConditions','AB','bifurcationAfterLPF',Inf,'Jmax',20,'LPFmax',6.5,'consoleOutput',true);

DOF1 = AQUINAS_DOF_Tracker_Object('name','wc','segment',S1,'activeEnd','top','typeOfDOF','w');

A1 = AQUINAS_Analysis_Object('type','MNA');

O1 = AQUINAS_Output_Mode_Object('mode','VarArgOut');

MNA_out = AQUINAS_Protocol(M1,S1,C1,C2,DOF1,F1,SOL1,A1,O1);

%% Post Processing
eqpath = zeros(length(MNA_out.Step)+1,2);
for I = 1:length(MNA_out.Step)
    eqpath(I+1,:) = [MNA_out.Step{I}.DOF_Trackers('wc'), MNA_out.Step{I}.LPF];
end

S = zeros(nel*ngps,1);
Mphi = zeros(nel*ngps,length(plotLPFs));
Mtheta = zeros(nel*ngps,length(plotLPFs));
for I = 1:length(plotLPFs)
    [~,istep] = min(abs(eqpath(2:end,2)-plotLPFs(I)));
    for E = 1:nel
        for J = 1:ngps
            if I == 1; S((E-1)*ngps+J) = MNA_out.Shell_Geom.seg{1}.el{E}.gp{J}.s; end
            Mphi((E-1)*ngps+J,I) = MNA_out.Step{istep}.seg{1}.el{E}.gp{J}.Stress_Reslts.Bending.Mphi;
            Mtheta((E-1)*ngps+J,I) = MNA_out.Step{istep}.seg{1}.el{E}.gp{J}.Stress_Reslts.Bending.Mtheta;
        end
    end
end


%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%
%% Plate's centroid deflection equilibrium path plot
figure('Name','Plate''s centroid deflection equilibrium path plot of cylinder under axial compression','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(1,1);

plot([-eqpath(1,1), -eqpath(end,1)/t],[2*pi, 2*pi],'Color','k','LineStyle','--','LineWidth',LW,'Marker','none'); grid on; hold on;
text(-0.5*eqpath(end,1)/t,6.6,'$P = 6.283M_p$','FontSize',FS,'Color','k','Interpreter','latex');

% AQUINAS solution
ph1 = plot(-eqpath(:,1)/t, eqpath(:,2),'Color','k','LineStyle','-','LineWidth',LW,'Marker','.','MarkerSize',MS,'MarkerEdgeColor','k','MarkerFaceColor','k');
lgnd{1} = 'AQUINAS FE solution';

xlim([-eqpath(1,1),-eqpath(end,1)/t]);
legend([ph1],lgnd,'interpreter','latex','fontsize',FS,'location','best');
xlabel('Dimensionless Deflection $w_c /t$','interpreter','latex','fontsize',FS);
ylabel('Dimensionless Load $P/M_p$','interpreter','latex','fontsize',FS);
set(gca,'ticklabelinterpreter','latex','fontsize',FS,'tickdir','out');

%% Development of meridional bending stress resultants plot
figure('Name','Development of meridional bending stress resultants plot','NumberTitle','off','WindowState','Maximized','color','w');
phs = []; lgnd = cell(length(plotLPFs),1); hold on; grid on;

for I = 1:length(plotLPFs)
    phs(I) = plot(S,Mphi(:,I),'Color',ones(3,1)*blackcales(I),'LineStyle','-','LineWidth',LW,'Marker',markers(I));
    lgnd{I} = strcat('$P/M_p = ',num2str(plotLPFs(I)),'$');
end

xlim([0.0,r]);
legend(phs,lgnd,'interpreter','latex','fontsize',FS,'location','best');
xlabel('Radial coordinate $r$ [mm]','interpreter','latex','fontsize',FS);
ylabel('Meridional Bending Moment $M_\phi$ [Nmm/mm]','interpreter','latex','fontsize',FS);
set(gca,'ticklabelinterpreter','latex','fontsize',FS,'tickdir','out');


%% Development of circumferential bending stress resultants plot
figure('Name','Development of circumferential bending stress resultants plot','NumberTitle','off','WindowState','Maximized','color','w');
phs = []; lgnd = cell(length(plotLPFs),1); hold on; grid on;

for I = 1:length(plotLPFs)
    phs(I) = plot(S,Mtheta(:,I),'Color',ones(3,1)*blackcales(I),'LineStyle','-','LineWidth',LW,'Marker',markers(I));
    lgnd{I} = strcat('$P/M_p = ',num2str(plotLPFs(I)),'$');
end

xlim([0.0,r]);
legend(phs,lgnd,'interpreter','latex','fontsize',FS,'location','best');
xlabel('Radial coordinate $r$ [mm]','interpreter','latex','fontsize',FS);
ylabel('Circumferential Bending Moment $M_\theta$ [Nmm/mm]','interpreter','latex','fontsize',FS);
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