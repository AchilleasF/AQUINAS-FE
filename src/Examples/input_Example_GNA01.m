%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example GNA 01: GNA of unpressurised thin cylindrical shells of varying
% length under uniform axial compression. A BC1f-S1 boundary is considered at
% the base of the cylinder and a BC2f-S3 at the top.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% [1] Sadowski A.J., Pototschnig L. & Constantinou P. (2018) "The 'panel
% analysis' technique in the computational study of axisymmetric
% thin-walled shell systems" Finite Elements in Analysis and Design, 152,
% 55-68.
% https://doi.org/10.1016/j.finel.2018.07.004
% The above contains many references to classical texts such as those of
% Koiter which are the origin of the Koiter circle theory used here.
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
Z = [1 1.5 2 3 4 5 6 7 8 9 ...
    10 15 20 30 40 50 60 70 80 90 ...
    100 200 300 400 500 600 700 800 900 ...
    1000 2000 5000 ...
    10000 20000 50000 ...
    100000 200000 500000 ...
    1000000 1500000 2000000 2500000 3000000 4000000 5000000 7000000 ...
    10000000 50000000 ...
    100000000]; % Batdorf parameters to be considered

% Material
MatName = 'Steel'; % Material name
E = 2.0e5; % [N/mm2] - Young's modulus of steel
nu = 0.3; % [-] - Poisson ratio of steel

% Loading
Ncl = E*t*t/(r*sqrt(3.0*(1.0-nu*nu))); % [N/mm] - Distributed axisymmetric axial load along top edge (acting vertically downwards)

% Enable Surrogate Optimisation for the minimisation of the critical buckling load with respect to the circumferential mode
enableSurOpt = false;

% Plot controls
FS = 24; % [-] - Font Size (pts)
MS = 6; % [-] - Markser size (pts)
LW = 3; % [-] - Line Width (pts)

%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%
nmax = 0.5*(12.0*(1.0-nu*nu))^0.25 * sqrt(r/t); % Max. circumferential mode number on the Koiter circle
Nmax = ceil(1.5 * nmax);

crit_LPFs = nan(length(Z),1); crit_circ_modes = nan(length(Z),1);
LPFs = nan(length(Z),Nmax+1); circ_modes = nan(length(Z),Nmax+1);
timerGNAs = tic;
for I = 1:length(Z)

    disp(strcat("Computing Batdorf parameter Z = ",num2str(Z(I))));
    L = sqrt(Z(I))*sqrt(r*t)/sqrt(sqrt(1-nu^2));
    if Z(I)>1e6; bifurcationCheckLPF = 0.0; else; bifurcationCheckLPF = 0.6; end

    M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)));

    S1 = AQUINAS_Segment_Object('type','Cone','material',MatName,'thickness',t,'rzbot',[r 0.0],'rztop',[r L],'els',200);

    F1 = AQUINAS_Nodal_Force_Object('type','Fw','rcoord',r,'zcoord',L,'magnitude', -Ncl);

    % BC1f-S1 boundary condition at base of cylinder
    C1 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'u'});
    C2 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'v'});
    C3 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'w'});

    % BC2f-S3 boundary condition at top of cylinder
    C5 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',L,'dofs',{'u'});
    C6 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',L,'dofs',{'v'});

    A1 = AQUINAS_Analysis_Object('type','GNA','circumferentialModes',0:Nmax);

    SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',6,'constrLagrMultF',0.0001,'dLPF',0.01,'maxAttempts',10,'terminationConditions','AC',...
    'bifurcationAfterLPF',bifurcationCheckLPF,'noMaxSteps',1000,'Jmax',30,'ksi',1e-3,'NonlinearSolver','ArcLength','simultaneousBifurcationTreatment','pickSmallestWithinRange',...
    'surrogate_optimisation',enableSurOpt,'so_Auto_Bounds',false,'so_circWave_ub',Nmax);

    O1 = AQUINAS_Output_Mode_Object('mode','VarArgOut');

    GNA_out = AQUINAS_Protocol(M1,S1,C1,C2,C3,C5,C6,F1,SOL1,A1,O1);

    % Post-processing
    if isempty(GNA_out.Step{end}.BuckledIntoMode)
        warning(strcat("'AQUINAS Warning: The cylinder has not buckled for a Batdorf parameter of ",num2str(Z(I))));
    else
        crit_LPFs(I) = GNA_out.Step{end}.LPF;
        crit_circ_modes(I) = GNA_out.Step{end}.BuckledIntoMode;
        if ~enableSurOpt
            for J = 1:(Nmax+1)
                LPFs(I,J) = crit_LPFs(I)*GNA_out.Step{end}.BifurcationMode(J-1).EigenValue{1};
                circ_modes(I,J) = GNA_out.Step{end}.BifurcationMode(J-1).Circumferential_Wave_No;
            end
        end
    end

end
timerGNAs = toc(timerGNAs);

%% Load data generated/digitised by Sadowski et al. in [1]
Data2D = readtable('input_Example_GNA01_2D.csv'); Data2D = table2array(Data2D);
Panel_ABA_Axi_Z = Data2D(:,1);
Panel_ABA_Axi_LPF = Data2D(:,2);
Panel_ABA_Z = Data2D(:,3);
Panel_ABA_LPFcr = Data2D(:,4);
Panel_ABA_ncr = Data2D(:,5);
Panel_ABA_Z_2 = Data2D(:,6);
Panel_ABA_ncr_2 = Data2D(:,7);
Panel_Yamaki_S1_Z = Data2D(:,8);
Panel_Yamaki_S1_LPF = Data2D(:,9);
Panel_Yamaki_S3_Z = Data2D(:,10);
Panel_Yamaki_S3_LPF = Data2D(:,11);
Panel_Yamaki_S1_Z_S = Data2D(:,12);
Panel_Yamaki_S1_ncr_S = Data2D(:,13);
Panel_Yamaki_S1_Z_L = Data2D(:,14);
Panel_Yamaki_S1_ncr_L = Data2D(:,15);
Panel_Yamaki_S3_Z_S = Data2D(:,16);
Panel_Yamaki_S3_ncr_S = Data2D(:,17);
Panel_Yamaki_S3_Z_L = Data2D(:,18);
Panel_Yamaki_S3_ncr_L = Data2D(:,19);
Data3D = readtable('input_Example_GNA01_3D.csv'); Data3D = table2array(Data3D);
[Panel_ABA_Z_3D,Panel_ABA_n_3D] = meshgrid(Data3D(1:end-2,1),Data3D(:,2));
Panel_ABA_LPF_3D = Data3D(:,3:end);

%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%
%% Geometrically nonlinear critical LPFs and circumferential wave-numbers plot
figure('Name','Geometrically nonlinear critical LPFs and circumferential wave-numbers','NumberTitle','off','WindowState','Maximized','color','w');
subplot(2,1,1); hold all;
ph1 = plot(Panel_ABA_Z,Panel_ABA_LPFcr,'MarkerFaceColor',[0 0 0],'Marker','o','LineWidth',2,'Color',[0 0 0]);
ph2 = plot(Panel_ABA_Axi_Z,Panel_ABA_Axi_LPF,'MarkerFaceColor',[0.4 0.4 0.4],'Marker','o','LineWidth',2,'Color',[0.4 0.4 0.4]);
ph3 = plot(Panel_Yamaki_S1_Z,Panel_Yamaki_S1_LPF,'LineStyle','--','LineWidth',3,'Color',[0.7 0.7 0.7]);
plot(Panel_Yamaki_S3_Z,Panel_Yamaki_S3_LPF,'LineStyle','--','LineWidth',3,'Color',[0.7 0.7 0.7]);
ph4 = plot(Z,crit_LPFs,'Marker','d','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','r','LineStyle','none','Color',[1 0 0]);
set(gca,'FontName','Times New Roman','FontSize',20,'XScale','log'); grid on; axis([1 1e8 0 1.8]);
set(gca,'YTick',[0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8],'YTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8'});
xlabel('Dimensionless length Z');
ylabel('Normalised buckling load N_G_N_A / N_c_l');
legend([ph1 ph2 ph3 ph4],{'FE 3D Panel Analysis','FE 2D Axisymmetric Analysis','Yamaki S1 & S3','AQUINAS'});

subplot(2,1,2); hold all;
ph1 = plot(Panel_ABA_Z,Panel_ABA_ncr,'MarkerFaceColor',[0 0 0],'Marker','o','LineWidth',2,'Color',[0 0 0]);
plot(Panel_ABA_Z_2,Panel_ABA_ncr_2,'MarkerFaceColor',[0 0 0],'Marker','o','LineWidth',2,'Color',[0 0 0],'LineStyle','--');
ph2 = plot(Panel_ABA_Axi_Z,zeros(size(Panel_ABA_Axi_Z)),'MarkerFaceColor',[0.4 0.4 0.4],'Marker','o','LineWidth',2,'Color',[0.4 0.4 0.4]);
ph3 = plot([1 1e8],0.909*sqrt(1000)*ones(1,2),'k','LineStyle','--','LineWidth',2);
ph4 = plot(Panel_Yamaki_S1_Z_S,Panel_Yamaki_S1_ncr_S,'LineStyle','--','LineWidth',3,'Color',[0.7 0.7 0.7]);
plot(Panel_Yamaki_S1_Z_L,Panel_Yamaki_S1_ncr_L,'LineStyle','--','LineWidth',3,'Color',[0.7 0.7 0.7]);
plot(Panel_Yamaki_S3_Z_S,Panel_Yamaki_S3_ncr_S,'LineStyle','--','LineWidth',3,'Color',[0.7 0.7 0.7]);
plot(Panel_Yamaki_S3_Z_L,Panel_Yamaki_S3_ncr_L,'LineStyle','--','LineWidth',3,'Color',[0.7 0.7 0.7]);
ph5 = plot(Z,crit_circ_modes,'Marker','d','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','r','LineStyle','none','Color',[1 0 0]);
xlabel('Dimensionless length Z');
ylabel('No. of circumferential full waves n');
axis([1 1e8 0 30]);
set(gca,'FontName','Times New Roman','FontSize',FS,'XScale','log'); grid on; hold on;
set(gca,'YTick',[0 5 10 15 20 25 30]);
legend([ph1 ph2 ph3 ph4 ph5],{'FE 3D Panel Analysis','FE 2D Axisymmetric Analysis','Koiter Circle Upper Bound','Yamaki S1 & S3','AQUINAS'});

%% Surfaces of geometrically nonlinear LPFs and circumferential wave-numbers
if ~enableSurOpt
    [Panel_3D_crit_LPF,Panel_3D_crit_n] = min(Panel_ABA_LPF_3D);
    [Z_3D,~] = meshgrid(Z,0:(Nmax));
    figure('Name','Surfaces of geometrically nonlinear LPFs and circumferential wave-numbers','NumberTitle','off','WindowState','Maximized','color','w');
    ax1 = subplot(2,1,1); hold on; grid on;
    ph1 = surface(Z_3D(2:end,:),circ_modes(:,2:end)',LPFs(:,2:end)','FaceAlpha',0.5); map = [linspace(0.0,0.4)',linspace(0.0,0.4)',linspace(1.0,0.4)']; colormap(ax1,map);
    ph2 = plot3(Panel_ABA_Z_3D(1,:),Panel_3D_crit_n,Panel_3D_crit_LPF,'ro--','LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor','w');
    ph3 = plot3(Z,crit_circ_modes,crit_LPFs,'kd-','LineWidth',0.8*LW,'MarkerSize',0.8*MS,'MarkerFaceColor','k');
    xlabel('Dimensionless length Z','interpreter','latex','FontName','TimesNewRoman','FontSize',FS);
%     ylabel('No. of circumferential full waves $n$','interpreter','latex','FontName','TimesNewRoman','FontSize',FS);
%     zlabel('Normalised buckling load $N_{GNA}/N_{cl}$','interpreter','latex','FontName','TimesNewRoman','FontSize',FS);
%     legend([ph2 ph3],{'Panel analysis (ABAQUS)','AQUINAS'},'interpreter','latex','FontName','TimesNewRoman','FontSize',FS);
    title('Surfaces of GNA buckling resistances $R_{GNA}/R_{cl}$ per candidate circumferential mode n','interpreter','latex','FontName','TimesNewRoman','FontSize',FS);
%     title('AQUINAS generated surface of GNA obtained buckling ratios $N_{GNA}/N_{cl}$ per circumferential mode n','interpreter','latex','FontName','TimesNewRoman','FontSize',FS);
    set(ax1,'xscale','log'); xlim([Panel_ABA_Z_3D(1,1),Panel_ABA_Z_3D(1,end)]); zlim([0 1.5]);
    hold off;
    ax2 = subplot(2,1,2); hold on; grid on;
    ph1 = surface(Panel_ABA_Z_3D,Panel_ABA_n_3D,Panel_ABA_LPF_3D,'FaceAlpha',0.5); map = [linspace(1.0,0.4)',linspace(0.0,0.4)',linspace(0.0,0.4)']; colormap(ax2,map);
    ph2 = plot3(Panel_ABA_Z_3D(1,:),Panel_3D_crit_n,Panel_3D_crit_LPF,'ro--','LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor','w');
    ph3 = plot3(Z,crit_circ_modes,crit_LPFs,'kd-','LineWidth',0.8*LW,'MarkerSize',0.8*MS,'MarkerFaceColor','k');
    xlabel('Dimensionless length Z','interpreter','latex','FontName','TimesNewRoman','FontSize',FS);
%     ylabel('No. of circumferential full waves $n$','interpreter','latex','FontName','TimesNewRoman','FontSize',FS);
%     zlabel('Normalised buckling load $N_{GNA}/N_{cl}$','interpreter','latex','FontName','TimesNewRoman','FontSize',FS);
    legend([ph2 ph3],{'Panel analysis (ABAQUS)','AQUINAS'},'interpreter','latex','FontName','TimesNewRoman','FontSize',FS);

%     title('Panel analysis (ABAQUS) generated surface of GNA obtained buckling ratios $N_{GNA}/N_{cl}$ per circumferential mode n','interpreter','latex','FontName','TimesNewRoman','FontSize',FS);
    set(ax2,'xscale','log'); xlim([Panel_ABA_Z_3D(1,1),Panel_ABA_Z_3D(1,end)]);
zlabel('test','Units','normalized','Position',[0.8 0.4],'Rotation',90)
end


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