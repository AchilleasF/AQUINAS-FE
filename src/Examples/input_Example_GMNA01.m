%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example GMNA 01: GMNA of unpressurised thin cylindrical shells under
% uniform axial compression. A BC1r-C1 boundary is considered at the
% base of the cylinder and a BC2r-C3 at the top.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.

%%%%%%%%%%%%%%%%%
% Preliminaries %
%%%%%%%%%%%%%%%%%
clearvars, clc, close all, warning('off','all')
addpath('..','../AQUINAS_Analysis','../AQUINAS_Objects','../Maple_Source');

%%%%%%%%%%
% Inputs %
%%%%%%%%%%
% Geometry
Omega = 1; % [-] - Dimensionless cylinder length
t = 1.0; % [mm] - Cylinders thickness
rts = [10:10:200 220:20:500 550:50:1000]'; % [-] - Cylinders r/t ratios

% Material
MatName = 'Steel'; % Material name
E = 2.0e5; % [N/mm2] - Young's modulus of steel
nu = 0.3; % [-] - Poisson ratio of steel
fy = 250; % [N/mm2] - Yield stress

% Loading
N = fy*t; % [N/mm] - Distributed axisymmetric axial load along top edge (acting vertically upwards)

% Plateau tolerance
plat_tol = 1e-2; % tolerance to check for plateauing portion of equilibrium path

% Plot controls
FS = 24; % [-] - Font Size (pts)
FSL = 14; % [-] - Fond Size for Legend
MS = 8; % [-] - Markser size (pts)
LW = 3; % [-] - Line Width (pts)

%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%
Nys = nan(length(rts),1); Ncls = nan(length(rts),1); Ncrs = nan(length(rts),1); Ncrs_n0 = nan(length(rts),1); Ncrs_gn = nan(length(rts),1);
GNA_CritLPFs = nan(length(rts),1); GMNA_CritLPFs = nan(length(rts),1); GMNA_CritLPFs_n0 = nan(length(rts),1);
GNA_CritNModes = nan(length(rts),1); GMNA_CritNModes = nan(length(rts),1); GMNA_modes = nan(length(rts),ceil(1.5*0.5*(12.0*(1.0-nu*nu))^0.25 * sqrt(rts(end)))+1); GNA_modes = nan(length(rts),ceil(1.5*0.5*(12.0*(1.0-nu*nu))^0.25 * sqrt(rts(end)))+1);
termination_cause_GMNAs = cell(length(rts),1);  termination_cause_GMNAs_n0 = cell(length(rts),1);

for I = 1:length(rts)
    disp(['Computing r/t = ',num2str(rts(I))]);
    r = rts(I)*t;
    h = (Omega*r*sqrt(r/t));
    lambda = pi*sqrt(r*t)/((3*(1-nu^2))^0.25);
    Nys(I) = fy*t;
    Ncls(I) = E*t*t/(r*sqrt(3.0*(1.0-nu*nu)));
    if Nys(I) < Ncls(I); N = Nys(I); else; N = Ncls(I); end
    nmax = 0.5*(12.0*(1.0-nu*nu))^0.25 * sqrt(r/t); % Max. circumferential mode number on the Koiter circle
    Nmax = ceil(1.5 * nmax);

    M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)),'sy',fy,'ep',0,'curveDef','True');

    % Segment generation. Coarser discretisation for longer cylinders (in order to avoid running out of memory)
    if r/t < 500; nelem = ceil(8*h/lambda); else; nelem = ceil(6*h/lambda); end
    S1 = AQUINAS_Segment_Object('type','Cone','material',MatName,'thickness',t,'rzbot',[r 0],'rztop',[r h],'els',nelem);

    F1 = AQUINAS_Nodal_Force_Object('type','Fw','rcoord',r,'zcoord',h,'magnitude', -N);
    F2 = AQUINAS_Nodal_Force_Object('type','Fw','rcoord',r,'zcoord',h,'magnitude', -Ncls(I));

    A1 = AQUINAS_Analysis_Object('type','GMNA','circumferentialModes',0:Nmax);
    A2 = AQUINAS_Analysis_Object('type','GMNA','circumferentialModes',0);
    A3 = AQUINAS_Analysis_Object('type','GNA','circumferentialModes',0:Nmax);

    DOF1 = AQUINAS_DOF_Tracker_Object('name','wtop','segment',S1,'activeEnd','top','typeOfDOF','w');

    % Reduced increment step for materially nonlinear solver
    SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',6,'constrLagrMultF',0.000001,'noGaussStations',4,'noSimpsonStations',9,'dLPF',0.005,'dLPFmax',0.005,'maxAttempts',10,'Jd',4,'checkAxisymStability',true,...
        'terminationConditions','ACD','noMaxSteps',1000,'Jmax',50,'ksi',1e-3,'zetaT',2e-4,'consoleOutput',false,'simultaneousBifurcationTreatment','pickSmallestWithinRange','bifurcationAfterLPF',0.7,'visualiseNonlinear',false);
    SOL2 = AQUINAS_Solver_Object('compiler','C++','noThreads',6,'constrLagrMultF',0.000001,'noGaussStations',4,'noSimpsonStations',9,'dLPF',0.02,'dLPFmax',0.02,'maxAttempts',10,'Jd',4,'checkAxisymStability',true,...
        'terminationConditions','ACD','noMaxSteps',1000,'Jmax',50,'ksi',1e-3,'zetaT',2e-3,'consoleOutput',false,'simultaneousBifurcationTreatment','pickSmallestWithinRange','bifurcationAfterLPF',0.7,'visualiseNonlinear',false);

    O1 = AQUINAS_Output_Mode_Object('mode','VarArgOut','level','basic');

    % Boundary conditions at base of cylinder
    C1 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'u'});
    C2 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'v'});
    C3 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'w'});
    C4 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'b'});

    % Boundary conditions at top of cylinder
    C5 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',h,'dofs',{'u'});
    C6 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',h,'dofs',{'v'});
    C7 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',h,'dofs',{'b'});

    % Submit for analysis
    GMNA_out = AQUINAS_Protocol(M1,S1,C1,C2,C3,C4,C5,C6,C7,F1,SOL1,DOF1,A1,O1);
    % Results post-processing
    termination_cause_GMNAs{I} = GMNA_out.Step{end}.termination_cause;
    GMNA_CritLPFs(I) = max(GMNA_out.Step{end}.LPF,GMNA_out.Step{end-1}.LPF); % Considering the case of a limit point, the cylinder's strength should be taken as the highest LPF achieved
    if ~isempty(GMNA_out.Step{end}.BuckledIntoMode)
        GMNA_CritNModes(I) = GMNA_out.Step{end}.BuckledIntoMode;
    else
        % Check for plateauing region of equilibrium path, and accept result if encountered
        slopeEnd   = abs((GMNA_out.Step{end}.LPF - GMNA_out.Step{end-1}.LPF)/(GMNA_out.Step{end}.l-GMNA_out.Step{end-1}.l));
        slopeEndm1 = abs((GMNA_out.Step{end-1}.LPF - GMNA_out.Step{end-2}.LPF)/(GMNA_out.Step{end-1}.l-GMNA_out.Step{end-2}.l));
        slopeEndm2 = abs((GMNA_out.Step{end-2}.LPF - GMNA_out.Step{end-3}.LPF)/(GMNA_out.Step{end-2}.l-GMNA_out.Step{end-3}.l));
        if slopeEnd < plat_tol && slopeEndm1 < plat_tol && slopeEndm2 < plat_tol
            termination_cause_GMNAs{I} = strcat(termination_cause_GMNAs{I}," Plateauing response detected.");
        end
    end
    Ncrs(I) = N*GMNA_CritLPFs(I);
    % Write output to file and generate 3D surface of LPFs corresponding to modes
    for key = GMNA_out.Step{end}.BifurcationMode.keys(); GMNA_modes(I,key{1}+1) = GMNA_out.Step{end}.BifurcationMode(key{1}).EigenValue{1}; end

    if isnan(GMNA_CritNModes(I)) || GMNA_CritNModes(I) == 0
        % Copy general GMNA results for cylinders that showcase an axisymmetric failure mode
        termination_cause_GMNAs_n0{I} = termination_cause_GMNAs{I};
        GMNA_CritLPFs_n0(I) = GMNA_CritLPFs(I);
        Ncrs_n0(I) = Ncrs(I);
    else
        % Submit for analysis
        GMNA_n0_out = AQUINAS_Protocol(M1,S1,C1,C2,C3,C4,C5,C6,C7,F1,SOL1,DOF1,A2,O1);
        % Results post-processing
        termination_cause_GMNAs_n0{I} = GMNA_n0_out.Step{end}.termination_cause;
        GMNA_CritLPFs_n0(I) = max(GMNA_n0_out.Step{end}.LPF,GMNA_n0_out.Step{end-1}.LPF); % Considering the case of a limit point, the cylinder's strength should be taken as the highest LPF achieved
        if contains(termination_cause_GMNAs_n0{I},"Analysis aborted.")
            % Check for plateauing region of equilibrium path, and accept result if encountered
            slopeEnd   = abs((GMNA_n0_out.Step{end}.LPF - GMNA_n0_out.Step{end-1}.LPF)/(GMNA_n0_out.Step{end}.l-GMNA_n0_out.Step{end-1}.l));
            slopeEndm1 = abs((GMNA_n0_out.Step{end-1}.LPF - GMNA_n0_out.Step{end-2}.LPF)/(GMNA_n0_out.Step{end-1}.l-GMNA_n0_out.Step{end-2}.l));
            slopeEndm2 = abs((GMNA_n0_out.Step{end-2}.LPF - GMNA_n0_out.Step{end-3}.LPF)/(GMNA_n0_out.Step{end-2}.l-GMNA_n0_out.Step{end-3}.l));
            if slopeEnd < plat_tol && slopeEndm1 < plat_tol && slopeEndm2 < plat_tol
                termination_cause_GMNAs_n0{I} = strcat(termination_cause_GMNAs_n0{I}," Plateauing response detected.");
                GMNA_CritLPFs_n0(I) = GMNA_n0_out.Step{end}.LPF;
            end
        end
        Ncrs_n0(I) = N*GMNA_CritLPFs_n0(I);
    end
    % Submit for analysis
    GNA_out = AQUINAS_Protocol(M1,S1,C1,C2,C3,C4,C5,C6,C7,F2,SOL2,DOF1,A3,O1);
    % Results post-processing
    GNA_CritLPFs(I) = GNA_out.Step{end}.LPF;
    if ~isempty(GNA_out.Step{end}.BuckledIntoMode)
        GNA_CritNModes(I) = GNA_out.Step{end}.BuckledIntoMode;
    end
    Ncrs_gn(I) = Ncls(I)*GNA_CritLPFs(I);
    for key = GNA_out.Step{end}.BifurcationMode.keys(); GNA_modes(I,key{1}+1) = GNA_out.Step{end}.BifurcationMode(key{1}).EigenValue{1}; end

end

%% Save results
% The results are saved at this point, and loaded afterwards, in order to avoid recomputation of the 
save('input_Example_GMNA01.mat','GNA_CritLPFs','GMNA_CritLPFs_n0','GMNA_CritLPFs','GNA_CritNModes','GMNA_CritNModes','rts','Ncrs','Ncrs_n0','Ncrs_gn','Nys','Ncls','termination_cause_GMNAs','termination_cause_GMNAs_n0','GMNA_modes','GNA_modes');

%% Load results
ABAQUS_GMNA_n0 = readtable('input_Example_GMNA01_ABAQUS.csv'); ABAQUS_GMNA_n0 = table2array(ABAQUS_GMNA_n0);
Cylinders_GMNA_struct = load('input_Example_GMNA01.mat');
rts = Cylinders_GMNA_struct.rts;
GMNA_CritLPFs = Cylinders_GMNA_struct.GMNA_CritLPFs;
GMNA_CritNModes = Cylinders_GMNA_struct.GMNA_CritNModes;
GNA_CritNModes = Cylinders_GMNA_struct.GNA_CritNModes;
Ncrs = Cylinders_GMNA_struct.Ncrs;
Nys = Cylinders_GMNA_struct.Nys;
Ncls = Cylinders_GMNA_struct.Ncls;
Ncrs_gn = Cylinders_GMNA_struct.Ncrs_gn;
Ncrs_n0 = Cylinders_GMNA_struct.Ncrs_n0;

%% Post Processing
yielded_indices = find(isnan(GMNA_CritNModes));
buckled_indices = find(~isnan(GMNA_CritNModes));
Ncrs_ab = nan(length(rts),1);
for i=1:length(rts)
    if Ncls(i)>Nys(i); Ncrs_ab(i) = ABAQUS_GMNA_n0(i,2)*Nys(i); else; Ncrs_ab(i) = ABAQUS_GMNA_n0(i,2)*Ncls(i); end
end
RkRpl = Ncrs./Nys;
RkRpl_gn = Ncrs_gn./Nys;
RkRpl_n0 = Ncrs_n0./Nys;
RkRpl_Ab = Ncrs_ab./Nys;
lambda = sqrt(Nys./Ncls);
RkRcr = Ncrs./Ncls;
RkRcr_gn = Ncrs_gn./Ncls;
RkRcr_n0 = Ncrs_n0./Ncls;
RkRcr_Ab = Ncrs_ab./Ncls;
RclRpl = Ncls./Nys;

%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%
%% Capacity curves
figure('Name','Dimensionless resistance of cylinders against relative slenderness and relative reistance','NumberTitle','off','WindowState','Maximized','color','w','PaperOrientation','landscape','PaperType','A3');

% Generalised capacity curve
subplot(1,2,1);
plot([lambda(1) lambda(end)],[1 1],'Color','k','LineStyle',':','LineWidth',LW); grid on; hold on;
plot(lambda, RclRpl,'Color','k','LineStyle','--','LineWidth',LW);
plot(lambda, RkRpl_Ab,'Color',[1.0,0.6471,0.0],'LineStyle','-','LineWidth',0.5*LW,'Marker','s','MarkerSize',1.25*MS,'MarkerEdgeColor',[1.0,0.6471,0.0],'MarkerFaceColor',[1.0,0.6471,0.0]);
plot(lambda, RkRpl,'Color','k','LineStyle','-','LineWidth',LW,'Marker','.','MarkerSize',5*MS,'MarkerEdgeColor','k','MarkerFaceColor','k');
plot(lambda, RkRpl_gn,'Color','r','LineStyle','-','LineWidth',0.5*LW,'Marker','.','MarkerSize',3*MS,'MarkerEdgeColor','r','MarkerFaceColor','r');
plot(lambda, RkRpl_n0,'Color',[0,0.5,1.0000],'LineStyle','-','LineWidth',0.5*LW,'Marker','.','MarkerSize',3*MS,'MarkerEdgeColor',[0,0.5,1.0000],'MarkerFaceColor',[0,0.5,1.0000]);
ax = gca; ax.FontSize = 14;
xlabel('Relative slenderness ','fontsize',FS, 'FontName','Verdana');
xlabel('Relative slenderness $\lambda = \sqrt{R_{s}/R_{cl}}$','interpreter','latex','fontsize',FS, 'FontName','Verdana');
ylabel('Dimensionless resistance $R_{GMNA}$/$R_{s}$','interpreter','latex','fontsize',FS);
xlim([min(lambda) max(lambda)]);
ylim([min(RkRpl) 1.05]);

% Modified capacity curve
subplot(1,2,2); lgnd = cell(6,1);
ph1 = plot([0 1.05],[1 1],'Color','k','LineStyle',':','LineWidth',LW); grid on; hold on;
lgnd{1} = 'Plastic collapse';
ph2 = plot([1 1],[0 1.05],'Color','k','LineStyle','--','LineWidth',LW);
lgnd{2} = 'Classical elastic bifurcation load';
ph5 =  plot(RkRcr_Ab, RkRpl_Ab,'Color',[1.0,0.6471,0.0],'LineStyle','-','LineWidth',0.5*LW,'Marker','s','MarkerSize',1.25*MS,'MarkerEdgeColor',[1.0,0.6471,0.0],'MarkerFaceColor',[1.0,0.6471,0.0]);
lgnd{5} = 'Axisymmetric ABAQUS GMNA';
ph3 = plot(RkRcr, RkRpl,'Color','k','LineStyle','-','LineWidth',LW,'Marker','.','MarkerSize',5*MS,'MarkerEdgeColor','k','MarkerFaceColor','k');
lgnd{3} = 'AQUINAS GMNA';
ph4 = plot(RkRcr_gn, RkRpl_gn,'Color','r','LineStyle','-','LineWidth',0.5*LW,'Marker','.','MarkerSize',3*MS,'MarkerEdgeColor','r','MarkerFaceColor','r');
lgnd{4} = 'AQUINAS GNA';
ph6 = plot(RkRcr_n0, RkRpl_n0,'Color',[0,0.5,1.0000],'LineStyle','-','LineWidth',0.5*LW,'Marker','.','MarkerSize',3*MS,'MarkerEdgeColor',[0,0.5,1.0000],'MarkerFaceColor',[0,0.5,1.0000]);
lgnd{6} = 'AQUINAS GMNA (axisymmetric failure mode)';
isplottedGMNAmodes = nan(size(GMNA_CritNModes)); isplottedGNAmodes = nan(size(GNA_CritNModes));
for I = 1:length(GMNA_CritNModes)
    if ~isnan(GMNA_CritNModes(I)) && GMNA_CritNModes(I)~=GMNA_CritNModes(I-1) && GMNA_CritNModes(I)<36 && isnan(isplottedGMNAmodes(I-1)) && GMNA_CritNModes(I)~=0
        text(RkRcr(I)-0.05,RkRpl(I),num2str(GMNA_CritNModes(I)),'FontSize',10,'Color','k')
    end
end
for I = 2:length(GNA_CritNModes)
    if GNA_CritNModes(I)~=GNA_CritNModes(I-1) && GNA_CritNModes(I)>17 && GNA_CritNModes(I)<36 && isnan(isplottedGNAmodes(I-1))
        text(RkRcr_gn(I)+0.02,RkRpl_gn(I),num2str(GNA_CritNModes(I)),'FontSize',10,'Color','r')
    end
end
ax = gca; ax.FontSize = 14;
xlabel('Relative strength $R_{GMNA}$/$R_{cl}$','interpreter','latex','fontsize',FS);
ylabel('Dimensionless resistance $R_{GMNA}$/$R_{s}$','interpreter','latex','fontsize',FS);
legend([ph1 ph2 ph3 ph4 ph5 ph6],lgnd,'interpreter','latex','fontsize',FSL,'Location','best');
xlim([min(RkRcr) 1.05]);
ylim([min(RkRpl) 1.05]);


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