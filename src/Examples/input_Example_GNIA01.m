%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example GNIA 01: GNIA of unpressurised thin cylindrical shell under axial
% compression with weld depression imperfection at mid-height. A BC1f-S1 boundary
% is considered at the base of the cylinder and a BC2f-S3 at the top.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% [1] Sadowski A.J., Pototschnig L. & Constantinou P. (2018) "The 'panel
% analysis' technique in the computational study of axisymmetric
% thin-walled shell systems" Finite Elements in Analysis and Design, 152,
% 55-68.
% https://doi.org/10.1016/j.finel.2018.07.004
% (The above contains many references to classical texts such as those of
% Koiter which are the origin of the Koiter circle theory used here.)
% [2] J.M. Rotter, J.G. Teng - "Elastic stability of cylindrical
% shells with weld depressions" ASCE J. Struct. Eng., 115 (5) (1989), pp. 1244-1263
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.

%%%%%%%%%%%%%%%%%
% Preliminaries %
%%%%%%%%%%%%%%%%%
clearvars, clc, close all
addpath('..','../AQUINAS_Analysis','../AQUINAS_Objects','../Maple_Source');
warning('off','all') % warnings related to complex eigenvalues can be rather distracting in the command window

%%%%%%%%%%
% Inputs %
%%%%%%%%%%
% Geometry
r = 1000.0; % [mm] - Cylinder radius
t = 1.0; % [mm] - Cylinder thickness
h = 3*r; % [mm] - Cylinder height
deltas = [0.0 0.1 0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0 2.5 3.0 3.5 4.0 5.0 6.0 8.0 10.0]; % [mm] - Weld depression imperfection amplitude at mid-height of cylinder

hw = h/2; % [mm] - Weld depression height

% Material
MatName = 'Steel'; % Material name
E = 2.0e5; % [N/mm2] - Young's modulus of steel
nu = 0.3; % [-] - Poisson ratio of steel

% Loading
Ncl = E*t*t/(r*sqrt(3.0*(1.0-nu*nu))); % [N/mm] - Distributed axisymmetric axial load along top edge (acting vertically downwards)

% Discretization
nel = 400; % Number of elements along the meridian of the shell

% Weld imperfection funcion
weldtype = 'A';

% Imperfections for which the critical mode will be plotted
pl_ds = [0.0 0.5 1.0 1.5 2.0];
scale = 10;

% Plot controls
FS = 24; % [-] - Font Size (pts)
FSs = 20; % [-] - Font Size for sublot titles/axes (pts)
MS = 8; % [-] - Markser size (pts)
LW = 3; % [-] - Line Width (pts)


%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%
% Initializations
nmax = 0.5*(12.0*(1.0-nu*nu))^0.25 * sqrt(r/t); % Max. circumferential mode number on the Koiter circle
Nmax = ceil(1.5 * nmax);
lambda = pi*sqrt(r*t)/((3*(1-nu^2))^0.25);
crit_circ_modes = nan(size(deltas));
crit_LPFs = nan(size(deltas));
eigenmodes = cell(length(deltas),1);
geometries = cell(length(deltas),1);

for I = 1:length(deltas)
    disp(strcat("Computing cylinder with weld imperfection amplitude delta = ",num2str(deltas(I))));
    Zb = linspace(0,hw,ceil(nel/2)+1); Zt = linspace(hw,h,ceil(nel/2)+1);
    if strcmp(weldtype,'A'); k = 1; elseif strcmp(weldtype,'B'); k = 0; end
    Rb = r - deltas(I)*exp(-pi*abs(Zb-hw)/lambda).*(cos(pi*abs(Zb-hw)/lambda) + k*sin(pi*abs(Zb-hw)/lambda));
    Rt = r - deltas(I)*exp(-pi*abs(Zt-hw)/lambda).*(cos(pi*abs(Zt-hw)/lambda) + k*sin(pi*abs(Zt-hw)/lambda));
    Rb = flipud(Rb'); Zb = flipud(Zb'); Rt = flipud(Rt'); Zt = flipud(Zt'); % Top-down order for R-Z point array

    M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)));

    S1 = AQUINAS_Segment_Object('type','Interp','material',MatName,'thickness',t,'rzpts',[Rt, Zt],'els',nel/2);
    S2 = AQUINAS_Segment_Object('type','Interp','material',MatName,'thickness',t,'rzpts',[Rb, Zb],'els',nel/2);
    S3 = AQUINAS_Segment_Object('type','Cone','material',MatName,'thickness',t,'rzbot',[r 0.0],'rztop',[r h],'els',nel);

    F1 = AQUINAS_Nodal_Force_Object('type','Fw','rcoord',Rt(1),'zcoord',Zt(1),'magnitude', -Ncl);

    A1 = AQUINAS_Analysis_Object('type','GNA','circumferentialModes',0:Nmax,'normalizeEigs','MaxValue');

    SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',6,'constrLagrMultF',0.00001,'dLPF',0.01,'dLPFmax',0.01,'maxAttempts',10,...
        'terminationConditions','AC','noMaxSteps',2000,'Jmax',100,'simultaneousBifurcationTreatment','pickSmallestWithinRange','ksi',1e-4);

    O1 = AQUINAS_Output_Mode_Object('mode','VarArgOut');

    % Boundary conditions at base of cylinder
    C1 = AQUINAS_Constraint_Object('rcoord',Rb(end),'zcoord',Zb(end),'dofs',{'u'},'state','buckling');
    C2 = AQUINAS_Constraint_Object('rcoord',Rb(end),'zcoord',Zb(end),'dofs',{'v'});
    C3 = AQUINAS_Constraint_Object('rcoord',Rb(end),'zcoord',Zb(end),'dofs',{'w'},'state','prebuckling');
    C4 = AQUINAS_Constraint_Object('rcoord',Rb(end),'zcoord',Zb(end),'dofs',{'b'},'state','prebuckling');
    % Boundary conditions at weld position
    C5 = AQUINAS_Constraint_Object('rcoord',Rt(end),'zcoord',Zt(end),'dofs',{'w'},'state','buckling');
    % Boundary conditions at top of cylinder
    C6 = AQUINAS_Constraint_Object('rcoord',Rt(1),'zcoord',Zt(1),'dofs',{'u'},'state','buckling');
    C7 = AQUINAS_Constraint_Object('rcoord',Rt(1),'zcoord',Zt(1),'dofs',{'v'});
    C8 = AQUINAS_Constraint_Object('rcoord',Rt(1),'zcoord',Zt(1),'dofs',{'b'},'state','prebuckling');
    % Boundary condition for axial displacement at base of cylinder, only for perfect model
    C9 = AQUINAS_Constraint_Object('rcoord',Rb(end),'zcoord',Zb(end),'dofs',{'w'});

    % Analysis
    if deltas(I) > 1e-3
        GNIA_out = AQUINAS_Protocol(M1,S1,S2,C1,C2,C3,C4,C5,C6,C7,C8,F1,SOL1,A1,O1);
    else
        GNIA_out = AQUINAS_Protocol(M1,S1,S2,C1,C2,C4,C6,C7,C8,C9,F1,SOL1,A1,O1);
    end

    % Post-processing
    geometries{I} = [[Rt;Rb],[Zt;Zb]];
    if isempty(GNIA_out.Step{end}.BuckledIntoMode)
        warning(strcat("'AQUINAS Warning: The cylinder has not buckled for an imperfection amplitude of ",num2str(deltas(I))));
    else
        crit_LPFs(I) = GNIA_out.Step{end}.LPF;
        crit_circ_modes(I) = GNIA_out.Step{end}.BuckledIntoMode;
        eigenmodes{I}.nwave = containers.Map('KeyType','int32','ValueType','any');
        for key = GNIA_out.Step{end}.BifurcationMode.keys
            eigenmodes{I}.nwave(key{1}) = struct('u',[GNIA_out.Step{end}.BifurcationMode(key{1}).EigenMode{1}.u{1};GNIA_out.Step{end}.BifurcationMode(key{1}).EigenMode{1}.u{2}],'w',[GNIA_out.Step{end}.BifurcationMode(key{1}).EigenMode{1}.w{1};GNIA_out.Step{end}.BifurcationMode(key{1}).EigenMode{1}.w{2}]);
        end
        disp(strcat("For a cylinder with weld imperfection amplitude delta = ",num2str(deltas(I)),", the critical LPF is ",num2str(crit_LPFs(I))," and the critical circumferential wavenumber is : ",num2str(crit_circ_modes(I))));
    end

end

%%%%%%%%%%%%%
% Load Data %
%%%%%%%%%%%%%
%% Load data generated/digitised by Sadowski et al. in [1] and data digitised from the corresponding Figure in [2]
Data = readtable('input_Example_GNIA01.csv'); Data = table2array(Data);
RotterTeng_delta = Data(:,1);
RotterTeng_n = Data(:,2);
RotterTeng_LPF = Data(:,3);
Panel_ABA_delta = Data(:,4);
Panel_ABA_n = Data(:,5);
Panel_ABA_LPF = Data(:,6);
BlueBook_delta = Data(:,7);
BlueBook_n = Data(:,8);
BlueBook_LPF = Data(:,9);


%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%
%% Critical LPFs plot
figure('Name','Critical LPF per imperfection amplitude','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(4,1);

% Rotter and Teng (1989)
ph1 = plot(RotterTeng_delta, RotterTeng_LPF,...
    'LineStyle','--','LineWidth',LW,'Color','r','Marker','o','MarkerSize',1.3*MS,'MarkerEdgeColor','r','MarkerFaceColor','r'); grid on; hold on;
lgnd{1} = 'Rotter and Teng (1989) FE [Digitised]';

% Panel technique
ph2 = plot(Panel_ABA_delta, Panel_ABA_LPF,...
    'LineStyle','--','LineWidth',LW,'Color',[0.6 0.1 0.6],'Marker','o','MarkerSize',1.3*MS,'MarkerEdgeColor',[0.6 0.1 0.6],'MarkerFaceColor',[0.6 0.1 0.6]);
lgnd{2} = '3D FE - Panel Analysis';

% Rotter and Teng (2004)
ph3 = plot(BlueBook_delta, BlueBook_LPF,...
    'LineStyle','--','LineWidth',LW,'Color',[0.6 0.6 1.0],'Marker','o','MarkerSize',1.3*MS,'MarkerEdgeColor',[0.6 0.6 1.0],'MarkerFaceColor',[0.6 0.6 1.0]);
lgnd{3} = 'Rotter and Teng (2004) [Digitised]';

% AQUINAS solution
ph4 = plot(deltas, crit_LPFs,...
    'LineStyle','-','Color','k','LineWidth',LW,'Marker','d','MarkerSize',MS,'MarkerEdgeColor','k','MarkerFaceColor','w');
lgnd{4} = 'Aquinas';

xlim([deltas(1) deltas(end)]);
legend([ph1,ph2,ph3,ph4],lgnd,'interpreter','latex','fontsize',FS,'location','best');
xlabel('Depression amplitude $\frac{\delta_{imp}}{t}$ [mm]','interpreter','latex','fontsize',FS);
ylabel('Normalised buckling load $\frac{N_{GNIA}}{N_{cl}}$ [-]','interpreter','latex','fontsize',FS);
set(gca,'ticklabelinterpreter','latex','fontsize',FS,'tickdir','out');

%% Critical circumferential wavenumbers plot
figure('Name','Critical circumferential wavenumber per imperfection amplitude','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(4,1);

% Rotter and Teng (1989)
ph1 = plot(RotterTeng_delta, RotterTeng_n,...
    'LineStyle','--','LineWidth',LW,'Color','r','Marker','o','MarkerSize',1.3*MS,'MarkerEdgeColor','r','MarkerFaceColor','r'); grid on; hold on;
lgnd{1} = 'Rotter and Teng (1989) FE [Digitised]';

% Panel technique
ph2 = plot(Panel_ABA_delta, Panel_ABA_n,...
    'LineStyle','--','LineWidth',LW,'Color',[0.6 0.1 0.6],'Marker','o','MarkerSize',1.3*MS,'MarkerEdgeColor',[0.6 0.1 0.6],'MarkerFaceColor',[0.6 0.1 0.6]);
lgnd{2} = '3D FE - Panel Analysis';

% Rotter and Teng (2004)
ph3 = plot(BlueBook_delta, BlueBook_n,...
    'LineStyle','--','LineWidth',LW,'Color',[0.6 0.6 1.0],'Marker','o','MarkerSize',1.3*MS,'MarkerEdgeColor',[0.6 0.6 1.0],'MarkerFaceColor',[0.6 0.6 1.0]);
lgnd{3} = 'Rotter and Teng (2004) [Digitised]';

% AQUINAS solution
ph4 = plot(deltas, crit_circ_modes,...
    'LineStyle','-','Color','k','LineWidth',LW,'Marker','d','MarkerSize',MS,'MarkerEdgeColor','k','MarkerFaceColor','w');
lgnd{4} = 'Aquinas';

xlim([deltas(1) deltas(end)]);
legend([ph1,ph2,ph3,ph4],lgnd,'interpreter','latex','fontsize',FS,'location','best');
xlabel('Depression amplitude $\frac{\delta_{imp}}{t}$ [mm]','interpreter','latex','fontsize',FS);
ylabel('No. of circumferential full waves n','interpreter','latex','fontsize',FS);
set(gca,'ticklabelinterpreter','latex','fontsize',FS,'tickdir','out');

%% Critical modes plot
figure('Name','Critical eigenmodes for selected imperfections','NumberTitle','off','WindowState','Maximized','color','w','PaperType','A3','PaperOrientation','landscape');
ds = zeros(size(pl_ds));
for i = 1:length(pl_ds)
    ind = find(deltas==pl_ds(i));
    subplot(1,length(pl_ds),i); hold on
    plot(geometries{ind}(:,1),geometries{ind}(:,2),'k--')
    plot(geometries{ind}(:,1)+scale*eigenmodes{ind}.nwave(crit_circ_modes(ind)).u,geometries{ind}(:,2)+scale*eigenmodes{ind}.nwave(crit_circ_modes(ind)).w,'k-');
    title(strcat("Critical mode ",num2str(crit_circ_modes(ind))," for $\delta_0$ = ",num2str(deltas(ind))),'Interpreter','latex','FontSize',FSs);
    set(gca,'XColor','none','Ycolor','none');
    axis tight; hold off
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