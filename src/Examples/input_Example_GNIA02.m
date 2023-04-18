%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example GNIA 02: GNIA of thin cylindrical shell under internal pressure and
% axial compression with weld depression imperfection at mid-height. A BC1f-S1 boundary
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
% [3] Teng, Jin-Guang, and J. Michael Rotter. "Buckling of pressurized
% axisymmetrically imperfect cylinders under axial loads."
% Journal of engineering mechanics 118.2 (1992): 229-247.
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
h = 3*r; % [mm] - Cylinder height
deltas = [0.5 1.0 1.5]; % [mm] - Weld depression imperfection amplitude at mid-height of cylinder

hw = h/2; % [mm] - Weld depression height

% Material
MatName = 'Steel'; % Material name
E = 2.0e5; % [N/mm2] - Young's modulus of steel
nu = 0.3; % [-] - Poisson ratio of steel

% Loading
Ncl = E*t*t/(r*sqrt(3.0*(1.0-nu*nu))); % [N/mm] - Distributed axisymmetric axial load along top edge (acting vertically downwards)
pbars = [0.0 0.2 0.5 1.0 1.5 2.0 2.5 3.0 4.0 6.0 8.0]; % Dimensionless internal pressure to be considered, given as pbar = p*r/Ncl

% Discretization
nel = 400; % Number of elements along the meridian of the shell

% Weld imperfection funcion
weldtype = 'A';

% Imperfections and dimensionless pressures for which the critical mode will be plotted
pl_d = 1.5;
pl_pbars = [0.2 0.5 1.0 4.0 8.0];
scale = 10;

% Plot controls
FS = 24; % [-] - Font Size (pts)
FSs = 14; % [-] - Font Size for subplots (pts)
MS = 8; % [-] - Markser size (pts)
LW = 3; % [-] - Line Width (pts)


%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%
% Initializations
nmax = 0.5*(12.0*(1.0-nu*nu))^0.25 * sqrt(r/t); % Max. circumferential mode number on the Koiter circle
Nmax = ceil(1.5 * nmax);
lambda = pi*sqrt(r*t)/((3*(1-nu^2))^0.25);
crit_circ_modes = nan(length(deltas),length(pbars));
crit_LPFs = nan(length(deltas),length(pbars));
eigenmodes = cell(length(deltas),length(pbars));
geometries = cell(length(deltas),1);

for I = 1:length(deltas)
    % Definition of geometry and boundary conditions, does not depend on the dimeensionless pressures
    Zb = linspace(0,hw,ceil(nel/2)+1); Zt = linspace(hw,h,ceil(nel/2)+1);
    if strcmp(weldtype,'A'); k = 1; elseif strcmp(weldtype,'B'); k = 0; end
    Rb = r - deltas(I)*exp(-pi*abs(Zb-hw)/lambda).*(cos(pi*abs(Zb-hw)/lambda) + k*sin(pi*abs(Zb-hw)/lambda));
    Rt = r - deltas(I)*exp(-pi*abs(Zt-hw)/lambda).*(cos(pi*abs(Zt-hw)/lambda) + k*sin(pi*abs(Zt-hw)/lambda));
    Rb = flipud(Rb'); Zb = flipud(Zb'); Rt = flipud(Rt'); Zt = flipud(Zt'); % Top-down order for R-Z point array
    geometries{I} = [[Rt;Rb],[Zt;Zb]];

    M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)));

    S1 = AQUINAS_Segment_Object('type','Interp','material',MatName,'thickness',t,'rzpts',[Rt, Zt],'els',nel/2);
    S2 = AQUINAS_Segment_Object('type','Interp','material',MatName,'thickness',t,'rzpts',[Rb, Zb],'els',nel/2);

    F1 = AQUINAS_Nodal_Force_Object('type','Fw','rcoord',Rt(1),'zcoord',Zt(1),'magnitude', -Ncl);

    A1 = AQUINAS_Analysis_Object('type','GNA','circumferentialModes',0:Nmax,'normalizeEigs','MaxValue');

    SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',6,'constrLagrMultF',0.00001,'dLPF',0.02,'dLPFmax',0.02,'maxAttempts',10,...
        'terminationConditions','AC','noMaxSteps',2000,'Jmax',100,'simultaneousBifurcationTreatment','pickSmallestWithinRange','ksi',1e-4,...
        'eigsNegativeTreatment','Drop','eigsSigma',1.0,'Jd',5);

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

    for J = 1:length(pbars)
        disp(strcat("Computing cylinder with weld imperfection amplitude delta = ",num2str(deltas(I))," and a dimensionless pressure pbar = ",num2str(pbars(J))));
        p = pbars(J)*Ncl/r;

        % Create pressure object based on the dimensionless pressure
        P1 = AQUINAS_Distributed_Pressure_Object('segments',{S1,S2},'type','pn','functionHandle',@() p);

        GNIA_out = AQUINAS_Protocol(M1,S1,S2,C1,C2,C3,C4,C5,C6,C7,C8,F1,P1,SOL1,A1,O1);

        % Post-processing
        if isempty(GNIA_out.Step{end}.BuckledIntoMode)
            warning(strcat("'AQUINAS Warning: The cylinder has not buckled for an imperfection amplitude delta = ",num2str(deltas(I))," and a dimensionless pressure pbar = ",num2str(pbars(J))));
        else
            crit_LPFs(I,J) = GNIA_out.Step{end}.LPF;
            crit_circ_modes(I,J) = GNIA_out.Step{end}.BuckledIntoMode;
            eigenmodes{I,J}.nwave = containers.Map('KeyType','int32','ValueType','any');
            for key = GNIA_out.Step{end}.BifurcationMode.keys
                eigenmodes{I,J}.nwave(key{1}) = struct('u',[GNIA_out.Step{end}.BifurcationMode(key{1}).EigenMode{1}.u{1};GNIA_out.Step{end}.BifurcationMode(key{1}).EigenMode{1}.u{2}],'w',[GNIA_out.Step{end}.BifurcationMode(key{1}).EigenMode{1}.w{1};GNIA_out.Step{end}.BifurcationMode(key{1}).EigenMode{1}.w{2}]);
            end
            disp(strcat("For a cylinder with weld imperfection amplitude delta = ",num2str(deltas(I))," and a dimensionless pressure pbar = ",num2str(pbars(J)),", the critical LPF is ",num2str(crit_LPFs(I,J))," and the critical circumferential wavenumber is : ",num2str(crit_circ_modes(I,J))));
        end

    end
end

%%%%%%%%%%%%%
% Load Data %
%%%%%%%%%%%%%
%% Load data digitised data from [2]
Data = readtable('input_Example_GNIA02.csv'); Data = table2array(Data);
RotterTeng_deltas = [0.0 1.0 1.5];
RotterTeng_pbars = [Data(:,1),Data(:,3),Data(:,5)];
RotterTeng_NcrNclRatios = [Data(:,2),Data(:,4),Data(:,6)];

%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%
%% Critical LPFs plot
figure('Name','Effect of internal pressures on buckling strengths','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(4,1);

% Rotter and Teng (1989)
ph1 = plot(RotterTeng_pbars(:,1), RotterTeng_NcrNclRatios(:,1),...
    'LineStyle','--','LineWidth',LW,'Color',[0.6 0.6 0.6],'Marker','x','MarkerSize',MS,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6]); grid on; hold on;
lgnd{1} = 'Rotter and Teng (1989) $\delta_0/t = 0.5$ [Digitised]';
ph2 = plot(RotterTeng_pbars(:,2), RotterTeng_NcrNclRatios(:,2),...
    'LineStyle','--','LineWidth',LW,'Color',[0.6 0.6 0.6],'Marker','^','MarkerSize',MS,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor','w');
lgnd{2} = 'Rotter and Teng (1989) $\delta_0/t = 1.0$ [Digitised]';
ph3 = plot(RotterTeng_pbars(:,3), RotterTeng_NcrNclRatios(:,3),...
    'LineStyle','--','LineWidth',LW,'Color',[0.6 0.6 0.6],'Marker','o','MarkerSize',MS,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor','w');
lgnd{3} = 'Rotter and Teng (1989) $\delta_0/t = 1.5$ [Digitised]';

% AQUINAS solution
ph4 = plot(crit_LPFs(1,:).*pbars, crit_LPFs(1,:),...
    'LineStyle','-','Color','k','LineWidth',LW,'Marker','x','MarkerSize',MS,'MarkerEdgeColor','k','MarkerFaceColor','w');
lgnd{4} = 'Aquinas $\delta_0/t = 0.5$';
ph5 = plot(crit_LPFs(2,:).*pbars, crit_LPFs(2,:),...
    'LineStyle','-','Color','k','LineWidth',LW,'Marker','^','MarkerSize',MS,'MarkerEdgeColor','k','MarkerFaceColor','w');
lgnd{5} = 'Aquinas $\delta_0/t = 1.0$';
ph6 = plot(crit_LPFs(3,:).*pbars, crit_LPFs(3,:),...
    'LineStyle','-','Color','k','LineWidth',LW,'Marker','o','MarkerSize',MS,'MarkerEdgeColor','k','MarkerFaceColor','w');
lgnd{6} = 'Aquinas $\delta_0/t = 1.5$';

xlim([pbars(1) pbars(end)]);
legend([ph1,ph2,ph3,ph4,ph5,ph6],lgnd,'interpreter','latex','fontsize',FS,'location','best');
xlabel('$\bar{p} = pr/\sigma_{cl}t$ [-]','interpreter','latex','fontsize',FS);
ylabel('Buckling strength $\sigma_{cr}/\sigma_{cl}$ [-]','interpreter','latex','fontsize',FS);
set(gca,'ticklabelinterpreter','latex','fontsize',FS,'tickdir','out');

%% Selected modes plot
figure('Name','Critical eigenmodes for selected imperfections','NumberTitle','off','WindowState','Maximized','color','w','PaperOrientation','landscape','PaperType','A3');
ind_d = find(deltas == pl_d);
for i = 1:length(pl_pbars)
    ind = find(pbars==pl_pbars(i));
    subplot(1,length(pl_pbars),i);
    hold on
    plot(geometries{ind_d}(:,1),geometries{ind_d}(:,2),'k--')
    plot(geometries{ind_d}(:,1)+scale*eigenmodes{ind_d,ind}.nwave(crit_circ_modes(ind_d,ind)).u,geometries{ind_d}(:,2)+scale*eigenmodes{ind_d,ind}.nwave(crit_circ_modes(ind_d,ind)).w,'k-')
    xlabel('Radial coordinate R [mm]','Interpreter','latex','FontSize',FSs);
    if i == 1; ylabel('Axial coordinate Z [mm]','interpreter','latex','fontsize',FS); end
    title(strcat("Critical mode ",num2str(crit_circ_modes(ind_d,ind))," for $\bar{p}$ = ",num2str(pbars(ind))),'Interpreter','latex','FontSize',FSs);
    set(gca,'XColor','none','Ycolor','none');
    axis tight;
    hold off
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