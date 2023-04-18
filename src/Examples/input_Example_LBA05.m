%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example LBA 05: LBA of closed spherical  shells under external pressure.
% An S3 boundary condition is applied at the bottom apex of the spheres.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% [1] : Hutchinson JW. (2016) - 'Buckling of spherical shells revisited'
% Proc. R. Soc. A 472: 20160577.
% http://dx.doi.org/10.1098/rspa.2016.0577
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.

%%%%%%%%%%%%%%%%%
% Preliminaries %
%%%%%%%%%%%%%%%%%
clearvars, clc, close all force
addpath('..','../AQUINAS_Analysis','../AQUINAS_Objects','../Maple_Source');
warning('off','all')

% Geometry
rt = 50:50:2000; % [-] - Ratios of spherical shell radii over the thickness of the shell wall to be examined
t = 1; % [mm] - Sphere thickness to be considered for all R over t ratios

% Material
MatName = 'Steel'; % Material name
E = 2.0e5; % [N/mm2] - Young's modulus of steel
nu = 0.3; % [-] - Poisson ratio of steel

% Shell segment FE discretization
numelem = 200;

% Parameter to enable the Surrogate Optimisation Algorithm
SurrOptOn = true;

% Number of eigenvalues requested for each circumferential mode
noEigs = 10;

% Plot controls
FS = 20; % [-] - Font Size (pts)
FSL = 18; % [-] - Font Size for Legends (pts)
MS = 8; % [-] - Markser size (pts)
LW = 2; % [-] - Line Width (pts)

pcl =@(r) 2*E/sqrt(3*(1-nu^2))*((t/r)^2); % [N/mm2] - Distributed normal pressure (classical solution, see [1])

%%%%%%%%%%%%%%%%%%%%%%%
% Analytical solution %
%%%%%%%%%%%%%%%%%%%%%%%
ns = round(roots([1 1 -sqrt(12*(1-nu^2))*rt(end)])); ns(ns<0) = nan;
nmax = ceil(1.5*min(ns));
phis = [linspace(pi/2,0,numelem+1)'; linspace(0,-pi/2,numelem+1)']; % Phi angle at each node along the shell's meridian, with zero being at the equator. This will only work if the mesh used for AQUINAS is not graded.
Pnm = cell(nmax,1);

for n = 1:nmax
   Pnm{n} = legendre(n,sin(phis))';
end

%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%
LBA_LPFs = cell(length(rt),1); LBA_n_degrees = cell(length(rt),1); LBA_errors = cell(length(rt),1);
LBA_CritLPFs = zeros(length(rt),1); LBA_CritOrder_m = zeros(length(rt),1); LBA_CritDegree_n = zeros(length(rt),1);
ncrits = zeros(length(rt),1);

for I = 1:length(rt)
    disp(['Computing R/t = ',num2str(rt(I))]);
    ns = round(roots([1 1 -sqrt(12*(1-nu^2))*rt(I)])); ns(ns<0) = nan;
    ncrits(I) = min(ns);

    r = rt(I)*t;

    M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)));

    S1 = AQUINAS_Segment_Object('type','Ellipse','material','Steel','formulation','thin','thickness',t,'rzbot',[r   0],'rztop',[0 r],'geom',[0 0 r r],'els',numelem);
    S2 = AQUINAS_Segment_Object('type','Ellipse','material','Steel','formulation','thin','thickness',t,'rzbot',[0  -r],'rztop',[r 0],'geom',[0 0 r r],'els',numelem);

    C1 = AQUINAS_Constraint_Object('rcoord',0,'zcoord',r,'dofs',{'v'});
    C2 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0,'dofs',{'w'});
    C3 = AQUINAS_Constraint_Object('rcoord',0,'zcoord',-r,'dofs',{'v'});

    P1 = AQUINAS_Distributed_Pressure_Object('segments',{S1,S2},'type','pn','functionHandle',@() -pcl(r));

    if SurrOptOn
        SOL = AQUINAS_Solver_Object('compiler','C++','noThreads',6,'surrogate_optimisation',true,'so_Auto_Bounds',false,'so_circWave_lb',1,'so_circWave_ub',ceil(1.5*ncrits(I)),'constrLagrMultF',1e-3);
    else
        SOL = AQUINAS_Solver_Object('compiler','C++','noThreads',6,'constrLagrMultF',1e-3);
    end

    A1 = AQUINAS_Analysis_Object('type','LBA','circumferentialModes',0:ceil(1.5*ncrits(I)),'noEigenvalues',noEigs,'normalizeEigs','MaxValue');

    O1 = AQUINAS_Output_Mode_Object('mode','VarArgOut');

    LBA_out = AQUINAS_Protocol(M1,S1,S2,C1,C2,C3,P1,SOL,A1,O1);

    % Post-processing
    phis = zeros(LBA_out.Shell_Geom.node{end}(end),1);
    for S=1:length(LBA_out.Shell_Geom.node)
        phis(LBA_out.Shell_Geom.node{S}(1):LBA_out.Shell_Geom.node{S}(end)) = acos(LBA_out.Shell_Geom.z{S}/r)';
    end

    cwKeys = LBA_out.CircWave.keys;
    LBAs = zeros(length(cwKeys),noEigs); LBA_n_degrees{I} = zeros(length(cwKeys),noEigs); LBA_errors{I} = Inf*ones(length(cwKeys),noEigs);
    for CWN = 1:length(cwKeys)
        for EIG = 1:noEigs
            LBAs(CWN,EIG) = LBA_out.CircWave(cwKeys{CWN}).EigenValue{EIG};
            wbar = zeros(LBA_out.Shell_Geom.node{end}(end),1);
            for S=1:length(LBA_out.Shell_Geom.node)
                segIndices = (LBA_out.Shell_Geom.node{S}(1):LBA_out.Shell_Geom.node{S}(end));
                wbar(segIndices) = sin(phis(segIndices)).*LBA_out.CircWave(cwKeys{CWN}).EigenMode{EIG}.u{S} + cos(phis(segIndices)).*LBA_out.CircWave(cwKeys{CWN}).EigenMode{EIG}.w{S};
            end
            for n=max(1,cwKeys{CWN}):nmax
                Wnm = Pnm{n}(:,CWN);
                Wnm = (max(abs(wbar))/max(abs(Wnm)))*Wnm;
                error = min(sum((wbar-Wnm).^2),sum((wbar+Wnm).^2));
                if error < LBA_errors{I}(CWN,EIG)
                    LBA_errors{I}(CWN,EIG) = error; LBA_n_degrees{I}(CWN,EIG) = n; bestFit = Wnm;
                end
            end
        end
    end
    LBAs(LBAs < 0) = NaN; [minLPF,minCircWaveInd] = min(LBAs(:,1));
    LBA_CritLPFs(I) = minLPF; LBA_CritOrder_m(I) = cwKeys{minCircWaveInd}; LBA_CritDegree_n(I) = LBA_n_degrees{I}(minCircWaveInd,1);
    LBA_LPFs{I} = LBAs;

end

%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%
%% Normalized critical pressure
figure('Name','Normalized critical pressure plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(2,1);

h1 = plot([0 rt(end)], [1.0 1.0],...
        'Color',[0.3 0.3 0.3],'LineStyle','--','Marker','none','LineWidth',LW); grid on; hold on;
lgnd{1} = 'Hutchinson(2016)';

h2 = plot(rt, LBA_CritLPFs,...
    'LineStyle','none','Marker','d','MarkerSize',MS,'MarkerEdgeColor','k','MarkerFaceColor','k');
lgnd{2} = 'AQUINAS FE critical LPF solution';

xlim([0 rt(end)]); ylim([0 1.2]);
legend([h1 h2],lgnd,'interpreter','latex','fontsize',FSL,'location','southeast');
xlabel('$\frac{R}{t}$ [-]','interpreter','latex','fontsize',FS);
ylabel('Normalized pressure   $p_{cr}/\frac{2E}{\sqrt{3(1-\nu^{2})}}(\frac{t}{r})^{2}$','interpreter','latex','fontsize',FS);
set(gca,'ticklabelinterpreter','latex','fontsize',FS,'tickdir','out');

%% Degree of associated Legendre polynomial
figure('Name','Degree n of the associated Legendre function Pnm for the axisymmetric mode plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(2,1);

LBA_AxisCritDegrees_n = zeros(length(rt),1);
for I = 1:length(rt)
    LBA_AxisCritDegrees_n(I) = LBA_n_degrees{I}(1,1);
end
h1 = plot(rt, ncrits,...
    'Color',[0.3 0.3 0.3],'LineStyle','--','Marker','none','LineWidth',LW); grid on; hold on;
lgnd{1} = 'Hutchinson(2016) $n_{cr}(n_{cr}+1)=\sqrt{12(1-\nu^{2})}\frac{r}{t}$';

h2 = plot(rt, LBA_AxisCritDegrees_n,...
    'LineStyle','none','Marker','d','MarkerSize',MS,'MarkerEdgeColor','k','MarkerFaceColor','k');
lgnd{2} = 'AQUINAS FE axisymmetric solution';

xlim([0 rt(end)]);
legend([h1 h2],lgnd,'interpreter','latex','fontsize',FSL,'location','southeast');
xlabel('$\frac{r}{t}$ [-]','interpreter','latex','fontsize',FS);
ylabel('Degree $n$ of the associated Legendre function $P_n^m$','interpreter','latex','fontsize',FS);
set(gca,'ticklabelinterpreter','latex','fontsize',FS,'tickdir','out');

%% 3D surface plot (pcr/pcl)
if ~SurrOptOn
    figure('Name','3D surface R/t - m - n with pcr/pcl colorbar plot','NumberTitle','off','WindowState','Maximized','color','w','PaperOrientation','landscape','PaperType','A2');
    lgnd = cell(3,1);

    [RTs,Ms] = meshgrid(rt,0:nmax);
    Ns = nan*RTs; LPFs = nan*RTs;
    for I=1:length(rt)
        for M=1:size(LBA_n_degrees{I},1)
            Ns(M,I) = LBA_n_degrees{I}(M,1);
            LPFs(M,I) = LBA_LPFs{I}(M,1);
        end
    end

    surf(RTs,Ms,Ns,LPFs,'FaceAlpha',0.7,'FaceColor','texturemap'); grid on; hold on;
    lgnd{1} = 'AQUINAS FE solution';
    colormap(gca,'turbo')
    hcb = colorbar('Location','northoutside'); hcbTitle = get(hcb,'Title'); set(hcbTitle,'Interpreter','latex','String','Normalized critical pressure $\frac{p_{cr}}{p_{cl}}$','FontSize',FS)

    plot3(rt,LBA_CritOrder_m,LBA_CritDegree_n,'Color','b','LineStyle','-','LineWidth',LW,'Marker','o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',MS);
    lgnd{2} = 'Critical buckling mode';
    plot3(rt,ncrits,ncrits,'Color','r','LineStyle','-','LineWidth',LW,'Marker','none');
    lgnd{3} = 'Limit of the order $m$ of the associated Legendre function $P_n^m$ ($0\leq m\leq n_{cr}$ with $n_{cr}(n_{cr}+1)=\sqrt{12(1-\nu^{2})}\frac{r}{t}$ (Hutchinson 2016))';

    xlabel(' $\frac{r}{t}$ [-]','interpreter','latex','fontsize',FS);
    ylabel('Order $m$ of the associated Legendre function $P_n^m$','interpreter','latex','fontsize',FS);
    zlabel('Degree $n$ of the associated Legendre function $P_n^m$','interpreter','latex','fontsize',FS);
    legend(lgnd,'interpreter','latex','fontsize',FSL,'location','northwest');
    set(gca, 'ticklabelinterpreter','latex','fontsize',FS,'tickdir','out');
end

%% 3D surface plot (error)
if ~SurrOptOn
    figure('Name','3D surface r/t - m - n with error colorbar plot','NumberTitle','off','WindowState','Maximized');
    lgnd = cell(3,1);

    errors = nan*RTs;
    for I=1:length(rt)
        for M=1:size(LBA_n_degrees{I},1)
            errors(M,I) = LBA_errors{I}(M,1);
        end
    end

    surf(RTs,Ms,Ns,errors,'FaceAlpha',0.7,'FaceColor','interp'); grid on; hold on;
    lgnd{1} = 'AQUINAS FE solution';
    colormap turbo
    hcb = colorbar; hcbTitle = get(hcb,'Title'); set(hcbTitle,'Interpreter','latex','String','Squared error in fitting of associated Legendre function $P_n^m$ ','FontSize',FS)

    plot3(rt,LBA_CritOrder_m,LBA_CritDegree_n,'Color','b','LineStyle','-','LineWidth',LW,'Marker','o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',MS);
    lgnd{2} = 'Critical buckling mode';
    plot3(rt,ncrits,ncrits,'Color','r','LineStyle','-','LineWidth',LW,'Marker','none')
    lgnd{3} = 'Limit of the order $m$ of the associated Legendre function $P_n^m$ ($0\leq m\leq n$)';

    xlabel(' $\frac{r}{t}$ [-]','interpreter','latex','fontsize',FS);
    ylabel('Order of the associated Legendre function $P_n^m$ (circumferential wave number)','interpreter','latex','fontsize',FS);
    zlabel('Degree of the associated Legendre function $P_n^m$','interpreter','latex','fontsize',FS);
    legend(lgnd,'interpreter','latex','fontsize',FSL,'location','northwest');
    set(gca, 'ticklabelinterpreter','latex','fontsize',FS,'tickdir','out');
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