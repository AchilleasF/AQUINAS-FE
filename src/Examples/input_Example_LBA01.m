%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example LBA 01: LBA of unpressurised thin cylindrical shells of varying
% length under uniform axial compression. A BC1f / S1 condition is
% assumed at the base, and a BC2f / S3 condition at the top.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% [1] Sadowski A.J., Pototschnig L. & Constantinou P. (2018) "The 'panel
% analysis' technique in the computational study of axisymmetric
% thin-walled shell systems" Finite Elements in Analysis and Design, 152,
% 55-68.
% https://doi.org/10.1016/j.finel.2018.07.004
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.

%%%%%%%%%%%%%%%%%
% Preliminaries %
%%%%%%%%%%%%%%%%%
clearvars, clc, close all
addpath('..','../AQUINAS_Analysis','../AQUINAS_Objects','../Maple_Source');
warning('off','all') % Supress warnings in order to avoid messages relating to complex eigenvalues

%%%%%%%%%%
% Inputs %
%%%%%%%%%%
% Geometry
r = 500.0; % [mm] - Cylinder radius
t = 1.0; % [mm] - Cylinder thickness
OMEGAS = [0.001:0.0005:0.009 0.01:0.005:0.09 0.1:0.05:0.9 1:0.1:10]; % [-] - Array of dimensionless cylinder lengths

% Material
MatName = 'Steel'; % Material name
E = 2.0e5; % [N/mm2] - Young's modulus of steel
nu = 0.3; % [-] - Poisson ratio of steel

% Loading
Ncl = E*t*t/(r*sqrt(3.0*(1.0-nu*nu))); % [N/mm] - Distributed axisymmetric axial load along top edge (acting vertically downwards)

% Data generation option
enableSurOpt = false; % If enableSurOpt is true then the Surrogate Optimisation capability of AQUINAS will be employed in order to identify the critical
                      % circumferential wavenumber that the shell will buckle into, as well as the corresponding LPF. If enableSurOpt is set to false then
                      % there will be no Surrogate Optimisation, and instead the critical circumferential wavenumber will be identified through a brute force
                      % approach of evaluating critical pressures for a range of circumferential wavenumber. The critical wavenumber will then be found as the
                      % one which corresponds to the minimum buckling pressure.
                      % If enableSurOpt is set to false then an additional 3D surface plot of the critical pressures for the various Batdorf parameters
                      % and circumferential wavenumbers considered is presented. This surface is not shown if the Surrogate Optimisation algorithm is used (enableSurOpt = false),
                      % since only a selection of circumferential wavenumbers are actually evaluated, for each Batdorf parameter Z.

% Plot controls
FS = 24; % [-] - Font Size (pts)
FSL = 20; % [-] - Font Size for Legends (pts)
MS = 5; % [-] - Markser size (pts)
LW = 3; % [-] - Line Width (pts)


%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%
% Reference [1] contains many references to classical texts such as those of
% Koiter which are the origin of the Koiter circle theory used here.
nmax = 0.5*(12.0*(1.0-nu*nu))^0.25 * sqrt(r/t); % Max. circumferential mode number on the Koiter circle
Nmax = ceil(1.5 * nmax);
if ~enableSurOpt
    LBA_LPFs = zeros(length(OMEGAS),Nmax+1);
end
LBA_CritLPFs = zeros(length(OMEGAS),1);
LBA_CritNModes = zeros(length(OMEGAS),1);

for O = 1:length(OMEGAS)
    OMEGA = OMEGAS(O);
    L = OMEGA*r*sqrt(r/t);
    disp(['Computing Omega = ',num2str(OMEGA)]);

    M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)));

    S1 = AQUINAS_Segment_Object('type','Cone','material',MatName,'thickness',t,'rzbot',[r 0.0],'rztop',[r L],'els',200);

    % BC1f condition at base of cylinder
    C1 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'u'});
    C2 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'v'});
    C3 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'w'});

    % BC2f condition at top of cylinder
    C4 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',L,'dofs',{'u'});
    C5 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',L,'dofs',{'v'});

    F1 = AQUINAS_Nodal_Force_Object('type','Fw','rcoord',r,'zcoord',L,'magnitude',-Ncl);

    O1 = AQUINAS_Output_Mode_Object('mode','VarArgOut');

    if enableSurOpt

        SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',6,'surrogate_optimisation',true,'so_Auto_Bounds',false,'so_circWave_lb',1,'so_circWave_ub',Nmax,'constrLagrMultF',0.001);

        A1 = AQUINAS_Analysis_Object('type','LBA','noEigenvalues',10);

        LBA_out = AQUINAS_Protocol(M1,S1,C1,C2,C3,C4,C5,F1,SOL1,A1,O1);

        % Post-processing
        ModeKeys = LBA_out.CircWave.keys;
        critLPF = Inf; critMode = nan;
        for I = 1:length(ModeKeys)
            if LBA_out.CircWave(ModeKeys{I}).EigenValue{1} < critLPF
                critLPF = LBA_out.CircWave(ModeKeys{I}).EigenValue{1};
                critMode = ModeKeys{I};
            end
        end
        if critLPF <= 0; critLPF = nan; end
        LBA_CritLPFs(O) = critLPF; LBA_CritNModes(O) = critMode;

    else

        SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',6,'constrLagrMultF',0.001);

        A1 = AQUINAS_Analysis_Object('type','LBA','circumferentialModes',0:Nmax,'noEigenvalues',10);

        LBA_out = AQUINAS_Protocol(M1,S1,C1,C2,C3,C4,C5,F1,SOL1,A1,O1);

        % Post-processing
        LBAs = zeros(1,Nmax+1); I = 0;
        for N = 0:Nmax
            I = I + 1; LBAs(I) = LBA_out.CircWave(N).EigenValue{1};
        end
        LBAs(LBAs < 0) = NaN; [m,i] = min(LBAs);
        LBA_LPFs(O,:) = LBAs;
        LBA_CritLPFs(O) = m; LBA_CritNModes(O) = i-1;

    end
end


%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%
%% Dimensionless buckling load vs dimensionless length plot (2D)
figure('Name','Dimensionless buckling load vs dimensionless length','NumberTitle','off','WindowState','Maximized','color','w','PaperOrientation','landscape','PaperType','A3');
lgnd = cell(3,1);

% Plate buckling assuming a flat plate of width 2*pi*r, length l and
% thickness t under uniform in-plane uni-axial compression
NclP = @(l) E*pi*pi*t*t*t./(12.0*l.*l*(1.0-nu*nu));
hs = (0.001:0.001:0.009)*r*sqrt(r/t);
plot(0.001:0.001:0.009, NclP(hs)/Ncl,'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',LW); grid on; hold on;
lgnd{1} = 'Plate-like buckling';

% Euler column buckling assuming a CHS with an effective length of 0.7*l,
% appropriate for BC1f - BC2f (one end fixed, one end vertical roller)
NclE = @(l) E*pi*pi*r*r*t./(2.0*(0.7*0.7*l.*l));
hs = (1.5:10)*r*sqrt(r/t);
plot(1.5:10, NclE(hs)/Ncl,'Color',[0.5 0.5 0.5],'LineStyle','-','LineWidth',LW); grid on; hold on;
lgnd{2} = 'Euler column buckling';

% AQUINAS solution
plot(OMEGAS, LBA_CritLPFs,'Color','k','LineStyle','-','LineWidth',LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',MS);
lgnd{3} = 'AQUINAS FE solution';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Dimensionless length $\Omega = \frac{L}{r}\sqrt{\frac{t}{r}}$','interpreter','latex','fontsize',FS);
ylabel('Dimensionless critical buckling load $\frac{N}{N_{cl}}$','interpreter','latex','fontsize',FS);
set(gca,'XScale','log','YScale','log','ticklabelinterpreter','latex','fontsize',FS,'tickdir','out');

%% Dimensionless buckling load vs dimensionless length plot (3D)
if ~enableSurOpt
    figure('units','normalized','outerposition',[0 0 1 1],'WindowState','Maximized','color','w');

    [O,N] = meshgrid(OMEGAS, 0:Nmax);
    surf(O,N,LBA_LPFs','FaceAlpha',0.5,'FaceColor',[0.8 0.8 0.8]); shading faceted; grid on; hold on;

    [O,N] = meshgrid(0.001:0.001:0.009, 0:Nmax); hs = O*r*sqrt(r/t);
    mesh(O,N,NclP(hs)/Ncl,'FaceAlpha',0.5,'FaceColor',[0.8 0.8 0.8]);

    [O,N] = meshgrid(1.5:10, 0:Nmax); hs = O*r*sqrt(r/t);
    mesh(O,N,NclE(hs)/Ncl,'FaceAlpha',0.5,'FaceColor',[0.8 0.8 0.8]);

    plot3(OMEGAS,LBA_CritNModes,LBA_CritLPFs,'Color','k','LineStyle','-','LineWidth',LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',MS);
    xlabel('Dimensionless length $\Omega = \frac{l}{r}\sqrt{\frac{t}{r}}$','interpreter','latex','fontsize',FS);
    ylabel('Circumferential wave number $n$','interpreter','latex','fontsize',FS);
    zlabel('Dimensionless critical buckling load $\frac{N}{N_{cl}}$','interpreter','latex','fontsize',FS);
    set(gca,'XScale','log','ZScale','log','Zlim',[0,10],'ticklabelinterpreter','latex','fontsize',FS,'tickdir','out');
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