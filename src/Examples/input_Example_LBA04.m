%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example LBA 04: LBA of a truncated conical shell under the effect of
% meridional compression.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% [1] Spagnoli A. (2003). "Koiter circles in the buckling of axially compressed conical shells"
% International Journal of Solids and Structures, Volume 40, Pages 6095-6109
% https://doi.org/10.1016/S0020-7683(03)00369-X
% [2] Sadowski A.J. (2019) "On the advantages of hybrid beam-shell structural finite
% element models for the efficient analysis of metal wind turbine support towers"
% Finite Elements in Analysis and Design, Volume 162, Pages 19-33
% https://doi.org/10.1016/j.finel.2019.05.002.
% [3] Sadowski A.J., Pototschnig L. & Constantinou P. (2018) "The 'panel
% analysis' technique in the computational study of axisymmetric
% thin-walled shell systems" Finite Elements in Analysis and Design, 152,
% 55-68.
% https://doi.org/10.1016/j.finel.2018.07.004
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.

%%%%%%%%%%%%%%%%%
% Preliminaries %
%%%%%%%%%%%%%%%%%
clearvars, clc, close all force
addpath('..','../AQUINAS_Analysis','../AQUINAS_Objects','../Maple_Source');

%%%%%%%%%%
% Inputs %
%%%%%%%%%%
% Geometry
r1 = 1000.0; % [mm] - Cone top edge radius
alpha = 15*pi/180; % [rad] - Tapering angle (apex half angle)
t = 1.0; % [mm] - Cone thickness
L = 3000.0; % [mm] - Length of cone's meridian
r2 = r1 + L*sin(alpha);
z1 = L*cos(alpha);
z2 = 0.0;

% Material
MatName = 'Steel'; % Material name
E = 2.0e5; % [N/mm2] - Young's modulus of steel
nu = 0.3; % [-] - Poisson ratio of steel

% Loading
Ncone = -E*t*t*(cos(alpha)^2)/(r1*sqrt(3.0*(1.0-nu*nu))); % [N/mm] - Distributed axisymmetric axial load along top edge (positive acting vertically upwards)

% Tolerance for the inclusion of eigenvalue/eigenmode as sufficiently close to the classical reference load
lim = 0.02;


% Plot controls
FS = 20; % [-] - Font Size (pts)
FSL = 18; % [-] - Font Size for Legend (pts)
MS = 4; % [-] - Markser size (pts)
LW = 2; % [-] - Line Width (pts)


%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%
nmax = 32; % Max. circumferential wave number on the Koiter ellipse
mmax = ceil(sqrt(2)*sqrt(sqrt(3*(1-nu^2)))*L/(pi*sqrt(t)*sqrt((r1+r2)/2/cos(alpha)))); % Max. meridional half-wave number on the Koiter circle
Nmax = nmax; Mmax = ceil(1.5*mmax);
Emax = 20;
LBA_LPF = nan(Mmax,Nmax+1);

M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)));

% S1 = AQUINAS_Segment_Object('type','Cone','material',MatName,'thickness',t,'rzbot',[r2 Z2],'rztop',[r1 Z1],'els',10*Mmax);
S1 = AQUINAS_Segment_Object('type','Cone','material',MatName,'thickness',t,'rzbot',[r2 z2],'rztop',[r1 z1],'els',200);

% BC1f condition at base of cone
C1 = AQUINAS_Constraint_Object('rcoord',r2,'zcoord',z2,'dofs',{'u'});
C2 = AQUINAS_Constraint_Object('rcoord',r2,'zcoord',z2,'dofs',{'v'});
C3 = AQUINAS_Constraint_Object('rcoord',r2,'zcoord',z2,'dofs',{'w'});

% BC2f condition at top of cone
C4 = AQUINAS_Constraint_Object('rcoord',r1,'zcoord',z1,'dofs',{'u'});
C5 = AQUINAS_Constraint_Object('rcoord',r1,'zcoord',z1,'dofs',{'v'});

F1 = AQUINAS_Nodal_Force_Object('type','Fw','rcoord',r1,'zcoord',z1,'magnitude',Ncone);

SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',6,'constrLagrMultF',0.001);

A1 = AQUINAS_Analysis_Object('type','LBA','circumferentialModes',0:Nmax,'noEigenvalues',Emax);

O1 = AQUINAS_Output_Mode_Object('mode','VarArgOut');

LBA_out = AQUINAS_Protocol(M1,S1,C1,C2,C3,C4,C5,F1,SOL1,A1,O1);


%%%%%%%%%%%%%%%%%%%
% Post-processing %
%%%%%%%%%%%%%%%%%%%
% Interrogation of each computed LBA eigenmode, not just the most critical
% one, to identify the no. of meridional half-waves (by fitting a Fourier
% half-sine series) and thus identify its location on the 'Koiter valley'
for N = 0:Nmax % For each circumferential wave number
    for EIG = 1:Emax % For each requested eigenvalue
        amps = zeros(Mmax,1);

        x = 2*pi/L*LBA_out.Shell_Geom.s{1}'; % Scaling geometry to be on [0,2*pi]
        y = cos(alpha)*LBA_out.CircWave(N).EigenMode{EIG}.u{1} + sin(alpha)*LBA_out.CircWave(N).EigenMode{EIG}.w{1};

        for M = 1:Mmax % For each meridional half-wave number
            sum_amp = 0.0;

            % Using the trapezoid rule to fit a Fourier half-sine series
            % for a desired harmonic - adapted from Eq. 5 in Ref. [2]
            for P = 1:length(x)-1
               sum_amp = sum_amp + 0.5*(y(P)*sin(M*x(P)*0.5) + y(P+1)*sin(M*x(P+1))*0.5);
            end
            amps(M) = sum_amp*(x(2)-x(1))/pi;
        end
        [~, M] = max(abs(amps));
        if M <= Mmax
            LBA_LPF(M,N+1) = LBA_out.CircWave(N).EigenValue{EIG};
        end
    end
end

bm = nan((Nmax+1)*Mmax,1); bc = nan((Nmax+1)*Mmax,1);
rho_bar = (r1/cos(alpha) + r2/cos(alpha))/2;
lambda_bar = pi*sqrt(rho_bar*t)/sqrt(sqrt(3*(1-nu^2)));
for N = 0:Nmax
    for M = 1:Mmax
        if ~isnan(LBA_LPF(M,N+1))
            if LBA_LPF(M,N+1)<1+lim
                bm(N*Mmax+M) = lambda_bar/(L/M);
                bc(N*Mmax+M) = lambda_bar/(pi*rho_bar/N);
            end
        end
    end
end

abar = 1/sqrt(2);
bbar = cos(alpha)/sqrt(2);
bm_bar = 0:2*abar/1000:2*abar;
bc_bar = bbar*sqrt(1-((bm_bar-abar).^2)/(abar^2));
A_upper = (4*(1+lim)+sqrt(16*((1+lim)^2)-16))/8;
abar_upper = sqrt(A_upper);
bbar_upper = cos(alpha)*sqrt(A_upper);
bm_bar_upper = 0:2*abar_upper/1000:2*abar_upper;
bc_bar_upper = bbar_upper*sqrt(1-((bm_bar_upper-abar_upper).^2)/(abar_upper^2));

A_lower = (4*(1+lim)-sqrt(16*((1+lim)^2)-16))/8;
abar_lower = sqrt(A_lower);
bbar_lower = cos(alpha)*sqrt(A_lower);
bm_bar_lower = 0:2*abar_lower/1000:2*abar_lower;
bc_bar_lower = bbar_lower*sqrt(1-((bm_bar_lower-abar_lower).^2)/(abar_lower^2));


%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%
figure('Name','Koiter ellipse','NumberTitle','off','units','normalized','outerposition',[0 0 1 1],'WindowState','Maximized','color','w');
lgnd = cell(3,1);

h1 = plot(bm_bar, bc_bar,...
    'Color','k','LineStyle','-','Marker','none'); grid on; hold on;
lgnd{1} = 'Koiter ellipse ($N=N_{cone}$)';

h2 = plot(bm_bar_upper, bc_bar_upper,...
    'Color',[0.5 0.5 0.5],'LineStyle','--','Marker','none');
lgnd{2} = strcat('Koiter ellipse ($N=',num2str(1+lim),'N_{cone}$)');

plot(bm_bar_lower, bc_bar_lower,...
        'Color',[0.5 0.5 0.5],'LineStyle','--','Marker','none');

h3 = plot(bm, bc,...
    'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',MS); grid on;
lgnd{3} = 'AQUINAS FE solution';

axis equal; ylim([0 1]);
legend([h1 h2 h3],lgnd,'interpreter','latex','fontsize',FSL,'location','northeast');
xlabel('Dimensionless wavenumber $\bar\beta_{m}$','interpreter','latex','fontsize',FS);
ylabel('Dimensionless wavenumber $\bar\beta_{c}$','interpreter','latex','fontsize',FS);
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