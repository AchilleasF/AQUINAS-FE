%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example LBA 02: LBA of an unpressurised thin cylindrical shell under
% uniform axial compression and direct comparison with 'Koiter circle'.
% A BC1f / S1 condition is assumed at the base, and a BC2f / S3 condition at the top.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% [1] Sadowski A.J., Pototschnig L. & Constantinou P. (2018) "The 'panel
% analysis' technique in the computational study of axisymmetric
% thin-walled shell systems" Finite Elements in Analysis and Design, 152,
% 55-68. https://doi.org/10.1016/j.finel.2018.07.004
% [2] Lin X. & Teng J.G. (2003) "Iterative Fourier decomposition of
% imperfection measurements at non-uniformtly distributed sampling points"
% Thin-Walled Structures 41, 901-924.
% https://doi.org/10.1016/S0263-8231(03)00041-7
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
h = 1000.0; % [mm] - Cylinder length / height

% Material
MatName = 'Steel'; % Material name
E = 2.0e5; % [N/mm2] - Young's modulus of steel
nu = 0.3; % [-] - Poisson ratio of steel

% Loading
Ncl = E*t*t/(r*sqrt(3.0*(1.0-nu*nu))); % [N/mm] - Distributed axisymmetric axial load along top edge (acting vertically downwards)

% Plot controls
FS = 24; % [-] - Font Size (pts)
FSL = 20; % [-] - Font Size for Legend (pts)
MS = 4; % [-] - Markser size (pts)
LW = 2; % [-] - Line Width (pts)


%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%
% Ref. [1] contains many references to classical texts such as those of
% Koiter which are the origin of the Koiter circle theory used here.
nmax = 0.5*(12.0*(1.0-nu*nu))^0.25 * sqrt(r/t); % Max. circumferential wave number on the Koiter circle
mmax = 2.0*nmax*h/(pi*r); % Max. meridional half-wave number on the Koiter circle
Nmax = ceil(1.5 * nmax); Mmax = ceil(1.5 * mmax); Emax = 20;
[Ns,Ms] = meshgrid(0:Nmax, 1:Mmax); LBA_LPF = nan(size(Ns));

M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)));

S1 = AQUINAS_Segment_Object('type','Cone','material',MatName,'thickness',t,'rzbot',[r 0.0],'rztop',[r h],'els',200);

% BC1f condition at base of cylinder
C1 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'u'});
C2 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'v'});
C3 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'w'});

% BC2f condition at top of cylinder
C4 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',h,'dofs',{'u'});
C5 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',h,'dofs',{'v'});

F1 = AQUINAS_Nodal_Force_Object('type','Fw','rcoord',r,'zcoord',h,'magnitude',-Ncl);

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
        x = 2*pi/h*LBA_out.Shell_Geom.z{1}'; % Scaling geometry to be on [0,2*pi]
        y = LBA_out.CircWave(N).EigenMode{EIG}.u{1};

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

        if M <= size(LBA_LPF,1) && (isnan(LBA_LPF(M,N+1)) || LBA_out.CircWave(N).EigenValue{EIG} < LBA_LPF(M,N+1))
            LBA_LPF(M,N+1) = LBA_out.CircWave(N).EigenValue{EIG};
        end
    end
end

% Analytical 3D N/Ncl valley - see Eq.1 in Ref. [1]
[m,n] = meshgrid(0:Mmax,0:Nmax); k = m*pi*r/h;
D = E*t*t*t/(12*(1-nu*nu));
mbar = m*pi/h;
eta2 = (((mbar*r).^2 + n.^2)./(mbar*r^2)).^2;
NoNcl = (D*eta2 + E*(t/r^2)./eta2)/Ncl;


%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%
%% 3D Koiter surface plot
figure('Name','Koiter surface','NumberTitle','off','units','normalized','outerposition',[0 0 1 1],'WindowState','Maximized','color','w','PaperType','A3','PaperOrientation','landscape');
lgnd = cell(3,1);

h1 = surf(Ms,Ns,LBA_LPF,'FaceAlpha',1.0,'FaceColor',[1.0 0.0 0.0]); grid on; hold on;
lgnd{1} = 'AQUINAS FE LBA eigenmodes';

h2 = surf(m,n,NoNcl,'FaceAlpha',0.2,'FaceColor',[0.2 0.2 0.2]);
lgnd{2} = 'Analytical $\frac{N}{N_{cl}}$ 3D surface';

tt = 0:pi/100:pi;
h3 = plot3((nmax*cos(tt)+nmax)*h/(pi*r),nmax*sin(tt),1.01*ones(size(tt)),'k-','LineWidth',LW);
% h3 = scatter3((nmax*cos(tt)+nmax)*h/(pi*r),nmax*sin(tt),ones(size(tt)),10,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k');
lgnd{3} = 'Analytical critical 2D Koiter circle';

legend([h1 h2 h3],lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('No. of meridional half-waves $m$','interpreter','latex','fontsize',FS);
ylabel('No. of circumferential full waves $n$','interpreter','latex','fontsize',FS);
zlabel('Dimensionless critical buckling load $\frac{N}{N_{cl}}$','interpreter','latex','fontsize',FS);
set(gca,'ticklabelinterpreter','latex','fontsize',FS,'tickdir','out');
axis([0 Mmax 0 Nmax 0.9 3]); axis normal;

% Any AQUINAS FE points that do *not* fall on the correct analytical
% surface are due to an insufficient mesh resolution to unambiguously
% capture higher-order computed LBA eigenmodes.


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