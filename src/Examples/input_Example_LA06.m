%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example LA 06: LA of a spherical cap under uniform vertical pressure. The
% meridional displacement of the spherical cap's edge is considered restrained.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% [1] - W. Flügge, Stresses in Shells, 2nd edition, Springer-Verlag, Berlin Heidelberg, Germany, 1973
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
a = 500; % [mm] - Sphere radius
t = 1; % [mm] - Sphere thickness
alpha = 10*pi/180; % [rad] Half angle of spherical shell section (see Fig 15. of Chapter 6 of [1])

% Material
MatName = 'Steel'; % Material name
E = 2.0e5; % [N/mm2] - Young's modulus of steel
nu = 0.3; % [-] - Poisson ratio of steel

% Shell segment FE discretization
numelem = 50;

p = 1; % [N/mm2] - Distributed vertical pressure along middle surface (positive acting upwards)

% Plot controls
FS = 20; % [-] - Font Size (pts)
FSL = 18; % [-] - Font Size for Legends (pts)
MS = 3; % [-] - Markser size (pts)
LW = 3; % [-] - Line Width (pts)

%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%

M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)));

S1 = AQUINAS_Segment_Object('type','Ellipse','material','Steel','formulation','thin','thickness',t,'rzbot',[sin(alpha)*a  cos(alpha)*a],'rztop',[0 a],'geom',[0 0 a a],'els',numelem);

C1 = AQUINAS_Constraint_Object('rcoord',sin(alpha)*a,'zcoord',cos(alpha)*a,'dofs',{'w'});

P1 = AQUINAS_Distributed_Pressure_Object('segments',{S1},'type','pn','functionHandle',@(R) -p*(cos(asin(R/a))^2),'withRespectTo','r');
P2 = AQUINAS_Distributed_Pressure_Object('segments',{S1},'type','pt','functionHandle',@(R) p*(cos(asin(R/a))*sin(asin(R/a))),'withRespectTo','r');

SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',6);

A1 = AQUINAS_Analysis_Object('type','LA');

O1 = AQUINAS_Output_Mode_Object('mode','VarArgOut');

LA_out = AQUINAS_Protocol(M1,S1,C1,P1,P2,SOL1,A1,O1);


%%%%%%%%%%%%%%%%%%%%%%%
% Analytical solution %
%%%%%%%%%%%%%%%%%%%%%%%
% The following bending theory solution can be found in chapter 6.2.1 of [1],
% while the example itself can be found in pages 346-349 of [1].

% Please note a difference in notation below. In AQUINAS, 'u' and 'w' are the
% global displacements parallel to the r and z axes respectively. According to
% Flugge's notation in [1], however, 'w' is the displacement normal to the shell's
% wall, while 'v' is the meridional displacement. Consequently, in the below the
% meaning of 'w' and 'v' is the OPPOSITE of that used in AQUINAS for 'w' and 'u'.

numstan = 100*ceil(alpha*180/pi) + 1; % Number of Analytical Stations, necessary for desired accuracy of displacement field

% Membrane theory solution
Nphi_m = @(phi) -0.5*p*a*ones(1,length(phi)); % Meridional membrane stress resultant
Ntheta_m = @(phi) 0.5*p*a*(1-2*(cos(phi).^2)); % Circumferential membrane stress resultant
EpsPhi_m = @(phi) (Nphi_m(phi)-nu*Ntheta_m(phi))/(E*t); % Meridional midsurface strain for membrane stress state only
EpsTheta_m = @(phi) (Ntheta_m(phi)-nu*Nphi_m(phi))/(E*t); % Circumferential midsurface strain for membrane stress state only
SigmaPhi_m = @(phi) Nphi_m(phi)/t; % Meridional membrane stress
SigmaTheta_m = @(phi) Ntheta_m(phi)/t; % Circumferential membrane stress

% Bending Theory solution
kappa4 = 3*(1-(nu^2))*((a^2)/(t^2)) - (nu^2)/4;  kappa2 = sqrt(kappa4);  kappa = sqrt(kappa2); % kappa value, eq. 6.20 of [1]
D = E*t/(1-nu^2); % Membrane stiffness
K = (t^2)*D/12; % Bending stiffness

% Auxiliary definitions of Kelvin functions of the 1st and 2nd kind, together with their derivatives
e1pi = exp(pi*complex(0,1)/4);
e3pi = exp(3*pi*complex(0,1)/4);
ber = @(n,y) real(besselj(n,e3pi*y));
bei = @(n,y) imag(besselj(n,e3pi*y));
dber = @(n,y) (real(besselj(n+1,y*e3pi)) + imag(besselj(n+1,y*e3pi)))/sqrt(2) + ber(n,y)*n/y;
dbei = @(n,y) (imag(besselj(n+1,y*e3pi)) - real(besselj(n+1,y*e3pi)))/sqrt(2) + bei(n,y)*n/y;
ker = @(n,y) real(exp(-n*pi*complex(0,1)/2)*besselk(n,e1pi*y));
kei = @(n,y) imag(exp(-n*pi*complex(0,1)/2)*besselk(n,e1pi*y));
dker = @(n,y) (real(exp(-(n+1)*pi*complex(0,1)/2)*besselk(n+1,e1pi*y)) + imag(exp(-(n+1)*pi*complex(0,1)/2)*besselk(n+1,e1pi*y)))/sqrt(2) + ker(n,y)*n/y;
dkei = @(n,y) (imag(exp(-(n+1)*pi*complex(0,1)/2)*besselk(n+1,e1pi*y)) - real(exp(-(n+1)*pi*complex(0,1)/2)*besselk(n+1,e1pi*y)))/sqrt(2) + kei(n,y)*n/y;

% Factors for the top edge boundary condition, needed for the evaluation of A1 and A2 constants (only those two needed due to the apex of the spherical cap being closed)
x = @(phi) kappa*phi*sqrt(2);

% Functions multiplied with the A1 and A2 factors respectively to find the rotation of the meridian, according to eq. 6.42b of [1]
xA1f = @(phi) (2*kappa2*dbei(0,x(phi)) - nu*dber(0,x(phi)))/(D*(1-nu^2));
xA2f = @(phi) -(2*kappa2*dber(0,x(phi)) + nu*dbei(0,x(phi)))/(D*(1-nu^2));

% Functions multiplied with the A1 and A2 factors respectively to find the meridional transverse shear force, according to eq. 6.42a of [1]
QphiA1f = @(phi) dber(0,x(phi));
QphiA2f = @(phi) dbei(0,x(phi));

% Functions multiplied with the A1 and A2 factors respectively to find the meridional membrane stress resultant, according to eq. 6.42c of [1]
NphiA1f = @(phi) -QphiA1f(phi)./phi;
NphiA2f = @(phi) -QphiA2f(phi)./phi;

% Functions multiplied with the A1 and A2 factors respectively to find the circumferential membrane stress resultant, according to eq. 6.42d of [1]
NthetaA1f = @(phi) kappa*sqrt(2)*(bei(0,x(phi)) + dber(0,x(phi))./x(phi));
NthetaA2f = @(phi) -kappa*sqrt(2)*(ber(0,x(phi)) - dbei(0,x(phi))./x(phi));

% Functions multiplied with the A1 and A2 factors respectively to find the meridional bending stress resultant, according to eq. 6.42e of [1]
MphiA1f = @(phi) (K*kappa*sqrt(2)/D/a/(1-nu^2))*(2*kappa2*(ber(0,x(phi)) - (1-nu)*dbei(0,x(phi))./x(phi)) + nu*(bei(0,x(phi)) + (1-nu)*dber(0,x(phi))./x(phi)));
MphiA2f = @(phi) (K*kappa*sqrt(2)/D/a/(1-nu^2))*(2*kappa2*(bei(0,x(phi)) + (1-nu)*dber(0,x(phi))./x(phi)) - nu*(ber(0,x(phi)) - (1-nu)*dbei(0,x(phi))./x(phi)));

% Functions multiplied with the A1 and A2 factors respectively to find the circumferential bending stress resultant, according to eq. 6.42f of [1]
MthetaA1f = @(phi) (K*kappa*sqrt(2)/D/a/(1-nu^2))*(2*kappa2*(nu*ber(0,x(phi)) + (1-nu)*dbei(0,x(phi))./x(phi)) + nu*(nu*bei(0,x(phi)) - (1-nu)*dber(0,x(phi))./x(phi)));
MthetaA2f = @(phi) (K*kappa*sqrt(2)/D/a/(1-nu^2))*(2*kappa2*(nu*bei(0,x(phi)) - (1-nu)*dber(0,x(phi))./x(phi)) - nu*(nu*ber(0,x(phi)) + (1-nu)*dbei(0,x(phi))./x(phi)));

% Forming and solving matrix system to find A1 and A2
BC1A1 = QphiA1f(alpha);
BC1A2 = QphiA2f(alpha);
BC2A1 = MphiA1f(alpha);
BC2A2 = MphiA2f(alpha);

BCfactors = [BC1A1   BC1A2
             BC2A1   BC2A2];
V = [-0.5*p*a*cos(alpha)*sin(alpha) 0]';
A = BCfactors\V; A1 = A(1); A2 = A(2);

x_b = @(phi) A1*xA1f(phi) + A2*xA2f(phi); % Meridional rotation
Nphi_b = @(phi) Nphi_m(phi) + A1*NphiA1f(phi) + A2*NphiA2f(phi); % Bending theory meridional membrane stress resultant
Ntheta_b = @(phi) Ntheta_m(phi) + A1*NthetaA1f(phi) + A2*NthetaA2f(phi); % Bending theory circumferential membrane stress resultant
Mphi_b = @(phi) A1*MphiA1f(phi) + A2*MphiA2f(phi); % Bending theory bending membrane stress resultant
Mtheta_b = @(phi) A1*MthetaA1f(phi) + A2*MthetaA2f(phi); % Bending theory circumferential bending stress resultant
% The following stresses-strains are obtained from general shell theory relations
SigmaPhi_bi = @(phi) Nphi_b(phi)/t + 6.0*Mphi_b(phi)/(t*t); % Complete meridional membrane stress - inner surface
SigmaPhi_bm = @(phi) Nphi_b(phi)/t; % Complete meridional membrane stress - midfsurface
SigmaPhi_bo = @(phi) Nphi_b(phi)/t - 6.0*Mphi_b(phi)/(t*t); % Complete meridional membrane stress - outer surface
SigmaTheta_bi = @(phi) Ntheta_b(phi)/t + 6.0*Mtheta_b(phi)/(t*t); % Complete circumferential membrane stress - inner surface
SigmaTheta_bm = @(phi) Ntheta_b(phi)/t; % Complete circumferential membrane stress - midsurface
SigmaTheta_bo = @(phi) Ntheta_b(phi)/t - 6.0*Mtheta_b(phi)/(t*t); % Complete circumferential membrane stress - outer surface
KappaPhi = @(phi) (Mphi_b(phi)-nu*Mtheta_b(phi))/(K*(1-nu^2)); % Meridional curvature
KappaTheta = @(phi) (Mtheta_b(phi)-nu*Mphi_b(phi))/(K*(1-nu^2)); % Circumferential curvature
EpsPhi_bi = @(phi) (Nphi_b(phi)-nu*Ntheta_b(phi))/(D*(1-nu^2)) + 0.5*t*KappaPhi(phi); % Complete meridional (axial) strain - inner surface
EpsPhi_bm = @(phi) (Nphi_b(phi)-nu*Ntheta_b(phi))/(D*(1-nu^2)); % Complete meridional (axial) strain - midsurface
EpsPhi_bo = @(phi) ((Nphi_b(phi)-nu*Ntheta_b(phi))/(D*(1-nu^2))) - 0.5*t*KappaPhi(phi); % Complete meridional (axial) strain - outer surface
EpsTheta_bi = @(phi) ((Ntheta_b(phi)-nu*Nphi_b(phi))/(D*(1-nu^2))) + 0.5*t*KappaTheta(phi); % Complete circumferential strain - inner surface
EpsTheta_bm = @(phi) (Ntheta_b(phi)-nu*Nphi_b(phi))/(D*(1-nu^2)); % Complete circumferential strain - midsurface
EpsTheta_bo = @(phi) ((Ntheta_b(phi)-nu*Nphi_b(phi))/(D*(1-nu^2))) - 0.5*t*KappaTheta(phi); % Complete circumferential strain - outer surface

% Numerical scheme for the calculation of the midsurface displacements,
% that capitilizes on the boundary conditions at phi=alpha (the axial/vertical displacement is equal to zero there)
delta = @(phi) EpsTheta_bm(phi)*(a*sin(phi)); % radial displacement
w_next = delta(alpha)*sin(alpha); % normal displacement of next integration station
v_next = delta(alpha)*cos(alpha); % meridional displacement of next integration station
dvds = EpsPhi_bm(alpha)-w_next/a; % first derivative of meridional displacement with respect to meiridional coordinate of next integration station
vg = zeros(numstan,1);  wg = vg;
vg(1) = cos(alpha)*v_next + sin(alpha)*w_next; wg(1) = cos(alpha)*w_next - sin(alpha)*v_next;
for i=2:numstan
    iphi = alpha - (i-1)*alpha/(numstan-1); % phi angle of current integrations station
    vbp = v_next - dvds*a*alpha/(numstan-1); % numerical integration of meridional displacements, based on the finite difference method, assuming a finite length of ds=a*alpha/(numstan-1)
    wbp = a*EpsTheta_bm(iphi) - vbp*cos(iphi)/sin(iphi); % evaluation of normal displacement at current numerical station, according to eqs. 6.4 of [1], simplified for spherical shell case (r1=r2=a)
    dvds = EpsPhi_bm(iphi) - wbp/a; % evaluation of meridional displacement derivative at current numerical station, according to eqs. 6.4 of [1], simplified for spherical shell case (r1=r2=a)
    v_next = vbp;
    vg(i) = cos(iphi)*vbp + sin(iphi)*wbp; % transformation for the evaluation of global radial displacement v (local to global transformation), for comparison with AQUINAS
    wg(i) = cos(iphi)*wbp - sin(iphi)*vbp; % transformation for the evaluation of global axial displacement w (local to global transformation), for comparison with AQUINAS
end
wg = flipud(wg); vg = flipud(vg);

% %%%%%%%%%%%%
% % Plotting %
% %%%%%%%%%%%%
%% Radial displacement plot
figure('Name','Radial displacement plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(2,1);

% % Shell bending theory solution
plot(vg, 0:alpha/(numstan-1):alpha,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','none'); grid on; hold on;
lgnd{1} = 'Bending theory displacement';

% AQUINAS solution
plot(LA_out.DOFs.u{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',MS);
plot(LA_out.DOFs.u{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','w','MarkerSize',0.8*MS);
lgnd{2} = 'AQUINAS FE solution';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Radial displacement of the wall [$mm$]','interpreter','latex','fontsize',FS);
ylabel('Phi angle from the apex [$rad$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Axial displacement plot
figure('Name','Axial displacement plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(2,1);

% Shell bending theory solution
plot(wg, 0:alpha/(numstan-1):alpha,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','none'); grid on; hold on;
lgnd{1} = 'Bending theory displacement';

% AQUINAS solution
plot(LA_out.DOFs.w{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',MS);
plot(LA_out.DOFs.w{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','w','MarkerSize',0.8*MS);
lgnd{2} = 'AQUINAS FE solution';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Axial displacement of the wall [$mm$]','interpreter','latex','fontsize',FS);
ylabel('Phi angle from the apex [$rad$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Meridional strain plot
figure('Name','Meridional strain plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(7,1);

% Shell membrane theory solution
plot(EpsPhi_m(0:alpha/(numstan-1):alpha), 0:alpha/(numstan-1):alpha,...
    'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',LW); grid on; hold on;
lgnd{1} = 'Membrane theory midsurface strain';

% Shell bending theory solution
plot(EpsPhi_bi(0:alpha/(numstan-1):alpha), 0:alpha/(numstan-1):alpha,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','none');
plot(EpsPhi_bm(0:alpha/(numstan-1):alpha), 0:alpha/(numstan-1):alpha,...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','none');
plot(EpsPhi_bo(0:alpha/(numstan-1):alpha), 0:alpha/(numstan-1):alpha,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','none');
lgnd{2} = 'Bending theory strain (inner surface)';
lgnd{3} = 'Bending theory strain (midsurface)';
lgnd{4} = 'Bending theory strain (outer surface)';

% AQUINAS solution
plot(LA_out.Strains.Epsilon.Eps_phi.i{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_phi.m{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_phi.o{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_phi.i{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Strains.Epsilon.Eps_phi.m{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(LA_out.Strains.Epsilon.Eps_phi.o{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{5} = 'AQUINAS FE solution (inner surface)';
lgnd{6} = 'AQUINAS FE solution (midsurface)';
lgnd{7} = 'AQUINAS FE solution (outer surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Meridional strains in the wall [$-$]','interpreter','latex','fontsize',FS);
ylabel('Phi angle from the apex [$rad$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Circumferential strain plot
figure('Name','Circumferential strain plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(7,1);

% Shell membrane theory solution
plot(EpsTheta_m(0:alpha/(numstan-1):alpha), 0:alpha/(numstan-1):alpha,...
    'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',LW); grid on; hold on;
lgnd{1} = 'Membrane theory midsurface strain';

% Shell bending theory solution
plot(EpsTheta_bi(0:alpha/(numstan-1):alpha), 0:alpha/(numstan-1):alpha,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','none');
plot(EpsTheta_bm(0:alpha/(numstan-1):alpha), 0:alpha/(numstan-1):alpha,...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','none');
plot(EpsTheta_bo(0:alpha/(numstan-1):alpha), 0:alpha/(numstan-1):alpha,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','none');
lgnd{2} = 'Bending theory strain (inner surface)';
lgnd{3} = 'Bending theory strain (midsurface)';
lgnd{4} = 'Bending theory strain (outer surface)';

% AQUINAS solution
plot(LA_out.Strains.Epsilon.Eps_theta.i{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_theta.m{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_theta.o{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_theta.i{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Strains.Epsilon.Eps_theta.m{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(LA_out.Strains.Epsilon.Eps_theta.o{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{5} = 'AQUINAS FE solution (inner surface)';
lgnd{6} = 'AQUINAS FE solution (midsurface)';
lgnd{7} = 'AQUINAS FE solution (outer surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Circumferential strains in the wall [$-$]','interpreter','latex','fontsize',FS);
ylabel('Phi angle from the apex [$rad$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Meridional stress plot
figure('Name','Meridional stress plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(7,1);

% Shell membrane theory solution
plot(SigmaPhi_m(0:alpha/(numstan-1):alpha), 0:alpha/(numstan-1):alpha,...
    'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',LW); grid on; hold on;
lgnd{1} = 'Membrane theory midsurface stress';

% Shell bending theory solution
plot(SigmaPhi_bi(0:alpha/(numstan-1):alpha), 0:alpha/(numstan-1):alpha,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','none');
plot(SigmaPhi_bm(0:alpha/(numstan-1):alpha), 0:alpha/(numstan-1):alpha,...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','none');
plot(SigmaPhi_bo(0:alpha/(numstan-1):alpha), 0:alpha/(numstan-1):alpha,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','none');
lgnd{2} = 'Bending theory stress (inner surface)';
lgnd{3} = 'Bending theory stress (midsurface)';
lgnd{4} = 'Bending theory stress (outer surface)';

% AQUINAS solution
plot(LA_out.Stresses.Sig_phi.i{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stresses.Sig_phi.m{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(LA_out.Stresses.Sig_phi.o{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stresses.Sig_phi.i{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Stresses.Sig_phi.m{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(LA_out.Stresses.Sig_phi.o{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{5} = 'AQUINAS FE solution (inner surface)';
lgnd{6} = 'AQUINAS FE solution (midsurface)';
lgnd{7} = 'AQUINAS FE solution (outer surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Meridional stresses in the wall [$N/mm^2$]','interpreter','latex','fontsize',FS);
ylabel('Phi angle from the apex [$rad$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');



%% Circumferential stress plot
figure('Name','Circumferential stress plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(7,1);

% Shell membrane theory solution
plot(SigmaTheta_m(0:alpha/(numstan-1):alpha), 0:alpha/(numstan-1):alpha,...
    'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',LW); grid on; hold on;
lgnd{1} = 'Membrane theory midsurface stress';

% Shell bending theory solution
plot(SigmaTheta_bi(0:alpha/(numstan-1):alpha), 0:alpha/(numstan-1):alpha,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','none');
plot(SigmaTheta_bm(0:alpha/(numstan-1):alpha), 0:alpha/(numstan-1):alpha,...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','none');
plot(SigmaTheta_bo(0:alpha/(numstan-1):alpha), 0:alpha/(numstan-1):alpha,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','none');
lgnd{2} = 'Bending theory stress (inner surface)';
lgnd{3} = 'Bending theory stress (midsurface)';
lgnd{4} = 'Bending theory stress (outer surface)';

% AQUINAS solution
plot(LA_out.Stresses.Sig_theta.i{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stresses.Sig_theta.m{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(LA_out.Stresses.Sig_theta.o{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stresses.Sig_theta.i{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Stresses.Sig_theta.m{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(LA_out.Stresses.Sig_theta.o{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{5} = 'AQUINAS FE solution (inner surface)';
lgnd{6} = 'AQUINAS FE solution (midsurface)';
lgnd{7} = 'AQUINAS FE solution (outer surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Circumferential stresses in the wall [$N/mm^2$]','interpreter','latex','fontsize',FS);
ylabel('Phi angle from the apex [$rad$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Membrane stress resultants
figure('Name','Membrane stress resultant plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(6,1);

% Shell membrane theory solution
plot(Nphi_m(0:alpha/(numstan-1):alpha), 0:alpha/(numstan-1):alpha,...
    'Color',[1.0 0.7 0.7],'LineStyle','--','LineWidth',LW); grid on; hold on;
plot(Ntheta_m(0:alpha/(numstan-1):alpha), 0:alpha/(numstan-1):alpha,...
    'Color',[0.7 0.7 1.0],'LineStyle','--','LineWidth',LW);
lgnd{1} = 'Membrane theory meridional MSR';
lgnd{2} = 'Membrane theory circumferential MSR';

% Shell bending theory solution
plot(Nphi_b(0:alpha/(numstan-1):alpha), 0:alpha/(numstan-1):alpha,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','none');
plot(Ntheta_b(0:alpha/(numstan-1):alpha), 0:alpha/(numstan-1):alpha,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','none');
lgnd{3} = 'Bending theory meridional MSR';
lgnd{4} = 'Bending theory circumferential MSR';

% AQUINAS solution
plot(LA_out.Stress_Reslts.Membrane.Nphi{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stress_Reslts.Membrane.Ntheta{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stress_Reslts.Membrane.Nphi{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Stress_Reslts.Membrane.Ntheta{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{5} = 'AQUINAS FE solution meridional MSR';
lgnd{6} = 'AQUINAS FE solution circumferential MSR';

legend(lgnd(3:end),'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Membrane stress resultants (MSR) in the wall [$N/mm$]','interpreter','latex','fontsize',FS);
ylabel('Phi angle from the apex [$rad$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Bending stress resultants
figure('Name','Bending stress resultant plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(4,1);

% Shell bending theory solution
plot(Mphi_b(0:alpha/(numstan-1):alpha), 0:alpha/(numstan-1):alpha,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
plot(Mtheta_b(0:alpha/(numstan-1):alpha), 0:alpha/(numstan-1):alpha,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{1} = 'Bending theory meridional BSR';
lgnd{2} = 'Bending theory circumferential BSR';

% AQUINAS solution
plot(-LA_out.Stress_Reslts.Bending.Mphi{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(-LA_out.Stress_Reslts.Bending.Mtheta{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(-LA_out.Stress_Reslts.Bending.Mphi{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(-LA_out.Stress_Reslts.Bending.Mtheta{1}, asin(LA_out.Shell_Geom.r{1}/a),...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{3} = 'AQUINAS FE solution meridional BSR';
lgnd{4} = 'AQUINAS FE solution circumferential BSR';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Bending stress resultants (BSR) in the wall [$Nmm/mm$]','interpreter','latex','fontsize',FS);
ylabel('Phi angle from the apex [$rad$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


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