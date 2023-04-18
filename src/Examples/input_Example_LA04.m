%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example LA 04: LA of a thin conical shell under axial/vertical
% uniform pressure. The meridian of the cone is pinned at its base.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% [1] Timoshenko, Stephen, and S. Woinowsky-Krieger - Theory of Plates and Shells,
% 2nd ed, New York, McGraw-Hill, 1959.
% [2] W. Flügge - Stresses in Shells,
% 2nd edition, Berlin Heidelberg, Springer-Verlag, Germany, 1973
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
t = 5; % [mm] - Cone thickness

% Material
MatName = 'Steel'; % Material name
E = 2.0e5; % [N/mm2] - Young's modulus of steel
nu = 0.3; % [-] - Poisson ratio of steel

L = 200; % [mm] Total length of the cone's meridian
a = 30*pi/180; % [rad] Half angle of angle span of the apex of the cone

% Loading
q = 1; % [N/mm2] - Axial uniform distributed pressure (self weight)
pn = -q*sin(a); % [N/mm2] - Distributed normal pressure along inner surface (positive acting radially outwards)
pt =  q*cos(a); % [N/mm2] - Distributed tangential pressure (positive acting tangentially towards the base of the cone)

% Plot controls
FS = 20; % [-] - Font Size (pts)
FSL = 18; % [-] - Font Size for Legends (pts)
MS = 3; % [-] - Markser size (pts)
LW = 3; % [-] - Line Width (pts)


%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%

M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)));

S1 = AQUINAS_Segment_Object('type','Cone','material',MatName,'thickness',t,'RZbot',[sin(a)*L 0],'RZtop',[0 cos(a)*L],'els',100);

C1 = AQUINAS_Constraint_Object('rcoord',sin(a)*L,'zcoord',0,'dofs',{'u'});
C2 = AQUINAS_Constraint_Object('rcoord',sin(a)*L,'zcoord',0,'dofs',{'w'});

P1 = AQUINAS_Distributed_Pressure_Object('segments',{S1},'type','pn','functionHandle',@() pn);
P2 = AQUINAS_Distributed_Pressure_Object('segments',{S1},'type','pt','functionHandle',@() pt);

SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',6);

A1 = AQUINAS_Analysis_Object('type','LA');

O1 = AQUINAS_Output_Mode_Object('mode','VarArgOut');

LA_out = AQUINAS_Protocol(M1,S1,C1,C2,P1,P2,SOL1,A1,O1);

%%%%%%%%%%%%%%%%%%%%%%%
% Analytical solution %
%%%%%%%%%%%%%%%%%%%%%%%
% Please note a difference in notation below. In AQUINAS, u and w are the
% global displacements parallel to the r and z axes respectively. According to
% the notation used by Timoshenko, however, w is the displacement normal to the wall
% and v is the meridional displacement, while 'delta' is used for the radial displacement
% (equivalent to AQUINAS's global u). In addition, 'upsilon' is adopted here for the axial
% displacements, although the same notation is not used in Timoshenko's book (no reference
% on the axial displacements there).

%  The following problem is presented in article 133 of [1], pages 562-566.

D = E*t/(1-nu^2); % Membrane stiffness
K = (t^2)*D/12; % Bending stiffness

% Membrane theory solution
% The following membrane theory solution has been derived with the additional help of relevant equations presented in [2], or from classical shell theory equations,
% again for the problem of article 133 of [1].
% Notice that a is tha apex half-angle here and not the angle at the base of the cone, as presented in Fig. 2.12 of [2]
Ns_m = @(s) -q*s/2/cos(a); % Meridional membrane stress resultant, according to eq. 2.18 of [2]
Nth_m = @(s) -q*s*sin(a)*tan(a); % Meridional membrane stress resultant, according to eq. 2.17 of [2]
es_mm = @(s) (Ns_m(s) - nu*Nth_m(s))/E/t; % Midsurface meridional strain, according to eq. 2.54 of [2]
eth_mm = @(s) (Nth_m(s) - nu*Ns_m(s))/E/t; % Midsurface circumferential strain, according to eq. 2.54 of [2]
delta_m = @(s) sin(a)*s.*eth_mm(s); % Radial displacement, following classical shell theory (or through observation of equation (p) of article 133 of [1])
v_m = @(s) (q/E/t/cos(a)/2)*(nu*(sin(a)^2)- 1/2)*(s.^2) - (q/E/t/cos(a)/2)*(nu*(sin(a)^2)- 1/2)*(L^2); % Tangential displacement of the meridian, obtained through manipulation and integration of the above defined equations for Ns_m and Nth_m, taking into account eqs. (a) and (b) of article 108 of [1] (or 2.56a and 2.56b of [2])
w_m = @(s) -(q*tan(a)/E/t/cos(a))*(sin(a)^2 - nu/2 - 1/4 + nu*(sin(a)^2)/2)*(s.^2) + tan(a)*(q/E/t/cos(a)/2)*(nu*(sin(a)^2)- 1/2)*(L^2); % Normal displacement of the meridian, obtained through manipulation and integration of the above defined equations for Ns_m and Nth_m, taking into account eqs. (a) and (b) of article 108 of [1] (or 2.56a and 2.56b of [2])
dw_m = @(s) -(q*tan(a)/E/t/cos(a))*(sin(a)^2 - nu/2 - 1/4 + nu*(sin(a)^2)/2)*(2*s); % First derivative of the normal displacement of the meridian with respect to the meridional coordinate s, obtained through derivation of w_m presented above
d2w_m = @(s) -(q*tan(a)/E/t/cos(a))*(sin(a)^2 - nu/2 - 1/4 + nu*(sin(a)^2)/2)*2; % Second derivative of the normal displacement of the meridian with respect to the meridional coordinate s, obtained through derivation of dw_m presented above
upsilon_m = @(s) w_m(s)*sin(a) - v_m(s)*cos(a); % Axial displacement of the meridian, obtained through basic geometric transformation of the tangential and normal displacements of the meridian
SigmaS_mm = @(s) Ns_m(s)/t; % Midsurface meridional stress, following classical shell theory definitions
SigmaTh_mm = @(s) Nth_m(s)/t; % Midsurface circumferential stress, following classical shell theory definitions
% Please note that the following variables are only defined here for auxiliary purposes corresponding to the bending theory solution. These 'equivalent' curvatures and moments from membrane theory arise only because the first and second derivatives of the normal displacement do not vanish
ks_m = @(s) -d2w_m(s); % 'equivalent' membrane theory meridional curvature
kth_m = @(s) -dw_m(s)./s; % 'equivalent' membrane theory circumferential curvature
Ms_m = @(s) -K*(ks_m(s) + nu*kth_m(s)); % 'equivalent' membrane theory meridional moment
Mth_m = @(s) -K*(kth_m(s) + nu*ks_m(s)); % 'equivalent' membrane theory circumferential curvature

% Bending theory solution

numstan = 1000 + 1; % Number of Analytical Stations, necessary for desired accuracy of displacement field

% Auxiliary definitions of Kelvin functions of the 1st and 2nd kind, together with their derivatives
e1pi = exp(pi*complex(0,1)/4);
e3pi = exp(3*pi*complex(0,1)/4);
y = @(s) 2*sqrt(sqrt(3*(1-nu^2)))*sqrt(2*tan(a)/t)*sqrt(s);
dy = @(s) sqrt(sqrt(3*(1-nu^2)))*sqrt(2*tan(a)/t).*(1./sqrt(s));
d2y = @(s) -0.5*sqrt(sqrt(3*(1-nu^2)))*sqrt(2*tan(a)/t).*(1./(sqrt(s).^3));
ber = @(n,y) real(besselj(n,e3pi*y));
bei = @(n,y) imag(besselj(n,e3pi*y));
dber = @(n,y) (real(besselj(n+1,y*e3pi)) + imag(besselj(n+1,y*e3pi)))/sqrt(2) + ber(n,y)*n/(y);
dbei = @(n,y) (imag(besselj(n+1,y*e3pi)) - real(besselj(n+1,y*e3pi)))/sqrt(2) + bei(n,y)*n/(y);
ker = @(n,y) real(exp(-n*pi*complex(0,1)/2)*besselk(n,e1pi*y));
kei = @(n,y) imag(exp(-n*pi*complex(0,1)/2)*besselk(n,e1pi*y));
dker = @(n,y) (real(exp(-(n+1)*pi*complex(0,1)/2)*besselk(n+1,e1pi*y)) + imag(exp(-(n+1)*pi*complex(0,1)/2)*besselk(n+1,e1pi*y)))/sqrt(2) + ker(n,y)*n/(y);
dkei = @(n,y) (imag(exp(-(n+1)*pi*complex(0,1)/2)*besselk(n+1,e1pi*y)) - real(exp(-(n+1)*pi*complex(0,1)/2)*besselk(n+1,e1pi*y)))/sqrt(2) + kei(n,y)*n/(y);

lambda = sqrt(sqrt(12*(1-nu^2)*(cot(a)^2)/(t^2))); % lambda, as presented in eq. (b) of article 133 of [1]
ksi = @(s) 2*lambda*sqrt(s); % ksi, as given in the explanation of eq. (f) of article 133 of [1]
psi1 = @(x) ber(0,x); % psi1 of eq. (f) of article 133 of [1], given in terms of Kelvin function, obtained through comparison with values of Table 86
dpsi1 = @(x) dber(0,x); % dpsi1 of eq. (f) of article 133 of [1], given in terms of a Kelvin function, obtained through comparison with values of Table 86
psi2 = @(x) -bei(0,x); % psi2 of eq. (f) of article 133 of [1], given in terms of a Kelvin function, obtained through comparison with values of Table 86
dpsi2 = @(x) -dbei(0,x); % dpsi2 of eq. (f) of article 133 of [1], given in terms of a Kelvin function, obtained through comparison with values of Table 86

% Factors for the boundary conditions, needed for the evaluation of C1 and C2 constants
QsC1f = @(s) psi1(ksi(s)) + 2*dpsi2(ksi(s))./ksi(s); % Factor of C1 for transverse shear Qy, according to eq. (m) of article 133 of [1]
QsC2f = @(s) psi2(ksi(s)) - 2*dpsi1(ksi(s))./ksi(s); % Factor of C2 for transverse shear Qy, according to eq. (m) of article 133 of [1]

MsC1f = @(s) (2./(ksi(s).^2)).*(-ksi(s).*dpsi2(ksi(s)) + 2*(1-nu)*psi2(ksi(s)) - 4*(1-nu)*dpsi1(ksi(s))./ksi(s)); % Factor of C1 for meridional moment My, according to eq. (o) of article 133 of [1]
MsC2f = @(s) (2./(ksi(s).^2)).*(ksi(s).*dpsi1(ksi(s)) - 2*(1-nu)*psi1(ksi(s)) - 4*(1-nu)*dpsi2(ksi(s))./ksi(s)); % Factor of C2 for meridional moment My, according to eq. (o) of article 133 of [1]

deltaC1f = @(s) -(sin(a)*tan(a)/2/E/t)*(ksi(s).*dpsi1(ksi(s)) - 2*psi1(ksi(s)) - 4*dpsi2(ksi(s))./ksi(s)) + (nu*sin(a)*tan(a)/E/t)*(psi1(ksi(s)) + 2*dpsi2(ksi(s))./ksi(s)); % factor of C1 for radial displacement delta, according to eq. (p) of article 133 of [1]
deltaC2f = @(s) -(sin(a)*tan(a)/2/E/t)*(ksi(s).*dpsi2(ksi(s)) - 2*psi2(ksi(s)) + 4*dpsi1(ksi(s))./ksi(s)) + (nu*sin(a)*tan(a)/E/t)*(psi2(ksi(s)) - 2*dpsi1(ksi(s))./ksi(s)); % factor of C2 for radial displacement delta, according to eq. (p) of article 133 of [1]

% Construction and solution of boundary condition equations in matrix form
BC1C1c = MsC1f(L);
BC1C2c = MsC2f(L);
BC2C1c = deltaC1f(L);
BC2C2c = deltaC2f(L);

BCfactors = [BC1C1c  BC1C2c
             BC2C1c  BC2C2c];

V = [-Ms_m(L) -delta_m(L)]';
C = BCfactors\V;  C1 = C(1);  C2 = C(2);


delta_b = @(s) C1*deltaC1f(s) + C2*deltaC2f(s); % Radial displacement, according to eq. (p) of article 133 of [1], using the above defined equations for the factors of C1 and C2
Qs_b = @(s) (C1*QsC1f(s) + C2*QsC2f(s))./s;  % Transverse shear, according to eq. (m) of article 133 of [1], using the above defined equations for the factors of C1 and C2
Ns_b = @(s) -tan(a)*Qs_b(s); % Meridional membrane stress resultant, according to equation (g) of article 133 of [1]
Nth_b = @(s) E*t*delta_b(s)./s/sin(a) + nu*Ns_b(s); % Circumferential membrane stress resultant, according to equation (L) of article 133 of [1] (or classic shell theory equations)
es_bm = @(s) (Ns_b(s)-nu*Nth_b(s))/E/t; % Meridional membrane midsurface strain, according to classic shell theory equations
eth_bm = @(s) (Nth_b(s)-nu*Ns_b(s))/E/t; % Circumferential membrane midsurface strain, according to classic shell theory equations
ks_b = @(s) -(2./(ksi(s).^2)/K).*(C1*(-ksi(s).*dpsi2(ksi(s)) + 2*psi2(ksi(s)) - 4*dpsi1(ksi(s))./ksi(s)) + C2*(ksi(s).*dpsi1(ksi(s)) - 2*psi1(ksi(s)) - 4*dpsi2(ksi(s))./ksi(s))); % Meridional curvature, obtained through observation of eq. (o) of article 133 of [1]
kth_b = @(s) -(2./(ksi(s).^2)/K).*(C1*(-2*psi2(ksi(s)) + 4*dpsi1(ksi(s))./ksi(s)) + C2*(2*psi1(ksi(s)) + 4*dpsi2(ksi(s))./ksi(s))); % Circumferential curvature, obtained through observation of eq. (o) of article 133 of [1]
Ms_b = @(s) -K*(ks_b(s) + nu*kth_b(s)); % Meridional moment, obtained through classic shell theory equations (or through observation of eq.(i) of article 133 of [1])
Mth_b = @(s) -K*(kth_b(s) + nu*ks_b(s)); % circumferential moment, obtained through classic shell theory equations (or through observation of eq.(i) of article 133 of [1])

% Numerical scheme for the calculation of the midsurface displacements,
% that capitilizes on the boundary conditions at s=L (upsilon=0 => upsilon_b=-upsilon_m  for s=L)
v_b = zeros(numstan,1);  w_b = v_b;
v_next = sin(a)*delta_b(L) + cos(a)*upsilon_m(L); v_b(1) = v_next;
dvds = @(s) es_bm(s); % dvds equal to es_bm, according to eq. (a) of article 108 of [1] (for r1=Inf)
w_b(1) = (delta_b(L) - v_b(1)*sin(a))/cos(a); % Normal displacement of the meridian at the base of the cone, obtained through simple geometric transformation
for i=2:numstan
    is = L - (i-1)*L/(numstan-1);
    is_next = L - (i-2)*L/(numstan-1);
    v_b(i) = v_next - dvds(is_next)*L/(numstan-1); % evaluation of tangential displacement of the meridian by numerical addition of dvds times a finite ds length
    w_b(i) = (delta_b(is) - v_b(i)*sin(a))/cos(a); % Normal displacement of the meridian, obtained through simple geometric transformation
    v_next = v_b(i);
end
upsilon_b = flipud(w_b)*sin(a) - flipud(v_b)*cos(a);

% Superposition of membrane and bending theory solutions
Ns = @(s) Ns_m(s) + Ns_b(s); % Complete meridional membrane stress resultant
Nth = @(s) Nth_m(s) + Nth_b(s); % Complete circumferential membrane stress resultant
Ms = @(s) Ms_m(s) + Ms_b(s); % Complete meridional bending stress resultant
Mth = @(s) Mth_m(s) + Mth_b(s); % Complete circumferential bending stress resultant
ks = @(s) ks_b(s) + ks_m(s); % Complete meridional curvature
kth = @(s) kth_b(s) + kth_m(s); % Complete circumferential curvature
es_o = @(s) (Ns(s)-nu*Nth(s))/E/t + 0.5*t*ks(s);  % Complete meridional strain - outter surface
es_m = @(s) (Ns(s)-nu*Nth(s))/E/t; % Complete meridional strain - middle surface
es_i = @(s) (Ns(s)-nu*Nth(s))/E/t - 0.5*t*ks(s);  % Complete meridional strain - inner surface
eth_o = @(s) (Nth(s)-nu*Ns(s))/E/t + 0.5*t*kth(s); % Complete circumferential strain - outter surface
eth_m = @(s) (Nth(s)-nu*Ns(s))/E/t; % Complete circumferential strain - middle surface
eth_i = @(s) (Nth(s)-nu*Ns(s))/E/t - 0.5*t*kth(s); % Complete circumferential strain - inner surface
SigmaS_i = @(s) Ns(s)/t + 6.0*Ms(s)/(t*t); % Complete meridional membrane stress - inner surface
SigmaS_m = @(s) Ns(s)/t; % Complete meridional membrane stress - midfsurface
SigmaS_o = @(s) Ns(s)/t - 6.0*Ms(s)/(t*t); % Complete meridional membrane stress - outer surface
SigmaTh_i = @(s) Nth(s)/t + 6.0*Mth(s)/(t*t); % Complete circumferential membrane stress - inner surface
SigmaTh_m = @(s) Nth(s)/t; % Complete circumferential membrane stress - midsurface
SigmaTh_o = @(s) Nth(s)/t - 6.0*Mth(s)/(t*t); % Complete circumferential membrane stress - outer surface

delta = @(s) delta_b(s) + delta_m(s);
upsilon = upsilon_b + upsilon_m(0:L/(numstan-1):L)';

% %%%%%%%%%%%%
% % Plotting %
% %%%%%%%%%%%%
%% Radial displacement plot
figure('Name','Radial displacement plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(3,1);

% %% Shell membrane theory solution
plot(delta_m(0:L), 0:L,...
    'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',LW); grid on; hold on;
lgnd{1} = 'Membrane theory displacement';

% Shell bending theory solution
plot(delta(0:L), 0:L,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{2} = 'Bending theory displacement';

% AQUINAS solution
plot( LA_out.DOFs.u{1}, LA_out.Shell_Geom.s{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',MS);
plot( LA_out.DOFs.u{1}, LA_out.Shell_Geom.s{1},...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','w','MarkerSize',0.8*MS);
lgnd{3} = 'AQUINAS FE solution';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Radial displacement of the wall [$mm$]','interpreter','latex','fontsize',FS);
ylabel('Arc-length distance from apex of segment [$mm$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Axial displacement plot
figure('Name','Axial displacement plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(3,1);

% Shell membrane theory solution
plot(upsilon_m(0:L/(numstan-1):L), 0:L/(numstan-1):L,...
    'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',LW); grid on; hold on;
lgnd{1} = 'Membrane theory displacement';

% Shell bending theory solution
plot(upsilon, 0:L/(numstan-1):L,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{2} = 'Bending theory displacement';

% AQUINAS solution
plot( LA_out.DOFs.w{1}, LA_out.Shell_Geom.s{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',MS);
plot( LA_out.DOFs.w{1}, LA_out.Shell_Geom.s{1},...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','w','MarkerSize',0.8*MS);
lgnd{3} = 'AQUINAS FE solution';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Axial displacement of the wall [$mm$]','interpreter','latex','fontsize',FS);
ylabel('Arc-length distance from apex of segment [$mm$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Meridional strain plot
figure('Name','Meridional strain plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(7,1);

% Shell membrane theory solution
plot(es_mm(0:L), 0:L,...
    'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',LW); grid on; hold on;
lgnd{1} = 'Membrane theory midsurface strain';

% Shell bending theory solution
plot(es_i(0:L), 0:L,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(es_m(0:L), 0:L,...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(es_o(0:L), 0:L,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{2} = 'Bending theory strain (inner surface)';
lgnd{3} = 'Bending theory strain (midsurface)';
lgnd{4} = 'Bending theory strain (outer surface)';

% AQUINAS solution
plot(LA_out.Strains.Epsilon.Eps_phi.i{1}, LA_out.Shell_Geom.s{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_phi.m{1}, LA_out.Shell_Geom.s{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_phi.o{1}, LA_out.Shell_Geom.s{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_phi.i{1}, LA_out.Shell_Geom.s{1},...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Strains.Epsilon.Eps_phi.m{1}, LA_out.Shell_Geom.s{1},...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(LA_out.Strains.Epsilon.Eps_phi.o{1}, LA_out.Shell_Geom.s{1},...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{5} = 'AQUINAS FE solution (inner surface)';
lgnd{6} = 'AQUINAS FE solution (midsurface)';
lgnd{7} = 'AQUINAS FE solution (outer surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Meridional strains in the wall [$-$]','interpreter','latex','fontsize',FS);
ylabel('Arc-length distance from apex of segment [$mm$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Circumferential strain plot
figure('Name','Circumferential strain plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(7,1);

% Shell membrane theory solution
plot(eth_mm(0:L), 0:L,...
    'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',LW); grid on; hold on;
lgnd{1} = 'Membrane theory midsurface strain';

% Shell bending theory solution
plot(eth_i(0:L), 0:L,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(eth_m(0:L), 0:L,...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(eth_o(0:L), 0:L,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{2} = 'Bending theory strain (inner surface)';
lgnd{3} = 'Bending theory strain (midsurface)';
lgnd{4} = 'Bending theory strain (outer surface)';

% AQUINAS solution
plot(LA_out.Strains.Epsilon.Eps_theta.i{1}, LA_out.Shell_Geom.s{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_theta.m{1}, LA_out.Shell_Geom.s{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_theta.o{1}, LA_out.Shell_Geom.s{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Strains.Epsilon.Eps_theta.i{1}, LA_out.Shell_Geom.s{1},...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Strains.Epsilon.Eps_theta.m{1}, LA_out.Shell_Geom.s{1},...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(LA_out.Strains.Epsilon.Eps_theta.o{1}, LA_out.Shell_Geom.s{1},...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{5} = 'AQUINAS FE solution (inner surface)';
lgnd{6} = 'AQUINAS FE solution (midsurface)';
lgnd{7} = 'AQUINAS FE solution (outer surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Circumferential strains in the wall [$-$]','interpreter','latex','fontsize',FS);
ylabel('Arc-length distance from apex of segment [$mm$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Meridional stress plot
figure('Name','Meridional stress plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(7,1);

% Shell membrane theory solution
plot(SigmaS_mm(0:L), 0:L,...
    'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',LW); grid on; hold on;
lgnd{1} = 'Membrane theory midsurface stress';

% Shell bending theory solution
plot(SigmaS_i(0:L), 0:L,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(SigmaS_m(0:L), 0:L,...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(SigmaS_o(0:L), 0:L,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{2} = 'Bending theory stress (inner surface)';
lgnd{3} = 'Bending theory stress (midsurface)';
lgnd{4} = 'Bending theory stress (outer surface)';

% AQUINAS solution
plot(LA_out.Stresses.Sig_phi.i{1}, LA_out.Shell_Geom.s{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stresses.Sig_phi.m{1}, LA_out.Shell_Geom.s{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(LA_out.Stresses.Sig_phi.o{1}, LA_out.Shell_Geom.s{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stresses.Sig_phi.i{1}, LA_out.Shell_Geom.s{1},...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Stresses.Sig_phi.m{1}, LA_out.Shell_Geom.s{1},...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(LA_out.Stresses.Sig_phi.o{1}, LA_out.Shell_Geom.s{1},...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{5} = 'AQUINAS FE solution (inner surface)';
lgnd{6} = 'AQUINAS FE solution (midsurface)';
lgnd{7} = 'AQUINAS FE solution (outer surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Meridional stresses in the wall [$N/mm^2$]','interpreter','latex','fontsize',FS);
ylabel('Arc-length distance from apex of segment [$mm$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');



%% Circumferential stress plot
figure('Name','Circumferential stress plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(7,1);

% Shell membrane theory solution
plot(SigmaTh_mm(0:L), 0:L,...
    'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',LW); grid on; hold on;
lgnd{1} = 'Membrane theory midsurface stress';

% Shell bending theory solution
plot(SigmaTh_i(0:L), 0:L,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(SigmaTh_m(0:L), 0:L,...
    'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(SigmaTh_o(0:L), 0:L,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{2} = 'Bending theory stress (inner surface)';
lgnd{3} = 'Bending theory stress (midsurface)';
lgnd{4} = 'Bending theory stress (outer surface)';

% AQUINAS solution
plot(LA_out.Stresses.Sig_theta.i{1}, LA_out.Shell_Geom.s{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stresses.Sig_theta.m{1}, LA_out.Shell_Geom.s{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',MS);
plot(LA_out.Stresses.Sig_theta.o{1}, LA_out.Shell_Geom.s{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stresses.Sig_theta.i{1}, LA_out.Shell_Geom.s{1},...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Stresses.Sig_theta.m{1}, LA_out.Shell_Geom.s{1},...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',0.8*MS);
plot(LA_out.Stresses.Sig_theta.o{1}, LA_out.Shell_Geom.s{1},...
    'LineStyle','none','Marker','s','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{5} = 'AQUINAS FE solution (inner surface)';
lgnd{6} = 'AQUINAS FE solution (midsurface)';
lgnd{7} = 'AQUINAS FE solution (outer surface)';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Circumferential stresses in the wall [$N/mm^2$]','interpreter','latex','fontsize',FS);
ylabel('Arc-length distance from apex of segment [$mm$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Membrane stress resultants
figure('Name','Membrane stress resultant plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(6,1);

% Shell membrane theory solution
plot(Ns_m(0:L), 0:L,...
    'Color',[1.0 0.7 0.7],'LineStyle','--','LineWidth',LW); grid on; hold on;
plot(Nth_m(0:L), 0:L,...
    'Color',[0.7 0.7 1.0],'LineStyle','--','LineWidth',LW);
lgnd{1} = 'Membrane theory meridional MSR';
lgnd{2} = 'Membrane theory circumferential MSR';

% Shell bending theory solution
plot(Ns(0:L), 0:L,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
plot(Nth(0:L), 0:L,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{3} = 'Bending theory meridional MSR';
lgnd{4} = 'Bending theory circumferential MSR';

% AQUINAS solution
plot(LA_out.Stress_Reslts.Membrane.Nphi{1}, LA_out.Shell_Geom.s{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stress_Reslts.Membrane.Ntheta{1}, LA_out.Shell_Geom.s{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stress_Reslts.Membrane.Nphi{1}, LA_out.Shell_Geom.s{1},...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(LA_out.Stress_Reslts.Membrane.Ntheta{1}, LA_out.Shell_Geom.s{1},...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{5} = 'AQUINAS FE solution meridional MSR';
lgnd{6} = 'AQUINAS FE solution circumferential MSR';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Membrane stress resultants (MSR) in the wall [$N/mm$]','interpreter','latex','fontsize',FS);
ylabel('Arc-length distance from apex of segment [$mm$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');


%% Bending stress resultants
figure('Name','Bending stress resultant plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(4,1);

% Shell bending theory solution
plot(Ms(0:L), 0:L,...
    'Color','r','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
plot(Mth(0:L), 0:L,...
    'Color','b','LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS);
lgnd{1} = 'Bending theory meridional BSR';
lgnd{2} = 'Bending theory circumferential BSR';

% AQUINAS solution
plot(-LA_out.Stress_Reslts.Bending.Mphi{1}, LA_out.Shell_Geom.s{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(-LA_out.Stress_Reslts.Bending.Mtheta{1}, LA_out.Shell_Geom.s{1},...
    'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(-LA_out.Stress_Reslts.Bending.Mphi{1}, LA_out.Shell_Geom.s{1},...
    'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
plot(-LA_out.Stress_Reslts.Bending.Mtheta{1}, LA_out.Shell_Geom.s{1},...
    'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{3} = 'AQUINAS FE solution meridional BSR';
lgnd{4} = 'AQUINAS FE solution circumferential BSR';

legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Bending stress resultants (BSR) in the wall [$Nmm/mm$]','interpreter','latex','fontsize',FS);
ylabel('Arc-length distance from apex of segment [$mm$]','interpreter','latex','fontsize',FS);
set(gca,'YDir','reverse', 'ticklabelinterpreter','latex', 'fontsize',FS, 'tickdir','out');
