%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example LA 02: LA of a thin-walled two-segment isotropic cylindrical shell
% with a midsurface eccentricity. A BC1f boundary condition is considered at
% the base of the cylinder and a BC2f at the top.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% [1] Sadowski A.J. & Rotter J.M. (2012) "Cylindrical shell bending theory
% for orthotropic shells under general axisymmetric pressure distributions"
% Engineering Structures, 42, 258-265.
% http://dx.doi.org/10.1016/j.engstruct.2012.04.024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.

%%%%%%%%%%%%%%%%%
% Preliminaries %
%%%%%%%%%%%%%%%%%
clearvars, clc, close all force
addpath('..','../AQUINAS_Analysis','../AQUINAS_Objects','../Maple_Source');

% Geometry inputs for two-segment cylindrical shell
zB = 0; % [mm] - Meridional location of base of assembly
zM = 0.5*500; % [mm] - Meridional location of junction of both segments
zT = 500; % [mm] - Meridional location of top of assembly

rB = 100; % [mm] - Midsurface radius of bottom segment
tB = 1; % [mm] - Thickness of bottom segment
tT = 1; % [mm] - Thickness of top segment
e = 1; % [mm] - Eccentricity of midsurfaces between top and bottom segments (e > 0, top segment is more inward than bottom one)
rT = rB-e; % [mm] - Midsurface radius of top segment

% Material
MatName = 'Steel'; % Material name
E = 2.0e5; % [N/mm2] - Young's modulus of steel
nu = 0.3; % [-] - Poisson ratio of steel

% Loading inputs
N = 300; % [N/mm] - Meridional compressive shell edge load at zT
pn = 1; % [N/mm2] - Outwards acting normal pressure

% Display inputs
dz = 1; % [mm] - Sampling distance for plotting

% Plot controls
FS = 20; % [-] - Font Size (pts)
FSL = 18; % [-] - Font Size for Legends (pts)
MS = 3; % [-] - Markser size (pts)
LW = 3; % [-] - Line Width (pts)

%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%
M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)));

S1 = AQUINAS_Segment_Object('type','Cone','material',MatName,'thickness',tB,'rzbot',[rB zB],'rztop',[rB zM],'els',100);
S2 = AQUINAS_Segment_Object('type','Cone','material',MatName,'thickness',tT,'rzbot',[rT zM],'rztop',[rT zT],'els',100);

% BC1f condition at base of cylinder
C1 = AQUINAS_Constraint_Object('rcoord',rB,'zcoord',zB,'dofs',{'u'});
C2 = AQUINAS_Constraint_Object('rcoord',rB,'zcoord',zB,'dofs',{'w'});

% BC2f condition at top of cylinder
C3 = AQUINAS_Constraint_Object('rcoord',rT,'zcoord',zT,'dofs',{'u'});

% Eccentricity constraint condition between segments
C4 = AQUINAS_Constraint_Object('type','ecc','segments',{S1,S2},'ends',{'top','bot'});

F1 = AQUINAS_Nodal_Force_Object('type','Fw','rcoord',rT,'zcoord',zT,'magnitude',-N);

P1 = AQUINAS_Distributed_Pressure_Object('segments',{S1,S2},'type','pn','functionHandle',@() pn);

SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',6);

A1 = AQUINAS_Analysis_Object('type','LA');

O1 = AQUINAS_Output_Mode_Object('mode','VarArgOut');

LA_out = AQUINAS_Protocol(M1,S1,S2,C1,C2,C3,C4,F1,P1,SOL1,A1,O1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEYOND THIS PLEASE DO NOT MODIFY UNLESS YOU KNOW WHAT YOU ARE DOING %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%
% Analytical solution %
%%%%%%%%%%%%%%%%%%%%%%%
% Please note a difference in notation below. In AQUINAS, u and w are the
% global displacements parallel to the r and z axes respectively. In usual
% notation used in classical cylindrical shell theory, however, w is the
% displacement normal to the wall and u is the axial displacement.
% Consequently, in the below the meaning of 'w' and 'u' is the OPPOSITE of
% that used in AQUINAS.

if abs(e) > 0.5*(tB+tT); error('Eccentricity is too high.'); end
rT = rB - e; % [mm] - Midsurface radius of top segment
LB = zM - zB; % [mm] - Meridional length of bottom segment
LT = zT - zM; % [mm] - Meridional length of top segment
N0B = N*rT/rB; % [N/mm] - Meridional membrane stress resultant of bottom segment from simple equilibrium relations
lambdaB = pi*sqrt(rB*tB)*(3*(1 - nu*nu))^(-0.25); % [mm] - Linear bending half-wavelength of bottom segment
lambdaT = pi*sqrt(rT*tT)*(3*(1 - nu*nu))^(-0.25); % [mm] - Linear bending half-wavelength of top segment
CB = E*tB*(1-nu*nu)^(-1); CT = E*tT*(1-nu*nu)^(-1); % [N/mm] - Membrane stiffnesses of both segments
DB = E*tB*tB*tB*(12*(1-nu*nu))^(-1); DT = E*tT*tT*tT*(12*(1-nu*nu))^(-1); % [Nmm/mm] - Bending stiffnesses of both segments

% Generic w coefficients in bottom segment - z is a local coordinate on the segment
c1 = @(z,lambda) exp(-pi*z/lambda).*cos(pi*z/lambda);
c2 = @(z,lambda) exp(-pi*z/lambda).*sin(pi*z/lambda);
c3 = @(z,lambda) exp(pi*z/lambda).*cos(pi*z/lambda);
c4 = @(z,lambda) exp(pi*z/lambda).*sin(pi*z/lambda);
c1p = @(z,lambda) -(pi/lambda)*exp(-pi*z/lambda).*(cos(pi*z/lambda)+sin(pi*z/lambda));
c2p = @(z,lambda) (pi/lambda)*exp(-pi*z/lambda).*(cos(pi*z/lambda)-sin(pi*z/lambda));
c3p = @(z,lambda) (pi/lambda)*exp(pi*z/lambda).*(cos(pi*z/lambda)-sin(pi*z/lambda));
c4p = @(z,lambda) (pi/lambda)*exp(pi*z/lambda).*(cos(pi*z/lambda)+sin(pi*z/lambda));
c1pp = @(z,lambda) 2*(pi/lambda)^2*exp(-pi*z/lambda).*sin(pi*z/lambda);
c2pp = @(z,lambda) -2*(pi/lambda)^2*exp(-pi*z/lambda).*cos(pi*z/lambda);
c3pp = @(z,lambda) -2*(pi/lambda)^2*exp(pi*z/lambda).*sin(pi*z/lambda);
c4pp = @(z,lambda) 2*(pi/lambda)^2*exp(pi*z/lambda).*cos(pi*z/lambda);
c1ppp = @(z,lambda) 2*(pi/lambda)^3*exp(-pi*z/lambda).*(cos(pi*z/lambda)-sin(pi*z/lambda));
c2ppp = @(z,lambda) 2*(pi/lambda)^3*exp(-pi*z/lambda).*(cos(pi*z/lambda)+sin(pi*z/lambda));
c3ppp = @(z,lambda) -2*(pi/lambda)^3*exp(pi*z/lambda).*(cos(pi*z/lambda)+sin(pi*z/lambda));
c4ppp = @(z,lambda) 2*(pi/lambda)^3*exp(pi*z/lambda).*(cos(pi*z/lambda)-sin(pi*z/lambda));
c1i = @(z,lambda) 0.5*lambda/pi*(1-exp(-pi*z/lambda).*(cos(pi*z/lambda)-sin(pi*z/lambda)));
c2i = @(z,lambda) 0.5*lambda/pi*(1-exp(-pi*z/lambda).*(cos(pi*z/lambda)+sin(pi*z/lambda)));
c3i = @(z,lambda) 0.5*lambda/pi*(-1+exp(pi*z/lambda).*(sin(pi*z/lambda)+cos(pi*z/lambda)));
c4i = @(z,lambda) 0.5*lambda/pi*(1+exp(pi*z/lambda).*(sin(pi*z/lambda)-cos(pi*z/lambda)));

% BC1f condition at zB
a1_0 = c1(0,lambdaB); a2_0 = c2(0,lambdaB); a3_0 = c3(0,lambdaB); a4_0 = c4(0,lambdaB); wBm_0 = rB*(pn*rB+nu*N0B)/(E*tB);
a1pp_0 = c1pp(0,lambdaB); a2pp_0 = c2pp(0,lambdaB); a3pp_0 = c3pp(0,lambdaB); a4pp_0 = c4pp(0,lambdaB); wBmpp_0 = 0;

% BC2f condition at zT
b1_LT = c1(LT,lambdaT); b2_LT = c2(LT,lambdaT); b3_LT = c3(LT,lambdaT); b4_LT = c4(LT,lambdaT); wTm_LT = rT*(pn*rT+nu*N)/(E*tT);
b1pp_LT = c1pp(LT,lambdaT); b2pp_LT = c2pp(LT,lambdaT); b3pp_LT = c3pp(LT,lambdaT); b4pp_LT = c4pp(LT,lambdaT); wTmpp_LT = 0;

% Full continuity at zM
a1_LB = c1(LB,lambdaB); a2_LB = c2(LB,lambdaB); a3_LB = c3(LB,lambdaB); a4_LB = c4(LB,lambdaB);
a1p_LB = c1p(LB,lambdaB); a2p_LB = c2p(LB,lambdaB); a3p_LB = c3p(LB,lambdaB); a4p_LB = c4p(LB,lambdaB);
a1pp_LB = c1pp(LB,lambdaB); a2pp_LB = c2pp(LB,lambdaB); a3pp_LB = c3pp(LB,lambdaB); a4pp_LB = c4pp(LB,lambdaB);
a1ppp_LB = c1ppp(LB,lambdaB); a2ppp_LB = c2ppp(LB,lambdaB); a3ppp_LB = c3ppp(LB,lambdaB); a4ppp_LB = c4ppp(LB,lambdaB);
b1_0 = c1(0,lambdaT); b2_0 = c2(0,lambdaT); b3_0 = c3(0,lambdaT); b4_0 = c4(0,lambdaT);
b1p_0 = c1p(0,lambdaT); b2p_0 = c2p(0,lambdaT); b3p_0 = c3p(0,lambdaT); b4p_0 = c4p(0,lambdaT);
b1pp_0 = c1pp(0,lambdaT); b2pp_0 = c2pp(0,lambdaT); b3pp_0 = c3pp(0,lambdaT); b4pp_0 = c4pp(0,lambdaT);
b1ppp_0 = c1ppp(0,lambdaT); b2ppp_0 = c2ppp(0,lambdaT); b3ppp_0 = c3ppp(0,lambdaT); b4ppp_0 = c4ppp(0,lambdaT);

% Construction and solution of matrix system
M = [a1_0           a2_0           a3_0           a4_0           0              0              0              0             ;...
     a1pp_0         a2pp_0         a3pp_0         a4pp_0         0              0              0              0             ;...
     a1_LB          a2_LB          a3_LB          a4_LB          -b1_0          -b2_0          -b3_0          -b4_0         ;...
     a1p_LB         a2p_LB         a3p_LB         a4p_LB         -b1p_0         -b2p_0         -b3p_0         -b4p_0        ;...
     rB*DB*a1pp_LB  rB*DB*a2pp_LB  rB*DB*a3pp_LB  rB*DB*a4pp_LB  -rT*DT*b1pp_0  -DT*rT*b2pp_0  -DT*rT*b3pp_0  -DT*rT*b4pp_0 ;...
     a1ppp_LB       a2ppp_LB       a3ppp_LB       a4ppp_LB       -b1ppp_0       -b2ppp_0       -b3ppp_0       -b4ppp_0      ;...
     0              0              0              0              b1_LT          b2_LT          b3_LT          b4_LT         ;...
     0              0              0              0              b1pp_LT        b2pp_LT        b3pp_LT        b4pp_LT       ];
V = [-wBm_0; -wBmpp_0; 0; 0; -rT*N*e; 0; -wTm_LT; -wTmpp_LT];
disp('Ignore any warning relating to matrix ill-conditioning - it does not relate to AQUINAS.');
Cs = M\V; A = Cs(1:4); B = Cs(5:8);


% Post-processing
% Bottom segment
zloc = (zB:dz:zM) - zB;
wB = A(1)*c1(zloc,lambdaB) + A(2)*c2(zloc,lambdaB) + A(3)*c3(zloc,lambdaB) + A(4)*c4(zloc,lambdaB) + rB*(pn*rB+nu*N0B)/(E*tB);
wBp = A(1)*c1p(zloc,lambdaB) + A(2)*c2p(zloc,lambdaB) + A(3)*c3p(zloc,lambdaB) + A(4)*c4p(zloc,lambdaB);
wBpp = A(1)*c1pp(zloc,lambdaB) + A(2)*c2pp(zloc,lambdaB) + A(3)*c3pp(zloc,lambdaB) + A(4)*c4pp(zloc,lambdaB);
wmB = ones(size(wB)).*rB*(pn*rB+nu*N0B)/(E*tB);

uB = -zloc*N0B/CB - (nu/rB)*(A(1)*c1i(zloc,lambdaB) + A(2)*c2i(zloc,lambdaB) + A(3)*c3i(zloc,lambdaB) + A(4)*c4i(zloc,lambdaB) + zloc*rB*(pn*rB+nu*N0B)/(E*tB));
uBp = -N0B/CB-nu*wB/rB;
umB = -zloc*(N0B+nu*pn*rB)/(E*tB);

MzB = -wBpp*DB; sigmabzBout = 6*MzB/(tB*tB); sigmabzBin = -sigmabzBout;
MthB = nu*MzB; sigmabthBout = 6*MthB/(tB*tB); sigmabthBin = -sigmabthBout;

NzB = CB*(uBp+nu*wB/rB); sigmazBmid = NzB/tB; sigmazBout = sigmabzBout + sigmazBmid; sigmazBin = sigmabzBin + sigmazBmid;
NzmB = -ones(size(NzB)).*N0B; sigmazmBmid = NzmB/tB;

NthB = CB*(wB/rB+nu*uBp); sigmathBmid = NthB/tB; sigmathBout = sigmabthBout + sigmathBmid; sigmathBin = sigmabthBin + sigmathBmid;
NthmB = ones(size(NzB))*pn*rB; sigmathmBmid = NthmB/tB;

sigmavmBout = sqrt(sigmazBout.*sigmazBout - sigmazBout.*sigmathBout + sigmathBout.*sigmathBout);
sigmavmBmid = sqrt(sigmazBmid.*sigmazBmid - sigmazBmid.*sigmathBmid + sigmathBmid.*sigmathBmid);
sigmavmBin = sqrt(sigmazBin.*sigmazBin - sigmazBin.*sigmathBin + sigmathBin.*sigmathBin);

% Top segment
zloc = (zM:dz:zT) - zM;
wT = B(1)*c1(zloc,lambdaT) + B(2)*c2(zloc,lambdaT) + B(3)*c3(zloc,lambdaT) + B(4)*c4(zloc,lambdaT) + rT*(pn*rT+nu*N)/(E*tT);
wTp = B(1)*c1p(zloc,lambdaT) + B(2)*c2p(zloc,lambdaT) + B(3)*c3p(zloc,lambdaT) + B(4)*c4p(zloc,lambdaT);
wTpp = B(1)*c1pp(zloc,lambdaT) + B(2)*c2pp(zloc,lambdaT) + B(3)*c3pp(zloc,lambdaT) + B(4)*c4pp(zloc,lambdaT);
wmT = ones(size(wT)).*rT*(pn*rT+nu*N)/(E*tT);

uT = -zloc*N/CT - (nu/rT)*(B(1)*c1i(zloc,lambdaT) + B(2)*c2i(zloc,lambdaT) + B(3)*c3i(zloc,lambdaT) + B(4)*c4i(zloc,lambdaT) + zloc*rT*(pn*rT+nu*N)/(E*tT)) + uB(end) + wBp(end)*e;
uTp =- N/CT-nu*wT/rT;
umT = -LB*(N+nu*pn*rB)/(E*tB) - zloc*(N+nu*pn*rT)/(E*tT);

MzT = -wTpp*DT; sigmabzTout = 6*MzT/(tT*tT); sigmabzTin = -sigmabzTout;
MthT = nu*MzT; sigmabthTout = 6*MthT/(tT*tT); sigmabthTin = -sigmabthTout;

NzT = CT*(uTp+nu*wT/rT); sigmazTmid = NzT/tT; sigmazTout = sigmabzTout + sigmazTmid; sigmazTin = sigmabzTin + sigmazTmid;
NzmT = -ones(size(NzT)).*N; sigmazmTmid = NzmT/tT;

NthT = CT*(wT/rT+nu*uTp); sigmathTmid = NthT/tT; sigmathTout = sigmabthTout + sigmathTmid; sigmathTin = sigmabthTin + sigmathTmid;
NthmT = ones(size(NzT))*pn*rT; sigmathmTmid = NthmT/tT;

sigmavmTout = sqrt(sigmazTout.*sigmazTout - sigmazTout.*sigmathTout + sigmathTout.*sigmathTout);
sigmavmTmid = sqrt(sigmazTmid.*sigmazTmid - sigmazTmid.*sigmathTmid + sigmathTmid.*sigmathTmid);
sigmavmTin = sqrt(sigmazTin.*sigmazTin - sigmazTin.*sigmathTin + sigmathTin.*sigmathTin);


%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%

%% Radial / normal displacement plot
figure('Name','Radial displacement plot','NumberTitle','off','WindowState','Maximized','color','w');
pH = cell(6,1);
lgnd = cell(6,1);

% Shell membrane theory solution
pH{1} = plot(wmT,zM:dz:zT,'Color',[1.0 0.5 0.5],'LineStyle','--','LineWidth',2*LW);  hold on; grid on;
pH{2} = plot(wmB,zB:dz:zM,'Color',[0.5 0.5 1.0],'LineStyle','--','LineWidth',2*LW);
lgnd{1} = 'Top segment membrane theory';
lgnd{2} = 'Bottom segment membrane theory';

% Shell bending theory solution
pH{3} = plot(wT,zM:dz:zT,'Color','r','LineStyle','-','LineWidth',LW);
pH{4} = plot(wB,zB:dz:zM,'Color','b','LineStyle','-','LineWidth',LW);
lgnd{3} = 'Top segment bending theory';
lgnd{4} = 'Bottom segment bending theory';

% AQUINAS solution
pH{5} = plot(LA_out.DOFs.u{1}, LA_out.Shell_Geom.z{1},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.DOFs.u{1}, LA_out.Shell_Geom.z{1},...
     'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
pH{6} = plot(LA_out.DOFs.u{2}, LA_out.Shell_Geom.z{2},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.DOFs.u{2}, LA_out.Shell_Geom.z{2},...
     'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{5} = 'Top segment AQUINAS FE solution';
lgnd{6} = 'Bottom segment AQUINAS FE solution';

set(gca,'FontSize',FS,'FontName','Times');
xlabel('Wall radial displacement $w$ [mm]','interpreter','latex');
ylabel('Meridional position $z$ [mm]','interpreter','latex');
legend([pH{1} pH{2} pH{3} pH{4} pH{5} pH{6}],lgnd,'interpreter','latex','fontsize',FSL,'location','best');

%% Axial / meridional displacement plot
figure('Name','Meridional displacement plot','NumberTitle','off','WindowState','Maximized','color','w');
pH = cell(6,1);
lgnd = cell(6,1);

% Shell membrane theory solution
pH{1} = plot(umT,zM:dz:zT,'Color',[1.0 0.5 0.5],'LineStyle','--','LineWidth',2*LW);  hold on; grid on;
pH{2} = plot(umB,zB:dz:zM,'Color',[0.5 0.5 1.0],'LineStyle','--','LineWidth',2*LW);
lgnd{1} = 'Top segment membrane theory';
lgnd{2} = 'Bottom segment membrane theory';

% Shell bending theory solution
pH{3} = plot(uT,zM:dz:zT,'Color','r','LineStyle','-','LineWidth',LW);
pH{4} = plot(uB,zB:dz:zM,'Color','b','LineStyle','-','LineWidth',LW);
lgnd{3} = 'Top segment bending theory';
lgnd{4} = 'Bottom segment bending theory';

% AQUINAS solution
pH{5} = plot(LA_out.DOFs.w{1}, LA_out.Shell_Geom.z{1},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.DOFs.w{1}, LA_out.Shell_Geom.z{1},...
     'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
pH{6} = plot(LA_out.DOFs.w{2}, LA_out.Shell_Geom.z{2},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.DOFs.w{2}, LA_out.Shell_Geom.z{2},...
     'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{5} = 'Top segment AQUINAS FE solution';
lgnd{6} = 'Bottom segment AQUINAS FE solution';

set(gca,'FontSize',FS,'FontName','Times');
xlabel('Wall meridional displacement $u$ [mm]','interpreter','latex');
ylabel('Meridional position $z$ [mm]','interpreter','latex');
legend([pH{1} pH{2} pH{3} pH{4} pH{5} pH{6}],lgnd,'interpreter','latex','fontsize',FSL,'location','best');


%% Deformed shape plot
figure('Name','Deformed shape plot','NumberTitle','off','WindowState','Maximized','color','w');
pH = cell(7,1);
lgnd = cell(7,1);

% Original geometry of segments
pH{1} = plot([rT rT],[zM zT],'Color','k','LineStyle','-','LineWidth',LW); hold on; grid on;
plot([rB rB],[zB zM],'Color','k','LineStyle','-','LineWidth',LW);
lgnd{1} = 'Original geometry of shell segment';

% Shell membrane theory solution
pH{2} = plot(rT+wmT,(zM:dz:zT)+umT,'Color',[1.0 0.5 0.5],'LineStyle','--','LineWidth',2*LW);
pH{3} = plot(rB+wmB,(zB:dz:zM)+umB,'Color',[0.5 0.5 1.0],'LineStyle','--','LineWidth',2*LW);
lgnd{2} = 'Top segment membrane theory';
lgnd{3} = 'Bottom segment membrane theory';

% Shell bending theory solution
pH{4} = plot(rT+wT,(zM:dz:zT)+uT,'Color','r','LineStyle','-','LineWidth',LW);
pH{5} = plot(rB+wB,(zB:dz:zM)+uB,'Color','b','LineStyle','-','LineWidth',LW);
lgnd{4} = 'Top segment bending theory';
lgnd{5} = 'Bottom segment bending theory';

% AQUINAS solution
pH{6} = plot(LA_out.Shell_Geom.r{1}+LA_out.DOFs.u{1}, LA_out.Shell_Geom.z{1}+LA_out.DOFs.w{1},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Shell_Geom.r{1}+LA_out.DOFs.u{1}, LA_out.Shell_Geom.z{1}+LA_out.DOFs.w{1},...
     'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
pH{7} = plot(LA_out.Shell_Geom.r{2}+LA_out.DOFs.u{2}, LA_out.Shell_Geom.z{2}+LA_out.DOFs.w{2},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Shell_Geom.r{2}+LA_out.DOFs.u{2}, LA_out.Shell_Geom.z{2}+LA_out.DOFs.w{2},...
     'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{6} = 'Top segment AQUINAS FE solution';
lgnd{7} = 'Bottom segment AQUINAS FE solution';

set(gca,'FontSize',FS,'FontName','Times');
xlabel('Radial deformed wall position $r+w$ [mm]','interpreter','latex');
ylabel('Meridional deformed wall position $z+u$ [mm]','interpreter','latex');
legend([pH{1} pH{2} pH{3} pH{4} pH{5} pH{6} pH{7}],lgnd,'interpreter','latex','fontsize',FSL,'location','best');


%% Meridional rotation plot
figure('Name','Meridional rotation plot','NumberTitle','off','WindowState','Maximized','color','w');
pH = cell(4,1);
lgnd = cell(4,1);

% Shell bending theory solution
pH{1} = plot(wTp,zM:dz:zT,'Color','r','LineStyle','-','LineWidth',LW); hold on; grid on;
pH{2} = plot(wBp,zB:dz:zM,'Color','b','LineStyle','-','LineWidth',LW);
lgnd{1} = 'Top segment bending theory';
lgnd{2} = 'Bottom segment bending theory';

% AQUINAS solution
pH{3} = plot(LA_out.DOFs.beta{1}, LA_out.Shell_Geom.z{1},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.DOFs.beta{1}, LA_out.Shell_Geom.z{1},...
     'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
pH{4} = plot(LA_out.DOFs.beta{2}, LA_out.Shell_Geom.z{2},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.DOFs.beta{2}, LA_out.Shell_Geom.z{2},...
     'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{3} = 'Top segment AQUINAS FE solution';
lgnd{4} = 'Bottom segment AQUINAS FE solution';

set(gca,'FontSize',FS,'FontName','Times');
xlabel('Wall meridional rotation $\beta = \frac{dw}{dz}$ [-]','interpreter','latex');
ylabel('Meridional position $z$ [mm]','interpreter','latex');
legend([pH{1} pH{2} pH{3} pH{4}],lgnd,'interpreter','latex','fontsize',FSL,'location','best');


%% Meridional stress plot
figure('Name','Meridional stress plot','NumberTitle','off','WindowState','Maximized','color','w');
pH = cell(14,1);
lgnd = cell(14,1);

% Shell membrane theory solution
pH{1} = plot(sigmazmTmid,zM:dz:zT,'Color',[1.0 0.5 0.5],'LineStyle','--','LineWidth',2*LW); grid on; hold on;
pH{2} = plot(sigmazmBmid,zB:dz:zM,'Color',[0.5 0.5 1.0],'LineStyle','--','LineWidth',2*LW);
lgnd{1} = 'Top segment membrane theory midsurface';
lgnd{2} = 'Bottom segment membrane theory midsurface';

% Shell bending theory solution
pH{3} = plot(sigmazTin,zM:dz:zT,'Color','r','LineStyle','-.','LineWidth',LW);
pH{4} = plot(sigmazTmid,zM:dz:zT,'Color','r','LineStyle','-','LineWidth',LW);
pH{5} = plot(sigmazTout,zM:dz:zT,'Color','r','LineStyle',':','LineWidth',LW);
pH{6} = plot(sigmazBin,zB:dz:zM,'Color','b','LineStyle','-.','LineWidth',LW);
pH{7} = plot(sigmazBmid,zB:dz:zM,'Color','b','LineStyle','-','LineWidth',LW);
pH{8} = plot(sigmazBout,zB:dz:zM,'Color','b','LineStyle',':','LineWidth',LW);
lgnd{3} = 'Top segment bending theory inner surface stress';
lgnd{4} = 'Top segment bending theory midsurface stress';
lgnd{5} = 'Top segment bending theory outter surface stress';
lgnd{6} = 'Bottom segment bending theory inner surface stress';
lgnd{7} = 'Bottom segment bending theory midsurface stress';
lgnd{8} = 'Bottom segment bending theory outter surface stress';

% AQUINAS solution
pH{9} = plot(LA_out.Stresses.Sig_phi.i{1}, LA_out.Shell_Geom.z{1},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stresses.Sig_phi.i{1}, LA_out.Shell_Geom.z{1},...
     'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
pH{10} = plot(LA_out.Stresses.Sig_phi.m{1}, LA_out.Shell_Geom.z{1},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','p','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stresses.Sig_phi.m{1}, LA_out.Shell_Geom.z{1},...
     'LineStyle','none','Marker','p','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
pH{11} = plot(LA_out.Stresses.Sig_phi.o{1}, LA_out.Shell_Geom.z{1},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stresses.Sig_phi.o{1}, LA_out.Shell_Geom.z{1},...
     'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
pH{12} = plot(LA_out.Stresses.Sig_phi.i{2}, LA_out.Shell_Geom.z{2},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stresses.Sig_phi.i{2}, LA_out.Shell_Geom.z{2},...
     'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
pH{13} = plot(LA_out.Stresses.Sig_phi.m{2}, LA_out.Shell_Geom.z{2},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','p','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stresses.Sig_phi.m{2}, LA_out.Shell_Geom.z{2},...
     'LineStyle','none','Marker','p','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
pH{14} = plot(LA_out.Stresses.Sig_phi.o{2}, LA_out.Shell_Geom.z{2},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stresses.Sig_phi.o{2}, LA_out.Shell_Geom.z{2},...
     'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{9} = 'Top segment AQUINAS inner surface stress';
lgnd{10} = 'Top segment AQUINAS midsurface stress';
lgnd{11} = 'Top segment AQUINAS outter surface stress';
lgnd{12} = 'Bottom segment AQUINAS inner surface stress';
lgnd{13} = 'Bottom segment AQUINAS midsurface stress';
lgnd{14} = 'Bottom segment AQUINAS outter surface stress';

set(gca,'FontSize',FS,'FontName','Times');
xlabel('Meridional stresses $\sigma_{\phi}$ $[N/mm^2]$','interpreter','latex');
ylabel('Meridional position $z$ [mm]','interpreter','latex');
legend([pH{1} pH{2} pH{3} pH{4} pH{5} pH{6} pH{7} pH{8} pH{9} pH{10} pH{11} pH{12} pH{13} pH{14}],lgnd,'interpreter','latex','fontsize',FSL,'location','best');

%% Circumferential stress plot
figure('Name','Circumferential stress plot','NumberTitle','off','WindowState','Maximized','color','w');
pH = cell(14,1);
lgnd = cell(14,1);

% Shell membrane theory solution
pH{1} = plot(sigmathmTmid,zM:dz:zT,'Color',[1.0 0.5 0.5],'LineStyle','--','LineWidth',2*LW); grid on; hold on;
pH{2} = plot(sigmathmBmid,zB:dz:zM,'Color',[0.5 0.5 1.0],'LineStyle','--','LineWidth',2*LW);
lgnd{1} = 'Top segment membrane theory midsurface';
lgnd{2} = 'Bottom segment membrane theory midsurface';

% Shell bending theory solution
pH{3} = plot(sigmathTin,zM:dz:zT,'Color','r','LineStyle','-.','LineWidth',LW);
pH{4} = plot(sigmathTmid,zM:dz:zT,'Color','r','LineStyle','-','LineWidth',LW);
pH{5} = plot(sigmathTout,zM:dz:zT,'Color','r','LineStyle',':','LineWidth',LW);
pH{6} = plot(sigmathBin,zB:dz:zM,'Color','b','LineStyle','-.','LineWidth',LW);
pH{7} = plot(sigmathBmid,zB:dz:zM,'Color','b','LineStyle','-','LineWidth',LW);
pH{8} = plot(sigmathBout,zB:dz:zM,'Color','b','LineStyle',':','LineWidth',LW);
lgnd{3} = 'Top segment bending theory inner surface stress';
lgnd{4} = 'Top segment bending theory midsurface stress';
lgnd{5} = 'Top segment bending theory outter surface stress';
lgnd{6} = 'Bottom segment bending theory inner surface stress';
lgnd{7} = 'Bottom segment bending theory midsurface stress';
lgnd{8} = 'Bottom segment bending theory outter surface stress';

% AQUINAS solution
pH{9} = plot(LA_out.Stresses.Sig_theta.i{1}, LA_out.Shell_Geom.z{1},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stresses.Sig_theta.i{1}, LA_out.Shell_Geom.z{1},...
     'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
pH{10} = plot(LA_out.Stresses.Sig_theta.m{1}, LA_out.Shell_Geom.z{1},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','p','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stresses.Sig_theta.m{1}, LA_out.Shell_Geom.z{1},...
     'LineStyle','none','Marker','p','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
pH{11} = plot(LA_out.Stresses.Sig_theta.o{1}, LA_out.Shell_Geom.z{1},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stresses.Sig_theta.o{1}, LA_out.Shell_Geom.z{1},...
     'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
pH{12} = plot(LA_out.Stresses.Sig_theta.i{2}, LA_out.Shell_Geom.z{2},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stresses.Sig_theta.i{2}, LA_out.Shell_Geom.z{2},...
     'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
pH{13} = plot(LA_out.Stresses.Sig_theta.m{2}, LA_out.Shell_Geom.z{2},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','p','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stresses.Sig_theta.m{2}, LA_out.Shell_Geom.z{2},...
     'LineStyle','none','Marker','p','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
pH{14} = plot(LA_out.Stresses.Sig_theta.o{2}, LA_out.Shell_Geom.z{2},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stresses.Sig_theta.o{2}, LA_out.Shell_Geom.z{2},...
     'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{9} = 'Top segment AQUINAS inner surface stress';
lgnd{10} = 'Top segment AQUINAS midsurface stress';
lgnd{11} = 'Top segment AQUINAS outter surface stress';
lgnd{12} = 'Bottom segment AQUINAS inner surface stress';
lgnd{13} = 'Bottom segment AQUINAS midsurface stress';
lgnd{14} = 'Bottom segment AQUINAS outter surface stress';

set(gca,'FontSize',FS,'FontName','Times');
xlabel('Circumferential stresses $\sigma_{\theta}$ $[N/mm^2]$','interpreter','latex');
ylabel('Meridional position $z$ [mm]','interpreter','latex');
legend([pH{1} pH{2} pH{3} pH{4} pH{5} pH{6} pH{7} pH{8} pH{9} pH{10} pH{11} pH{12} pH{13} pH{14}],lgnd,'interpreter','latex','fontsize',FSL,'location','best');


%% Von Mises stress plot
figure('Name','Von Mises stress plot','NumberTitle','off','WindowState','Maximized','color','w');
pH = cell(12,1);
lgnd = cell(12,1);

% Shell bending theory solution
pH{1} = plot(sigmavmTin,zM:dz:zT,'Color','r','LineStyle','-.','LineWidth',LW); hold on; grid on;
pH{2} = plot(sigmavmTmid,zM:dz:zT,'Color','r','LineStyle','-','LineWidth',LW);
pH{3} = plot(sigmavmTout,zM:dz:zT,'Color','r','LineStyle',':','LineWidth',LW);
pH{4} = plot(sigmavmBin,zB:dz:zM,'Color','b','LineStyle','-.','LineWidth',LW);
pH{5} = plot(sigmavmBmid,zB:dz:zM,'Color','b','LineStyle','-','LineWidth',LW);
pH{6} = plot(sigmavmBout,zB:dz:zM,'Color','b','LineStyle',':','LineWidth',LW);
lgnd{1} = 'Top segment bending theory inner surface stress';
lgnd{2} = 'Top segment bending theory midsurface stress';
lgnd{3} = 'Top segment bending theory outter surface stress';
lgnd{4} = 'Bottom segment bending theory inner surface stress';
lgnd{5} = 'Bottom segment bending theory midsurface stress';
lgnd{6} = 'Bottom segment bending theory outter surface stress';

% AQUINAS solution
pH{7} = plot(LA_out.Stresses.Sig_vonMises.i{1}, LA_out.Shell_Geom.z{1},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stresses.Sig_vonMises.i{1}, LA_out.Shell_Geom.z{1},...
     'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
pH{8} = plot(LA_out.Stresses.Sig_vonMises.m{1}, LA_out.Shell_Geom.z{1},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','p','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stresses.Sig_vonMises.m{1}, LA_out.Shell_Geom.z{1},...
     'LineStyle','none','Marker','p','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
pH{9} = plot(LA_out.Stresses.Sig_vonMises.o{1}, LA_out.Shell_Geom.z{1},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stresses.Sig_vonMises.o{1}, LA_out.Shell_Geom.z{1},...
     'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
pH{10} = plot(LA_out.Stresses.Sig_vonMises.i{2}, LA_out.Shell_Geom.z{2},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stresses.Sig_vonMises.i{2}, LA_out.Shell_Geom.z{2},...
     'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
pH{11} = plot(LA_out.Stresses.Sig_vonMises.m{2}, LA_out.Shell_Geom.z{2},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','p','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stresses.Sig_vonMises.m{2}, LA_out.Shell_Geom.z{2},...
     'LineStyle','none','Marker','p','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
pH{12} = plot(LA_out.Stresses.Sig_vonMises.o{2}, LA_out.Shell_Geom.z{2},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stresses.Sig_vonMises.o{2}, LA_out.Shell_Geom.z{2},...
     'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{7} = 'Top segment AQUINAS inner surface stress';
lgnd{8} = 'Top segment AQUINAS midsurface stress';
lgnd{9} = 'Top segment AQUINAS outter surface stress';
lgnd{10} = 'Bottom segment AQUINAS inner surface stress';
lgnd{11} = 'Bottom segment AQUINAS midsurface stress';
lgnd{12} = 'Bottom segment AQUINAS outter surface stress';

set(gca,'FontSize',FS,'FontName','Times');
xlabel('Von Mises stresses $\sigma_{vM}$ $[N/mm^2]$','interpreter','latex');
ylabel('Meridional position $z$ [mm]','interpreter','latex');
legend([pH{1} pH{2} pH{3} pH{4} pH{5} pH{6} pH{7} pH{8} pH{9} pH{10} pH{11} pH{12}],lgnd,'interpreter','latex','fontsize',FSL,'location','best');

%% Membrane stress resultant plot
figure('Name','Membrane stress resultant plot','NumberTitle','off','WindowState','Maximized','color','w');
pH = cell(12,1);
lgnd = cell(12,1);

% Shell membrane theory solution
pH{1} = plot(NzmT,zM:dz:zT,'Color',[1.0 0.5 0.5],'LineStyle','-','LineWidth',LW); hold on; grid on;
pH{2} = plot(NthmT,zM:dz:zT,'Color',[1.0 0.5 0.5],'LineStyle','--','LineWidth',LW);
pH{3} = plot(NzmB,zB:dz:zM,'Color',[0.5 0.5 1.0],'LineStyle','-','LineWidth',LW);
pH{4} = plot(NthmB,zB:dz:zM,'Color',[0.5 0.5 1.0],'LineStyle','--','LineWidth',LW);
lgnd{1} = 'Top segment membrane theory meridional MSR';
lgnd{2} = 'Top segment membrane theory circumferential MSR';
lgnd{3} = 'Bottom segment membrane theory meridional MSR';
lgnd{4} = 'Bottom segment membrane theory circumferential MSR';

% Shell bending theory solution
pH{5} = plot(NzT,zM:dz:zT,'Color','r','LineStyle','-','LineWidth',LW);
pH{6} = plot(NthT,zM:dz:zT,'Color','r','LineStyle','--','LineWidth',LW);
pH{7} = plot(NzB,zB:dz:zM,'Color','b','LineStyle','-','LineWidth',LW);
pH{8} = plot(NthB,zB:dz:zM,'Color','b','LineStyle','--','LineWidth',LW);
lgnd{5} = 'Top segment bending theory meridional MSR';
lgnd{6} = 'Top segment bending theory circumferential MSR';
lgnd{7} = 'Bottom segment bending theory meridional MSR';
lgnd{8} = 'Bottom segment bending theory circumferential MSR';

% AQUINAS solution
pH{9} = plot(LA_out.Stress_Reslts.Membrane.Nphi{1}, LA_out.Shell_Geom.z{1},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stress_Reslts.Membrane.Nphi{1}, LA_out.Shell_Geom.z{1},...
     'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
pH{10} = plot(LA_out.Stress_Reslts.Membrane.Ntheta{1}, LA_out.Shell_Geom.z{1},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stress_Reslts.Membrane.Ntheta{1}, LA_out.Shell_Geom.z{1},...
     'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
pH{11} = plot(LA_out.Stress_Reslts.Membrane.Nphi{2}, LA_out.Shell_Geom.z{2},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stress_Reslts.Membrane.Nphi{2}, LA_out.Shell_Geom.z{2},...
     'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
pH{12} = plot(LA_out.Stress_Reslts.Membrane.Ntheta{2}, LA_out.Shell_Geom.z{2},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stress_Reslts.Membrane.Ntheta{2}, LA_out.Shell_Geom.z{2},...
     'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{9} = 'Top segment AQUINAS meridional MSR';
lgnd{10} = 'Top segment AQUINAS circumferential MSR';
lgnd{11} = 'Bottom segment AQUINAS meridional MSR';
lgnd{12} = 'Bottom segment AQUINAS circumferential MSR';

set(gca,'FontSize',FS,'FontName','Times');
xlabel('Membrane Stress Resultant (MSR) [N/mm]','interpreter','latex');
ylabel('Meridional position $z$ [mm]','interpreter','latex');
legend([pH{1} pH{2} pH{3} pH{4} pH{5} pH{6} pH{7} pH{8} pH{9} pH{10} pH{11} pH{12}],lgnd,'interpreter','latex','fontsize',FSL,'location','best');


%% Bending stress resultant plot
figure('Name','Bending stress resultant plot','NumberTitle','off','WindowState','Maximized','color','w');
pH = cell(8,1);
lgnd = cell(8,1);

% Shell bending theory solution
pH{1} = plot(MzT,zM:dz:zT,'Color','r','LineStyle','-','LineWidth',LW); hold on; grid on;
pH{2} = plot(MthT,zM:dz:zT,'Color','r','LineStyle','--','LineWidth',LW);
pH{3} = plot(MzB,zB:dz:zM,'Color','b','LineStyle','-','LineWidth',LW);
pH{4} = plot(MthB,zB:dz:zM,'Color','b','LineStyle','--','LineWidth',LW);
lgnd{1} = 'Top segment bending theory meridional BSR';
lgnd{2} = 'Top segment bending theory circumferential BSR';
lgnd{3} = 'Bottom segment bending theory meridional BSR';
lgnd{4} = 'Bottom segment bending theory circumferential BSR';


% AQUINAS solution
pH{5} = plot(LA_out.Stress_Reslts.Bending.Mphi{1}, LA_out.Shell_Geom.z{1},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stress_Reslts.Bending.Mphi{1}, LA_out.Shell_Geom.z{1},...
     'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
pH{6} = plot(LA_out.Stress_Reslts.Bending.Mtheta{1}, LA_out.Shell_Geom.z{1},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
plot(LA_out.Stress_Reslts.Bending.Mtheta{1}, LA_out.Shell_Geom.z{1},...
     'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',0.8*MS);
pH{7} = plot(LA_out.Stress_Reslts.Bending.Mphi{2}, LA_out.Shell_Geom.z{2},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stress_Reslts.Bending.Mphi{2}, LA_out.Shell_Geom.z{2},...
     'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
pH{8} = plot(LA_out.Stress_Reslts.Bending.Mtheta{2}, LA_out.Shell_Geom.z{2},...
     'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
plot(LA_out.Stress_Reslts.Bending.Mtheta{2}, LA_out.Shell_Geom.z{2},...
     'LineStyle','none','Marker','d','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',0.8*MS);
lgnd{5} = 'Top segment AQUINAS meridional BSR';
lgnd{6} = 'Top segment AQUINAS circumferential BSR';
lgnd{7} = 'Bottom segment AQUINAS meridional BSR';
lgnd{8} = 'Bottom segment AQUINAS circumferential BSR';

set(gca,'FontSize',FS,'FontName','Times');
xlabel('Bending Stress Resultant (BSR) [Nmm/mm]','interpreter','latex');
ylabel('Meridional position $z$ [mm]','interpreter','latex');
legend([pH{1} pH{2} pH{3} pH{4} pH{5} pH{6} pH{7} pH{8}],lgnd,'interpreter','latex','fontsize',FSL,'location','best');


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