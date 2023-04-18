function BDTB = AQUINAS_BDTB_matrix(eta,theta,n,dat,DT)
% Last generated on Thu Mar 09 11:44:28 2023

eta2 = eta*eta; eta3 = eta2*eta; L = dat.L;
cn = cos(n*theta); sn = sin(n*theta);
N01 = (2-3*eta+eta3)*0.25; N11 = L*(1-eta-eta2+eta3)*0.25; N02 = (2+3*eta-eta3)*0.25; N12 = L*(-1-eta+eta2+eta3)*0.25; % Eq. 2 in Rotter & Teng (Vol. 31)
r = N01*dat.r1 + N02*dat.r2 + N11*dat.drds1 + N12*dat.drds2; % Eqs. 1a & 2 in Rotter & Teng (Vol. 31)
phi = N01*dat.phi1 + N02*dat.phi2 + N11*dat.dphids1 + N12*dat.dphids2; % Cubic interpolation to compute the angle phi at the current Gauss station
dphids = N01*dat.dphids1 + N02*dat.dphids2 + N11*dat.d2phids2_1 + N12*dat.d2phids2_2; % Eqs 1c & 2 in Rotter & Teng (Vol. 31)
c = cos(phi); s = sin(phi);

U1 = dat.elDOFs(1);
dUdS_1 = dat.elDOFs(2);
W1 = dat.elDOFs(5);
dWdS_1 = dat.elDOFs(6);
U2 = dat.elDOFs(7);
dUdS_2 = dat.elDOFs(8);
W2 = dat.elDOFs(11);
dWdS_2 = dat.elDOFs(12);

B0 = zeros(6,12);
G = zeros(6,12);

B0_1_1 = 0.3e1 / 0.4e1 * c * cn * (-1 + eta2) / L;
B0_1_2 = (c * cn * (-2 * eta - 1 + 3 * eta2)) / 0.4e1;
B0_1_3 = 0;
B0_1_4 = 0;
B0_1_5 = -0.3e1 / 0.4e1 * s * cn * (-1 + eta2) / L;
B0_1_6 = -(s * cn * (-2 * eta - 1 + 3 * eta2)) / 0.4e1;
B0_1_7 = -0.3e1 / 0.4e1 * c * cn * (-1 + eta2) / L;
B0_1_8 = (c * cn * (2 * eta - 1 + 3 * eta2)) / 0.4e1;
B0_1_9 = 0;
B0_1_10 = 0;
B0_1_11 = 0.3e1 / 0.4e1 * s * cn * (-1 + eta2) / L;
B0_1_12 = -(s * cn * (2 * eta - 1 + 3 * eta2)) / 0.4e1;

B0_2_1 = (cn * (eta + 2) * (eta - 1) ^ 2 * (c ^ 2 + s ^ 2) / r) / 0.4e1;
B0_2_2 = (cn * L * (eta + 1) * (eta - 1) ^ 2 * (c ^ 2 + s ^ 2) / r) / 0.4e1;
B0_2_3 = (n * cn * (-3 * eta + 2 + eta3) / r) / 0.4e1;
B0_2_4 = (n * cn * L * (eta + 1) * (eta - 1) ^ 2 / r) / 0.4e1;
B0_2_5 = 0;
B0_2_6 = 0;
B0_2_7 = -(cn * (eta - 2) * (eta + 1) ^ 2 * (c ^ 2 + s ^ 2) / r) / 0.4e1;
B0_2_8 = (cn * L * (eta - 1) * (eta + 1) ^ 2 * (c ^ 2 + s ^ 2) / r) / 0.4e1;
B0_2_9 = -(n * cn * (eta - 2) * (eta + 1) ^ 2 / r) / 0.4e1;
B0_2_10 = (n * cn * L * (eta - 1) * (eta + 1) ^ 2 / r) / 0.4e1;
B0_2_11 = 0;
B0_2_12 = 0;

B0_3_1 = -(c * n * sn * (eta + 2) * (eta - 1) ^ 2 / r) / 0.4e1;
B0_3_2 = -(c * n * sn * L * (eta + 1) * (eta - 1) ^ 2 / r) / 0.4e1;
B0_3_3 = -(sn * (L * c * eta2 + (c * L - 3 * r) * eta - 2 * c * L - 3 * r) * (eta - 1) / r / L) / 0.4e1;
B0_3_4 = -(sn * ((-3 * eta - 1) * r + L * c * (eta - 1) * (eta + 1)) * (eta - 1) / r) / 0.4e1;
B0_3_5 = (s * n * sn * (eta + 2) * (eta - 1) ^ 2 / r) / 0.4e1;
B0_3_6 = (s * n * sn * L * (eta + 1) * (eta - 1) ^ 2 / r) / 0.4e1;
B0_3_7 = (c * n * sn * (eta - 2) * (eta + 1) ^ 2 / r) / 0.4e1;
B0_3_8 = -(c * n * sn * L * (eta - 1) * (eta + 1) ^ 2 / r) / 0.4e1;
B0_3_9 = (sn * (L * c * eta2 + (-c * L - 3 * r) * eta - 2 * c * L + 3 * r) * (eta + 1) / r / L) / 0.4e1;
B0_3_10 = -(((-3 * eta + 1) * r + L * c * (eta - 1) * (eta + 1)) * sn * (eta + 1) / r) / 0.4e1;
B0_3_11 = -(s * n * sn * (eta - 2) * (eta + 1) ^ 2 / r) / 0.4e1;
B0_3_12 = (s * n * sn * L * (eta - 1) * (eta + 1) ^ 2 / r) / 0.4e1;

B0_4_1 = -0.3e1 / 0.4e1 * (c * dphids * (eta - 1) * (eta + 1) * L + 2 * eta * s) * cn / (L ^ 2);
B0_4_2 = -0.3e1 / 0.4e1 * ((eta + 0.1e1 / 0.3e1) * c * dphids * (eta - 0.1e1) * L + 0.2e1 * s * (eta - 0.1e1 / 0.3e1)) * cn / L;
B0_4_3 = 0;
B0_4_4 = 0;
B0_4_5 = 0.3e1 / 0.4e1 * (s * dphids * (eta - 1) * (eta + 1) * L - 2 * c * eta) * cn / (L ^ 2);
B0_4_6 = 0.3e1 / 0.4e1 * (s * (eta + 0.1e1 / 0.3e1) * dphids * (eta - 0.1e1) * L - 0.2e1 * (eta - 0.1e1 / 0.3e1) * c) * cn / L;
B0_4_7 = 0.3e1 / 0.4e1 * (c * dphids * (eta - 1) * (eta + 1) * L + 2 * eta * s) * cn / (L ^ 2);
B0_4_8 = -0.3e1 / 0.4e1 * cn * ((eta + 1) * (eta - 0.1e1 / 0.3e1) * c * dphids * L + 0.2e1 * s * (eta + 0.1e1 / 0.3e1)) / L;
B0_4_9 = 0;
B0_4_10 = 0;
B0_4_11 = -0.3e1 / 0.4e1 * (s * dphids * (eta - 1) * (eta + 1) * L - 2 * c * eta) * cn / (L ^ 2);
B0_4_12 = 0.3e1 / 0.4e1 * (s * (eta + 1) * (eta - 0.1e1 / 0.3e1) * dphids * L - 0.2e1 * (eta + 0.1e1 / 0.3e1) * c) * cn / L;

B0_5_1 = (s * (L * n ^ 2 * eta2 + (n ^ 2 * L - 3 * c * r) * eta - 2 * n ^ 2 * L - 3 * c * r) * (eta - 1) * cn / r ^ 2 / L) / 0.4e1;
B0_5_2 = (s * (eta - 1) * cn * ((-3 * c * eta - c) * r + L * n ^ 2 * (eta - 1) * (eta + 1)) / r ^ 2) / 0.4e1;
B0_5_3 = (s * n * cn * (eta + 2) * (eta - 1) ^ 2 / r ^ 2) / 0.4e1;
B0_5_4 = (s * n * cn * L * (eta + 1) * (eta - 1) ^ 2 / r ^ 2) / 0.4e1;
B0_5_5 = (cn * (L * n ^ 2 * eta2 + (n ^ 2 * L - 3 * c * r) * eta - 2 * n ^ 2 * L - 3 * c * r) * (eta - 1) * c / r ^ 2 / L) / 0.4e1;
B0_5_6 = (cn * ((-3 * eta - 1) * c * r + L * n ^ 2 * (eta - 1) * (eta + 1)) * (eta - 1) * c / r ^ 2) / 0.4e1;
B0_5_7 = -(s * (eta + 1) * (-L * eta * n ^ 2 + L * n ^ 2 * eta2 - 2 * n ^ 2 * L - 3 * c * r * eta + 3 * c * r) * cn / r ^ 2 / L) / 0.4e1;
B0_5_8 = (s * ((-3 * c * eta + c) * r + L * n ^ 2 * (eta - 1) * (eta + 1)) * (eta + 1) * cn / r ^ 2) / 0.4e1;
B0_5_9 = -(s * n * cn * (eta - 2) * (eta + 1) ^ 2 / r ^ 2) / 0.4e1;
B0_5_10 = (s * n * cn * L * (eta - 1) * (eta + 1) ^ 2 / r ^ 2) / 0.4e1;
B0_5_11 = -((n ^ 2 * (eta + 1) * (eta - 2) * L - 3 * r * c * (eta - 1)) * (eta + 1) * c * cn / r ^ 2 / L) / 0.4e1;
B0_5_12 = (cn * (eta + 1) * ((-3 * eta + 1) * c * r + L * n ^ 2 * (eta - 1) * (eta + 1)) * c / r ^ 2) / 0.4e1;

B0_6_1 = (n * sn * (c * (eta + 2) * (eta - 1) * (dphids * r - 2 * s) * L + 6 * r * s * (eta + 1)) * (eta - 1) / r ^ 2 / L) / 0.4e1;
B0_6_2 = (n * ((L * c * dphids * eta2 - c * dphids * L + 6 * eta * s + 2 * s) * r - 2 * L * c * s * (eta - 1) * (eta + 1)) * (eta - 1) * sn / r ^ 2) / 0.4e1;
B0_6_3 = (sn * (c * (eta + 2) * (eta - 1) * (dphids * r - 2 * s) * L + 3 * r * s * (eta + 1)) * (eta - 1) / r ^ 2 / L) / 0.4e1;
B0_6_4 = (sn * ((L * c * dphids * eta2 - c * dphids * L + 3 * eta * s + s) * r - 2 * L * c * s * (eta - 1) * (eta + 1)) * (eta - 1) / r ^ 2) / 0.4e1;
B0_6_5 = -n * sn * (eta - 1) * ((eta + 2) * (eta - 1) * (s * dphids * r / 0.2e1 + c ^ 2) * L - 0.3e1 * r * c * (eta + 1)) / r ^ 2 / L / 0.2e1;
B0_6_6 = -n * sn * ((L * dphids * eta2 * s / 0.2e1 - s * dphids * L / 0.2e1 - (3 * c * eta) - c) * r + L * (c ^ 2) * (eta - 1) * (eta + 1)) * (eta - 1) / r ^ 2 / 0.2e1;
B0_6_7 = -(n * (eta + 1) * sn * (c * (eta + 1) * (eta - 2) * (dphids * r - 2 * s) * L + 6 * r * s * (eta - 1)) / r ^ 2 / L) / 0.4e1;
B0_6_8 = (n * ((L * c * dphids * eta2 - c * dphids * L + 6 * eta * s - 2 * s) * r - 2 * L * c * s * (eta - 1) * (eta + 1)) * (eta + 1) * sn / r ^ 2) / 0.4e1;
B0_6_9 = -((c * (eta + 1) * (eta - 2) * (dphids * r - 2 * s) * L + 3 * r * s * (eta - 1)) * (eta + 1) * sn / r ^ 2 / L) / 0.4e1;
B0_6_10 = (sn * ((L * c * dphids * eta2 - c * dphids * L + 3 * eta * s - s) * r - 2 * L * c * s * (eta - 1) * (eta + 1)) * (eta + 1) / r ^ 2) / 0.4e1;
B0_6_11 = n * (eta + 1) * sn * ((eta - 2) * (eta + 1) * (s * dphids * r / 0.2e1 + c ^ 2) * L - 0.3e1 * r * c * (eta - 1)) / r ^ 2 / L / 0.2e1;
B0_6_12 = -n * ((L * dphids * eta2 * s / 0.2e1 - s * dphids * L / 0.2e1 - (3 * c * eta) + c) * r + L * (c ^ 2) * (eta - 1) * (eta + 1)) * (eta + 1) * sn / r ^ 2 / 0.2e1;


B0(1,1) = B0_1_1;
B0(1,2) = B0_1_2;
B0(1,3) = B0_1_3;
B0(1,4) = B0_1_4;
B0(1,5) = B0_1_5;
B0(1,6) = B0_1_6;
B0(1,7) = B0_1_7;
B0(1,8) = B0_1_8;
B0(1,9) = B0_1_9;
B0(1,10) = B0_1_10;
B0(1,11) = B0_1_11;
B0(1,12) = B0_1_12;

B0(2,1) = B0_2_1;
B0(2,2) = B0_2_2;
B0(2,3) = B0_2_3;
B0(2,4) = B0_2_4;
B0(2,5) = B0_2_5;
B0(2,6) = B0_2_6;
B0(2,7) = B0_2_7;
B0(2,8) = B0_2_8;
B0(2,9) = B0_2_9;
B0(2,10) = B0_2_10;
B0(2,11) = B0_2_11;
B0(2,12) = B0_2_12;

B0(3,1) = B0_3_1;
B0(3,2) = B0_3_2;
B0(3,3) = B0_3_3;
B0(3,4) = B0_3_4;
B0(3,5) = B0_3_5;
B0(3,6) = B0_3_6;
B0(3,7) = B0_3_7;
B0(3,8) = B0_3_8;
B0(3,9) = B0_3_9;
B0(3,10) = B0_3_10;
B0(3,11) = B0_3_11;
B0(3,12) = B0_3_12;

B0(4,1) = B0_4_1;
B0(4,2) = B0_4_2;
B0(4,3) = B0_4_3;
B0(4,4) = B0_4_4;
B0(4,5) = B0_4_5;
B0(4,6) = B0_4_6;
B0(4,7) = B0_4_7;
B0(4,8) = B0_4_8;
B0(4,9) = B0_4_9;
B0(4,10) = B0_4_10;
B0(4,11) = B0_4_11;
B0(4,12) = B0_4_12;

B0(5,1) = B0_5_1;
B0(5,2) = B0_5_2;
B0(5,3) = B0_5_3;
B0(5,4) = B0_5_4;
B0(5,5) = B0_5_5;
B0(5,6) = B0_5_6;
B0(5,7) = B0_5_7;
B0(5,8) = B0_5_8;
B0(5,9) = B0_5_9;
B0(5,10) = B0_5_10;
B0(5,11) = B0_5_11;
B0(5,12) = B0_5_12;

B0(6,1) = B0_6_1;
B0(6,2) = B0_6_2;
B0(6,3) = B0_6_3;
B0(6,4) = B0_6_4;
B0(6,5) = B0_6_5;
B0(6,6) = B0_6_6;
B0(6,7) = B0_6_7;
B0(6,8) = B0_6_8;
B0(6,9) = B0_6_9;
B0(6,10) = B0_6_10;
B0(6,11) = B0_6_11;
B0(6,12) = B0_6_12;


G_1_1 = 0.3e1 / 0.4e1 * s * cos(n * theta) * (-1 + eta2) / L;
G_1_2 = s * cos(n * theta) * (-2 * eta - 1 + 3 * eta2) / 0.4e1;
G_1_3 = 0;
G_1_4 = 0;
G_1_5 = 0.3e1 / 0.4e1 * c * cos(n * theta) * (-1 + eta2) / L;
G_1_6 = c * cos(n * theta) * (-2 * eta - 1 + 3 * eta2) / 0.4e1;
G_1_7 = -0.3e1 / 0.4e1 * s * cos(n * theta) * (-1 + eta2) / L;
G_1_8 = s * cos(n * theta) * (2 * eta - 1 + 3 * eta2) / 0.4e1;
G_1_9 = 0;
G_1_10 = 0;
G_1_11 = -0.3e1 / 0.4e1 * c * cos(n * theta) * (-1 + eta2) / L;
G_1_12 = c * cos(n * theta) * (2 * eta - 1 + 3 * eta2) / 0.4e1;

G_2_1 = 0;
G_2_2 = 0;
G_2_3 = 0.3e1 / 0.4e1 * sin(n * theta) * (-1 + eta2) / L;
G_2_4 = sin(n * theta) * (-2 * eta - 1 + 3 * eta2) / 0.4e1;
G_2_5 = 0;
G_2_6 = 0;
G_2_7 = 0;
G_2_8 = 0;
G_2_9 = -0.3e1 / 0.4e1 * sin(n * theta) * (-1 + eta2) / L;
G_2_10 = sin(n * theta) * (2 * eta - 1 + 3 * eta2) / 0.4e1;
G_2_11 = 0;
G_2_12 = 0;

G_3_1 = 0.3e1 / 0.4e1 * c * cos(n * theta) * (-1 + eta2) / L;
G_3_2 = c * cos(n * theta) * (-2 * eta - 1 + 3 * eta2) / 0.4e1;
G_3_3 = 0;
G_3_4 = 0;
G_3_5 = -0.3e1 / 0.4e1 * s * cos(n * theta) * (-1 + eta2) / L;
G_3_6 = -s * cos(n * theta) * (-2 * eta - 1 + 3 * eta2) / 0.4e1;
G_3_7 = -0.3e1 / 0.4e1 * c * cos(n * theta) * (-1 + eta2) / L;
G_3_8 = c * cos(n * theta) * (2 * eta - 1 + 3 * eta2) / 0.4e1;
G_3_9 = 0;
G_3_10 = 0;
G_3_11 = 0.3e1 / 0.4e1 * s * cos(n * theta) * (-1 + eta2) / L;
G_3_12 = -s * cos(n * theta) * (2 * eta - 1 + 3 * eta2) / 0.4e1;

G_4_1 = -s * n * sin(n * theta) * (eta + 2) * ((eta - 1) ^ 2) / r / 0.4e1;
G_4_2 = -s * n * sin(n * theta) * L * (eta + 1) * ((eta - 1) ^ 2) / r / 0.4e1;
G_4_3 = -sin(n * theta) * (eta + 2) * ((eta - 1) ^ 2) * s / r / 0.4e1;
G_4_4 = -sin(n * theta) * L * (eta + 1) * ((eta - 1) ^ 2) * s / r / 0.4e1;
G_4_5 = -c * n * sin(n * theta) * (eta + 2) * ((eta - 1) ^ 2) / r / 0.4e1;
G_4_6 = -c * n * sin(n * theta) * L * (eta + 1) * ((eta - 1) ^ 2) / r / 0.4e1;
G_4_7 = s * n * sin(n * theta) * (eta - 2) * ((eta + 1) ^ 2) / r / 0.4e1;
G_4_8 = -s * n * sin(n * theta) * L * (eta - 1) * ((eta + 1) ^ 2) / r / 0.4e1;
G_4_9 = sin(n * theta) * (-3 * eta - 2 + eta3) * s / r / 0.4e1;
G_4_10 = -sin(n * theta) * L * (eta - 1) * ((eta + 1) ^ 2) * s / r / 0.4e1;
G_4_11 = c * n * sin(n * theta) * (eta - 2) * ((eta + 1) ^ 2) / r / 0.4e1;
G_4_12 = -c * n * sin(n * theta) * L * (eta - 1) * ((eta + 1) ^ 2) / r / 0.4e1;

G_5_1 = cos(n * theta) * (eta + 2) * ((eta - 1) ^ 2) * (c ^ 2 + s ^ 2) / r / 0.4e1;
G_5_2 = cos(n * theta) * L * (eta + 1) * ((eta - 1) ^ 2) * (c ^ 2 + s ^ 2) / r / 0.4e1;
G_5_3 = n * cos(n * theta) * (-3 * eta + 2 + eta3) / r / 0.4e1;
G_5_4 = n * cos(n * theta) * L * (eta + 1) * ((eta - 1) ^ 2) / r / 0.4e1;
G_5_5 = 0;
G_5_6 = 0;
G_5_7 = -cos(n * theta) * (eta - 2) * ((eta + 1) ^ 2) * (c ^ 2 + s ^ 2) / r / 0.4e1;
G_5_8 = cos(n * theta) * L * (eta - 1) * ((eta + 1) ^ 2) * (c ^ 2 + s ^ 2) / r / 0.4e1;
G_5_9 = -n * cos(n * theta) * (eta - 2) * ((eta + 1) ^ 2) / r / 0.4e1;
G_5_10 = n * cos(n * theta) * L * (eta - 1) * ((eta + 1) ^ 2) / r / 0.4e1;
G_5_11 = 0;
G_5_12 = 0;

G_6_1 = -c * n * sin(n * theta) * (eta + 2) * ((eta - 1) ^ 2) / r / 0.4e1;
G_6_2 = -c * n * sin(n * theta) * L * (eta + 1) * ((eta - 1) ^ 2) / r / 0.4e1;
G_6_3 = -sin(n * theta) * (eta + 2) * ((eta - 1) ^ 2) * c / r / 0.4e1;
G_6_4 = -sin(n * theta) * L * (eta + 1) * ((eta - 1) ^ 2) * c / r / 0.4e1;
G_6_5 = s * n * sin(n * theta) * (eta + 2) * ((eta - 1) ^ 2) / r / 0.4e1;
G_6_6 = s * n * sin(n * theta) * L * (eta + 1) * ((eta - 1) ^ 2) / r / 0.4e1;
G_6_7 = c * n * sin(n * theta) * (eta - 2) * ((eta + 1) ^ 2) / r / 0.4e1;
G_6_8 = -c * n * sin(n * theta) * L * (eta - 1) * ((eta + 1) ^ 2) / r / 0.4e1;
G_6_9 = sin(n * theta) * (-3 * eta - 2 + eta3) * c / r / 0.4e1;
G_6_10 = -sin(n * theta) * L * (eta - 1) * ((eta + 1) ^ 2) * c / r / 0.4e1;
G_6_11 = -s * n * sin(n * theta) * (eta - 2) * ((eta + 1) ^ 2) / r / 0.4e1;
G_6_12 = s * n * sin(n * theta) * L * (eta - 1) * ((eta + 1) ^ 2) / r / 0.4e1;


G(1,1) = G_1_1;
G(1,2) = G_1_2;
G(1,3) = G_1_3;
G(1,4) = G_1_4;
G(1,5) = G_1_5;
G(1,6) = G_1_6;
G(1,7) = G_1_7;
G(1,8) = G_1_8;
G(1,9) = G_1_9;
G(1,10) = G_1_10;
G(1,11) = G_1_11;
G(1,12) = G_1_12;

G(2,1) = G_2_1;
G(2,2) = G_2_2;
G(2,3) = G_2_3;
G(2,4) = G_2_4;
G(2,5) = G_2_5;
G(2,6) = G_2_6;
G(2,7) = G_2_7;
G(2,8) = G_2_8;
G(2,9) = G_2_9;
G(2,10) = G_2_10;
G(2,11) = G_2_11;
G(2,12) = G_2_12;

G(3,1) = G_3_1;
G(3,2) = G_3_2;
G(3,3) = G_3_3;
G(3,4) = G_3_4;
G(3,5) = G_3_5;
G(3,6) = G_3_6;
G(3,7) = G_3_7;
G(3,8) = G_3_8;
G(3,9) = G_3_9;
G(3,10) = G_3_10;
G(3,11) = G_3_11;
G(3,12) = G_3_12;

G(4,1) = G_4_1;
G(4,2) = G_4_2;
G(4,3) = G_4_3;
G(4,4) = G_4_4;
G(4,5) = G_4_5;
G(4,6) = G_4_6;
G(4,7) = G_4_7;
G(4,8) = G_4_8;
G(4,9) = G_4_9;
G(4,10) = G_4_10;
G(4,11) = G_4_11;
G(4,12) = G_4_12;

G(5,1) = G_5_1;
G(5,2) = G_5_2;
G(5,3) = G_5_3;
G(5,4) = G_5_4;
G(5,5) = G_5_5;
G(5,6) = G_5_6;
G(5,7) = G_5_7;
G(5,8) = G_5_8;
G(5,9) = G_5_9;
G(5,10) = G_5_10;
G(5,11) = G_5_11;
G(5,12) = G_5_12;

G(6,1) = G_6_1;
G(6,2) = G_6_2;
G(6,3) = G_6_3;
G(6,4) = G_6_4;
G(6,5) = G_6_5;
G(6,6) = G_6_6;
G(6,7) = G_6_7;
G(6,8) = G_6_8;
G(6,9) = G_6_9;
G(6,10) = G_6_10;
G(6,11) = G_6_11;
G(6,12) = G_6_12;


beta1 = (((((3 * dWdS_1 + 3 * dWdS_2) * eta2 + (-2 * dWdS_1 + 2 * dWdS_2) * eta - dWdS_1 - dWdS_2) * c) - 0.2e1 * s * ((-0.3e1 / 0.2e1 * dUdS_1 - 0.3e1 / 0.2e1 * dUdS_2) * eta2 + (dUdS_1 - dUdS_2) * eta + dUdS_1 / 0.2e1 + dUdS_2 / 0.2e1)) * L + 0.3e1 * (-1 + eta2) * (((W1 - W2) * c) + s * (U1 - U2))) / L / 0.4e1;
beta2 = 0;
beta3 = (((((3 * dUdS_1 + 3 * dUdS_2) * eta2 + (-2 * dUdS_1 + 2 * dUdS_2) * eta - dUdS_1 - dUdS_2) * c) + 0.2e1 * s * ((-0.3e1 / 0.2e1 * dWdS_1 - 0.3e1 / 0.2e1 * dWdS_2) * eta2 + (dWdS_1 - dWdS_2) * eta + dWdS_1 / 0.2e1 + dWdS_2 / 0.2e1)) * L + 0.3e1 * (((U1 - U2) * c) - s * (W1 - W2)) * (-1 + eta2)) / L / 0.4e1;
beta4 = 0;
beta5 = ((((dUdS_1 + dUdS_2) * L + U1 - U2) * eta ^ 3 - L * (dUdS_1 - dUdS_2) * eta ^ 2 + ((-dUdS_1 - dUdS_2) * L - 3 * U1 + 3 * U2) * eta + L * (dUdS_1 - dUdS_2) + 2 * U1 + 2 * U2) * (c ^ 2 + s ^ 2) / r) / 0.4e1;
beta6 = 0;

Omega = zeros(6,6);
Omega(1,1) = beta1;
Omega(1,2) = beta2;
Omega(1,3) = beta3;
Omega(2,4) = beta4;
Omega(2,5) = beta5;
Omega(2,6) = beta6;
Omega(3,1) = beta4;
Omega(3,2) = beta5;
Omega(3,3) = beta6;
Omega(3,4) = beta1;
Omega(3,5) = beta2;
Omega(3,6) = beta3;

BL = Omega*G;
B = B0 + BL;
BDTB = r*B'*DT*B;
end
