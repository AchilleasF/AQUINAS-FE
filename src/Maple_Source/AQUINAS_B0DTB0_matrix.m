function B0DTB0 = AQUINAS_B0DTB0_matrix(eta,theta,n,dat,DT)
% Last generated on Thu Mar 09 11:44:26 2023
 
eta2 = eta*eta; eta3 = eta2*eta; L = dat.L;
cn = cos(n*theta); sn = sin(n*theta);
N01 = (2-3*eta+eta3)*0.25; N11 = L*(1-eta-eta2+eta3)*0.25; N02 = (2+3*eta-eta3)*0.25; N12 = L*(-1-eta+eta2+eta3)*0.25; % Eq. 2 in Rotter & Teng (Vol. 31)
r = N01*dat.r1 + N02*dat.r2 + N11*dat.drds1 + N12*dat.drds2; % Eqs. 1a & 2 in Rotter & Teng (Vol. 31)
phi = N01*dat.phi1 + N02*dat.phi2 + N11*dat.dphids1 + N12*dat.dphids2; % Cubic interpolation to compute the angle phi at the current Gauss station
dphids = N01*dat.dphids1 + N02*dat.dphids2 + N11*dat.d2phids2_1 + N12*dat.d2phids2_2; % Eqs 1c & 2 in Rotter & Teng (Vol. 31)
c = cos(phi); s = sin(phi);

B0 = zeros(6,12);

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


B0DTB0 = r*B0'*DT*B0;

end