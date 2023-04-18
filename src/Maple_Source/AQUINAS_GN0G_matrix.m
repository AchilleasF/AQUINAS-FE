function GN0G = AQUINAS_GN0G_matrix(eta,theta,n,Ns,dat)
% Last generated on Thu Mar 09 11:44:18 2023

eta2 = eta*eta; eta3 = eta2*eta; L = dat.L;
cn = cos(n*theta); sn = sin(n*theta);
N01 = (2-3*eta+eta3)*0.25; N11 = L*(1-eta-eta2+eta3)*0.25; N02 = (2+3*eta-eta3)*0.25; N12 = L*(-1-eta+eta2+eta3)*0.25; % Eq. 2 in Rotter & Teng (Vol. 31)
r = N01*dat.r1 + N02*dat.r2 + N11*dat.drds1 + N12*dat.drds2; % Eqs. 1a & 2 in Rotter & Teng (Vol. 31)
phi = N01*dat.phi1 + N02*dat.phi2 + N11*dat.dphids1 + N12*dat.dphids2; % Cubic interpolation to compute the angle phi at the current Gauss station
dphids = N01*dat.dphids1 + N02*dat.dphids2 + N11*dat.d2phids2_1 + N12*dat.d2phids2_2; % Eqs 1c & 2 in Rotter & Teng (Vol. 31)
c = cos(phi); s = sin(phi);

Nphi = Ns(1);
Ntheta = Ns(2);
Nphitheta = Ns(3);

G = zeros(6,12);


N0 = zeros(6,6);
N0(1,1) = Nphi;
N0(1,4) = Nphitheta;
N0(2,2) = Nphi;
N0(2,5) = Nphitheta;
N0(3,3) = Nphi;
N0(3,5) = Nphitheta;
N0(4,4) = Ntheta;
N0(4,1) = Nphitheta;
N0(5,5) = Ntheta;
N0(5,2) = Nphitheta;
N0(6,6) = Ntheta;
N0(6,3) = Nphitheta;

G_1_1 = 0.3e1 / 0.4e1 * s * cn * (-1 + eta2) / L;
G_1_2 = (s * cn * (-2 * eta - 1 + 3 * eta2)) / 0.4e1;
G_1_3 = 0;
G_1_4 = 0;
G_1_5 = 0.3e1 / 0.4e1 * c * cn * (-1 + eta2) / L;
G_1_6 = (c * cn * (-2 * eta - 1 + 3 * eta2)) / 0.4e1;
G_1_7 = -0.3e1 / 0.4e1 * s * cn * (-1 + eta2) / L;
G_1_8 = (s * cn * (2 * eta - 1 + 3 * eta2)) / 0.4e1;
G_1_9 = 0;
G_1_10 = 0;
G_1_11 = -0.3e1 / 0.4e1 * c * cn * (-1 + eta2) / L;
G_1_12 = (c * cn * (2 * eta - 1 + 3 * eta2)) / 0.4e1;

G_2_1 = 0;
G_2_2 = 0;
G_2_3 = 0.3e1 / 0.4e1 * sn * (-1 + eta2) / L;
G_2_4 = (sn * (-2 * eta - 1 + 3 * eta2)) / 0.4e1;
G_2_5 = 0;
G_2_6 = 0;
G_2_7 = 0;
G_2_8 = 0;
G_2_9 = -0.3e1 / 0.4e1 * sn * (-1 + eta2) / L;
G_2_10 = (sn * (2 * eta - 1 + 3 * eta2)) / 0.4e1;
G_2_11 = 0;
G_2_12 = 0;

G_3_1 = 0.3e1 / 0.4e1 * c * cn * (-1 + eta2) / L;
G_3_2 = (c * cn * (-2 * eta - 1 + 3 * eta2)) / 0.4e1;
G_3_3 = 0;
G_3_4 = 0;
G_3_5 = -0.3e1 / 0.4e1 * s * cn * (-1 + eta2) / L;
G_3_6 = -(s * cn * (-2 * eta - 1 + 3 * eta2)) / 0.4e1;
G_3_7 = -0.3e1 / 0.4e1 * c * cn * (-1 + eta2) / L;
G_3_8 = (c * cn * (2 * eta - 1 + 3 * eta2)) / 0.4e1;
G_3_9 = 0;
G_3_10 = 0;
G_3_11 = 0.3e1 / 0.4e1 * s * cn * (-1 + eta2) / L;
G_3_12 = -(s * cn * (2 * eta - 1 + 3 * eta2)) / 0.4e1;

G_4_1 = -(s * n * sn * (eta + 2) * (eta - 1) ^ 2 / r) / 0.4e1;
G_4_2 = -(s * n * sn * L * (eta + 1) * (eta - 1) ^ 2 / r) / 0.4e1;
G_4_3 = -(sn * (eta + 2) * (eta - 1) ^ 2 * s / r) / 0.4e1;
G_4_4 = -(sn * L * (eta + 1) * (eta - 1) ^ 2 * s / r) / 0.4e1;
G_4_5 = -(c * n * sn * (eta + 2) * (eta - 1) ^ 2 / r) / 0.4e1;
G_4_6 = -(c * n * sn * L * (eta + 1) * (eta - 1) ^ 2 / r) / 0.4e1;
G_4_7 = (s * n * sn * (eta - 2) * (eta + 1) ^ 2 / r) / 0.4e1;
G_4_8 = -(s * n * sn * L * (eta - 1) * (eta + 1) ^ 2 / r) / 0.4e1;
G_4_9 = (sn * (-3 * eta - 2 + eta3) * s / r) / 0.4e1;
G_4_10 = -(sn * L * (eta - 1) * (eta + 1) ^ 2 * s / r) / 0.4e1;
G_4_11 = (c * n * sn * (eta - 2) * (eta + 1) ^ 2 / r) / 0.4e1;
G_4_12 = -(c * n * sn * L * (eta - 1) * (eta + 1) ^ 2 / r) / 0.4e1;

G_5_1 = (cn * (eta + 2) * (eta - 1) ^ 2 * (c ^ 2 + s ^ 2) / r) / 0.4e1;
G_5_2 = (cn * L * (eta + 1) * (eta - 1) ^ 2 * (c ^ 2 + s ^ 2) / r) / 0.4e1;
G_5_3 = (n * cn * (-3 * eta + 2 + eta3) / r) / 0.4e1;
G_5_4 = (n * cn * L * (eta + 1) * (eta - 1) ^ 2 / r) / 0.4e1;
G_5_5 = 0;
G_5_6 = 0;
G_5_7 = -(cn * (eta - 2) * (eta + 1) ^ 2 * (c ^ 2 + s ^ 2) / r) / 0.4e1;
G_5_8 = (cn * L * (eta - 1) * (eta + 1) ^ 2 * (c ^ 2 + s ^ 2) / r) / 0.4e1;
G_5_9 = -(n * cn * (eta - 2) * (eta + 1) ^ 2 / r) / 0.4e1;
G_5_10 = (n * cn * L * (eta - 1) * (eta + 1) ^ 2 / r) / 0.4e1;
G_5_11 = 0;
G_5_12 = 0;

G_6_1 = -(c * n * sn * (eta + 2) * (eta - 1) ^ 2 / r) / 0.4e1;
G_6_2 = -(c * n * sn * L * (eta + 1) * (eta - 1) ^ 2 / r) / 0.4e1;
G_6_3 = -(sn * (eta + 2) * (eta - 1) ^ 2 * c / r) / 0.4e1;
G_6_4 = -(sn * L * (eta + 1) * (eta - 1) ^ 2 * c / r) / 0.4e1;
G_6_5 = (s * n * sn * (eta + 2) * (eta - 1) ^ 2 / r) / 0.4e1;
G_6_6 = (s * n * sn * L * (eta + 1) * (eta - 1) ^ 2 / r) / 0.4e1;
G_6_7 = (c * n * sn * (eta - 2) * (eta + 1) ^ 2 / r) / 0.4e1;
G_6_8 = -(c * n * sn * L * (eta - 1) * (eta + 1) ^ 2 / r) / 0.4e1;
G_6_9 = (sn * (-3 * eta - 2 + eta3) * c / r) / 0.4e1;
G_6_10 = -(sn * L * (eta - 1) * (eta + 1) ^ 2 * c / r) / 0.4e1;
G_6_11 = -(s * n * sn * (eta - 2) * (eta + 1) ^ 2 / r) / 0.4e1;
G_6_12 = (s * n * sn * L * (eta - 1) * (eta + 1) ^ 2 / r) / 0.4e1;


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



GN0G = r*G'*N0*G;

end
