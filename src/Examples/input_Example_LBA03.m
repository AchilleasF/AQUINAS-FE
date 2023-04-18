%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example LBA 03: LBAs of thin cylinders under uniform lateral pressure.
% A C1-C1 boundary condition case is to be examined for the buckling
% configuration, with the assumption of a membrane prebuckling stress state.
% An analytical solution based in Donnell's equations is also presented
% for comparison.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% [1] N. Yamaki, Elastic Stability of Circular Cylindrical Shells,
% Elsevier Science, North-Holland, 1984.
% [2] L. H. Sobel, Effects of boundary conditions on the stability of
% cylinders subject to lateral and axial pressures, AIAA Journal, 1964, 2:8, 1437-1440
% https://doi.org/10.2514/3.2572
% [3] T. Vodenitcharova, P. Ansourian, Buckling of circular cylindrical
% shells subject to uniform lateral pressure, Engineering Structures,
% Volume 18, Issue 8, 1996, 604-614
% https://doi.org/10.1016/0141-0296(95)00174-3.
% [4]  M.F. Beatty, F. Pan, On Determinants having Complex Conjugate Columns or Rows
% Journal of Elasticity 47, 69–72 (1997)
% https://doi.org/10.1023/A:1007405421507
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.

%%%%%%%%%%%%%%%%%
% Preliminaries %
%%%%%%%%%%%%%%%%%
clearvars, clc, close all force, format longg
addpath('..','../AQUINAS_Analysis','../AQUINAS_Objects','../Maple_Source');
warning('off','all'); % Supress warnings in order to avoid messages relating to complex eigenvalues


%%%%%%%%%%
% Inputs %
%%%%%%%%%%
% Geometry
r = 1000.0; % [mm] - Cylinder radius
t = 1.0; % [mm] - Cylinder thickness
Zs = [1:0.5:9.5 10:5:95 100:50:950 1000:500:9500 10000:5000:100000]; % [-] - Batdorf parameters to be considered (eqs. 2.7.14 of [1])

% Material
MatName = 'Steel'; % Material name
E = 2.0e5; % [N/mm2] - Young's modulus of steel
nu = 0.3; % [-] - Poisson ratio of steel
D = E*(t^3)/(12*(1-nu^2)); % Bending stiffness


tol = 1e-5; % Tolerance for comparison of determinant with zero
noFzeroStations = 10; % Number of stations, used as initial gueesses for the fzero command to search for a pressure where the determinant of Yamaki's solution
                      % (see function YamakiSolution_CylinderLateralPressure_Donnell below) is close enough to zero, suggesting a solution.

% Data generation option
onlyCriticalN = false; % If onlyCriticalN is true then the Surrogate Optimisation capability of AQUINAS will be employed in order to identify the critical
                       % circumferential wavenumber that the shell will buckle into, as well as the corresponding LPF. If onlyCriticalN is set to false then
                       % there will be no Surrogate Optimisation, and instead the critical circumferential wavenumber will be identified through a brute force
                       % approach of evaluating critical pressures for a range of circumferential wavenumber. The critical wavenumber will then be found as the
                       % one which corresponds to the minimum buckling pressure.
                       % If onlyCriticalN is set to false then an additional 3D surface plot of the critical pressures for the various Batdorf parameters
                       % and circumferential wavenumbers considered is presented. This surface is not shown if the Surrogate Optimisation algorithm is used (onlyCriticalN = false),
                       % since only a selection of circumferential wavenumbers are actually evaluated, for each Batdorf parameter Z.

% Plot controls
FS = 24; % [-] - Font Size (pts)
FSL = 20; % [-] - Font Size for Legends (pts)
MS = 5; % [-] - Markser size (pts)
LW = 3; % [-] - Line Width (pts)

% Initializations
if onlyCriticalN
    AQ_CritLPFs = zeros(length(Zs),1);
    AQ_CritNModes = zeros(length(Zs),1);
    YA_LPFs_allRoots = zeros(length(Zs),noFzeroStations); % Unnecessary array, only for validation of the obtained roots
    YA_CritLPFs = zeros(length(Zs),1);
    YA_CritNModes = zeros(length(Zs),1);
else
    Ns = 2:200; % The range of circumferential wavenumbers to be considered (has to be above 2, as explained in chapter 2.7 of [1])
    AQ_LPFs = zeros(length(Zs),length(Ns));
    AQ_CritLPFs = zeros(length(Zs),1);
    AQ_CritNModes = zeros(length(Zs),1);
    YA_LPFs = zeros(length(Zs),length(Ns));
    YA_LPFs_allRoots = zeros(length(Zs),length(Ns),noFzeroStations); % Unnecessary array, only for validation of the obtained roots
    YA_CritLPFs = zeros(length(Zs),1);
    YA_CritNModes = zeros(length(Zs),1);
end


for I = 1:length(Zs)
    disp(['Computing Z = ',num2str(Zs(I))]);

    if onlyCriticalN

        %%%%%%%%%%%%%%%%%%%%%
        % AQUINAS' solution %
        %%%%%%%%%%%%%%%%%%%%%
        [AQ_CritLPFs(I),AQ_CritNModes(I)] = AQUINAS_CylinderLateralPressure_Donnell_SurrogateOptimisation(Zs(I),r,t,nu,E,MatName);

        %%%%%%%%%%%%%%%%%%%%%
        % Yamaki's solution %
        %%%%%%%%%%%%%%%%%%%%%
        % The following solution is the one presented in chapter 2.7.2 of Yamaki's book [1], and which is based on the Donnell equation for the definition of the initial stability problem.
        f =@(p) YamakiSolution_CylinderLateralPressure_Donnell(Zs(I),AQ_CritNModes(I),p,r,t,nu,D);
        trial_ps = linspace(0,2*AQ_CritLPFs(I),noFzeroStations); % Looking for solutions between zero and up to two times the critical buckling pressure obtained from AQUINAS
        for K = 1:noFzeroStations
            try
                [pzero,fval,exitflag,output] = fzero(f,trial_ps(K));
            catch
                pzero = nan; fval = nan;
            end
            if ~isnan(pzero) && CheckForRootUniqueness(Zs(I),AQ_CritNModes(I),pzero,r,t,nu,D,1e-4,1e-3); pzero = nan; end % Basic check for validity of roots
            if ~isnan(pzero) && pzero<AQ_CritLPFs(I) && CheckForRootUniqueness(Zs(I),AQ_CritNModes(I),pzero,r,t,nu,D,1e-2,1e-3); pzero = nan; end % More relaxed check for validity of roots in case the obtained critical pressure is way lower compared to the one obtained from AQUINAS
            YA_LPFs_allRoots(I,K) = pzero;
        end
        pmin = min(YA_LPFs_allRoots(I,:));
        if pmin <= 0; pmin = nan; end
        YA_CritLPFs(I) = pmin; YA_CritNModes(I) = AQ_CritNModes(I);

    else

        %%%%%%%%%%%%%%%%%%%%%
        % AQUINAS' solution %
        %%%%%%%%%%%%%%%%%%%%%
        LBAs = zeros(1,length(Ns)); J = 0;
        for N = Ns
            J = J + 1;
            LBAs(J) = AQUINAS_CylinderLateralPressure_Donnell_perWavenumber(Zs(I),N,r,t,nu,E,MatName);
        end
        LBAs(LBAs <= 0) = NaN; [m,i] = min(LBAs);
        AQ_LPFs(I,:) = LBAs;
        AQ_CritLPFs(I) = m; AQ_CritNModes(I) = i+1;

        %%%%%%%%%%%%%%%%%%%%%
        % Yamaki's solution %
        %%%%%%%%%%%%%%%%%%%%%
        % The following solution is the one presented in chapter 2.7.2 of Yamaki's book [1], and which is based on the Donnell equation for the definition of the initial stability problem.
        LBAs = zeros(1,length(Ns)); J = 0;
        for N=Ns

            J = J + 1;
            f =@(p) YamakiSolution_CylinderLateralPressure_Donnell(Zs(I),N,p,r,t,nu,D);
            trial_ps = linspace(0,2*AQ_LPFs(I,J),noFzeroStations); % Looking for solutions between zero and up to two times the critical buckling pressure obtained from AQUINAS
            for K = 1:noFzeroStations
                try
                    [pzero,fval,exitflag,output] = fzero(f,trial_ps(K));
                catch
                    pzero = nan; fval = nan;
                end
                if ~isnan(pzero) && CheckForRootUniqueness(Zs(I),N,pzero,r,t,nu,D,1e-4,1e-3); pzero = nan; end % Basic check for validity of roots
                if ~isnan(pzero) && pzero<AQ_LPFs(I,J) && CheckForRootUniqueness(Zs(I),N,pzero,r,t,nu,D,1e-2,1e-3); pzero = nan; end % More relaxed check for validity of roots in case the obtained critical pressure is way lower compared to the one obtained from AQUINAS
                YA_LPFs_allRoots(I,J,K) = pzero;
            end
            pmin = min(YA_LPFs_allRoots(I,J,:));
            LBAs(J) = pmin;
        end
        LBAs(LBAs <= 0) = NaN; [m,i] = min(LBAs);
        YA_LPFs(I,:) = LBAs;
        YA_CritLPFs(I) = m; YA_CritNModes(I) = i+1;

    end

end

% Transformation of results in the kp - beta space for plotting, according to eqs. 2.7.14
AQ_kp = AQ_CritLPFs.*Zs'*(r*t)/sqrt(1-nu^2)*r/(D*(pi^2));
AQ_beta = (sqrt(r*t*Zs'/sqrt(1-nu^2))/pi/r).*AQ_CritNModes;
YA_kp = YA_CritLPFs.*Zs'*(r*t)/sqrt(1-nu^2)*r/(D*(pi^2));
YA_beta = (sqrt(r*t*Zs'/sqrt(1-nu^2))/pi/r).*YA_CritNModes;

pcl = @(Z) 0.92*E*(r./sqrt(r*t*Z/sqrt(1-nu^2)))*((t/r)^2.5);
omegas = @(Z) sqrt(r*t*Z/sqrt(1-nu^2))*sqrt(r/t)/r;

%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%
%% Buckling parameters kp and beta plot (2D)
figure('Name','Buckling parameters kp and beta plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(4,1);

h1 = plot(Zs,YA_kp,'Color',[1.0 0.8 0.8],'LineStyle','-','LineWidth',LW); hold on; grid on;
lgnd{1} = 'Yamaki $k_p$';

h2 = plot(Zs,YA_beta,'Color',[0.8 0.8 1.0],'LineStyle','-','LineWidth',LW);
lgnd{2} = 'Yamaki $\beta$';

h3 = plot(Zs,AQ_kp,'Color','r','LineStyle','None','LineWidth',0.5*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',MS);
lgnd{3} = 'AQUINAS $k_p$ (membrane)';

h4 = plot(Zs,AQ_beta,'Color','b','LineStyle','None','LineWidth',0.5*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',MS);
lgnd{4} = 'AQUINAS $\beta$ (membrane)';

xlim([1 100000]); ylim([0.1 1000]);
legend([h1 h2 h3 h4],lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('$Z$','interpreter','latex','fontsize',FS);
ylabel('$k_p,\beta$','interpreter','latex','fontsize',FS);
annotation('textbox','String','$k_p = prL^2/\pi^2(Et^3/12(1-\nu^2))$ \\ $\beta = (L/\pi r)n$ \\ $Z = \sqrt{1-\nu^2}L^2/rt$','Interpreter','latex','FontSize',FSL,'Units','normalized','Position',[0.2 0.65 0.24 0.15],'BackgroundColor','w','EdgeColor','k');
set(gca,'XScale','log','YScale','log','ticklabelinterpreter','latex','fontsize',FS,'tickdir','out');


%% Buckling pcr/pcl pressure ratio plot
figure('Name','Buckling pcr/pcl pressure ratio plot','NumberTitle','off','WindowState','Maximized','color','w');
lgnd = cell(2,1);

h1 = plot(omegas(Zs),YA_CritLPFs./pcl(Zs)','Color',[0.7 0.7 0.7],'LineStyle','--','LineWidth',LW); hold on; grid on;
lgnd{1} = 'Yamaki';

h2 = plot(omegas(Zs),AQ_CritLPFs./pcl(Zs)','Color','k','LineStyle','-','LineWidth',0.5*LW,'Marker','o','MarkerEdgeColor','k','MarkerSize',MS);
lgnd{2} = 'AQUINAS';

xlim([omegas(Zs(1)) omegas(Zs(end))]);
legend([h1 h2],lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Dimensionless length of cylinder $\omega = \frac{L}{r} \sqrt{\frac{r}{t}}$','interpreter','latex','fontsize',FS);
ylabel('Dimensionless critical buckling pressure  $\frac{p_cr}{p_cl}$','interpreter','latex','fontsize',FS);
set(gca,'XScale','log','YScale','log','ticklabelinterpreter','latex','fontsize',FS,'tickdir','out');

%% 3D surface for pcr/pcl ratios and critical circumferential wavenumber
if ~onlyCriticalN
    % Dimensionless buckling load vs dimensionless length plot (3D)
    figure('Name','Dimensionless buckling load vs dimensionless length plot','NumberTitle','off','WindowState','Maximized','color','w');
    lgnd = cell(3,1);
    AQ_LPFs_norm = zeros(size(AQ_LPFs)); YA_LPFs_norm = zeros(size(YA_LPFs));
    for I = 1:length(Zs)
        AQ_LPFs_norm(I,:) = AQ_LPFs(I,:)/pcl(Zs(I));
        YA_LPFs_norm(I,:) = YA_LPFs(I,:)/pcl(Zs(I));
    end


    [O,N] = meshgrid(omegas(Zs), Ns);
    h1 = surf(O,N,AQ_LPFs_norm','FaceAlpha',0.3,'FaceColor','r'); grid on; hold on;
    lgnd{1} = 'AQUINAS';
    h2 = surf(O,N,YA_LPFs_norm','FaceAlpha',0.3,'FaceColor',[0.3 0.3 1.0]);
    lgnd{2} = 'Yamaki';
    h3 = plot3(omegas(Zs),AQ_CritNModes,AQ_CritLPFs./pcl(Zs)','Color','k','LineStyle','-','LineWidth',LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',MS);
    lgnd{3} = 'AQUINAS critical modes';

    legend([h1 h2 h3],lgnd,'interpreter','latex','fontsize',FSL,'location','best');
    xlabel('Dimensionless length $\Omega = \frac{L}{r}\sqrt{\frac{t}{r}}$','interpreter','latex','fontsize',FS);
    ylabel('Circumferential wave number $n$','interpreter','latex','fontsize',FS);
    zlabel('Critical buckling pressure $p_{cr}$','interpreter','latex','fontsize',FS);
    set(gca,'XScale','log','Zlim',[0,8],'ticklabelinterpreter','latex','fontsize',FS,'tickdir','out');

end

function [detDelta] = YamakiSolution_CylinderLateralPressure_Donnell(Z,N,p,r,t,nu,D)

    L = sqrt(r*t*Z/sqrt(1-nu^2)); % Length of the cylinder corresponding to the current Batdorf parameter (eqs. 2.7.14 of [1])

    % The definition of the following parameters can be found in eqs. 2.7.2 of [1]
    lambda_p = p*(r^3)/D;
    kappa = (t^2)/12/(r^2);
    % The definition of the following parameters can be found in eqs. 2.7.14 of [1]
    gamma = (L/2/r)*N;

    % The factors of the polynomial of eq. 2.7.12 of [1]
    factorPoly_r4 = 1;
    factorPoly_r3 = 4*(N^2);
    factorPoly_r2 = 6*(N^4) + (1-nu^2)/kappa - (N^2)*lambda_p;
    factorPoly_r1 = (N^4)*(4*(N^2)-2*lambda_p);
    factorPoly_r0 = (N^6)*((N^2)-lambda_p);

    polyn = [factorPoly_r4, factorPoly_r3, factorPoly_r2, factorPoly_r1, factorPoly_r0];
    rs = sqrt(roots(polyn));

    % Give priority to the pairs of complex conjugates roots in the arrangement of the roots in rs
    real_roots = []; imag_roots = []; complex_roots = [];
    for I =1:length(rs)
        if abs(real(rs(I)))>0 && abs(imag(rs(I)))>0
            complex_roots(end+1) = rs(I);
        elseif abs(real(rs(I)))>0 && abs(imag(rs(I)))<1e-4
            real_roots(end+1) = rs(I);
        elseif abs(real(rs(I)))<1e-4 && abs(imag(rs(I)))>0
            imag_roots(end+1) = rs(I);
        end
    end
    if ~isempty(complex_roots); rs(1:length(complex_roots)) = cplxpair(complex_roots); end
    if ~isempty(real_roots); rs(length(complex_roots)+1:length(complex_roots)+length(real_roots)) = sort(real_roots,'ComparisonMethod','real'); end
    if ~isempty(imag_roots); rs(length(complex_roots)+length(real_roots)+1:length(complex_roots)+length(real_roots)+length(imag_roots)) = sort(imag_roots,'ComparisonMethod','abs'); end

    rhos = rs*(L/2/r); % transformation of the roots from r to rho, according to eqs 2.7.14 of [1]
    rho1 = rhos(1);
    rho2 = rhos(2);
    rho3 = rhos(3);
    rho4 = rhos(4);

    % aj factors according to eqs. 2.7.16 of [1]
    a1 = rho1*((nu*(rho1^2)-gamma^2)/(((rho1^2)+gamma^2).^2));
    a2 = rho2*((nu*(rho2^2)-gamma^2)/(((rho2^2)+gamma^2).^2));
    a3 = rho3*((nu*(rho3^2)-gamma^2)/(((rho3^2)+gamma^2).^2));
    a4 = rho4*((nu*(rho4^2)-gamma^2)/(((rho4^2)+gamma^2).^2));

    % bj factors according to eqs. 2.7.16 of [1]
    b1 = gamma*(((2+nu)*(rho1^2)+gamma^2)/(((rho1^2)+gamma^2).^2));
    b2 = gamma*(((2+nu)*(rho2^2)+gamma^2)/(((rho2^2)+gamma^2).^2));
    b3 = gamma*(((2+nu)*(rho3^2)+gamma^2)/(((rho3^2)+gamma^2).^2));
    b4 = gamma*(((2+nu)*(rho4^2)+gamma^2)/(((rho4^2)+gamma^2).^2));

    % In the following definitions of the coefficients of the Dj contants, the division with the hyperbolic cosine of the imaginary parts of the roots
    % aids in scaling down the overall magnitude of the obtained determinant and avoid numerical problem (this is a valid column operation for determinants)

    % Coefficients of the Dj constants for the normal displacement W1 in the buckling configuration, according to eqs. 2.7.17 of [1]
    W1D1 =@(x) cos(rho1*2*x/L);
    W1D2 =@(x) cos(rho2*2*x/L);
    W1D3 =@(x) cos(rho3*2*x/L);
    W1D4 =@(x) cos(rho4*2*x/L);

    % Coefficients of the Dj constants for the first derivative of W1 with respect to phi (eqs. 2.7.2 of [1]) in the buckling configuration, based on a derivation of eqs. 2.7.17 of [1]
    dW1D1 =@(x) -rs(1)*sin(rho1*2*x/L);
    dW1D2 =@(x) -rs(2)*sin(rho2*2*x/L);
    dW1D3 =@(x) -rs(3)*sin(rho3*2*x/L);
    dW1D4 =@(x) -rs(4)*sin(rho4*2*x/L);

    % Coefficients of the Dj constants for the tangential displacement U1 in the buckling configuration, according to eqs. 2.7.17 of [1]
    U1D1 =@(x) (L/2/r)*a1*sin(rho1*2*x/L);
    U1D2 =@(x) (L/2/r)*a2*sin(rho2*2*x/L);
    U1D3 =@(x) (L/2/r)*a3*sin(rho3*2*x/L);
    U1D4 =@(x) (L/2/r)*a4*sin(rho4*2*x/L);

    % Coefficients of the Dj constants for the circumferential displacement V1 in the buckling configuration, according to eqs. 2.7.17 of [1]
    V1D1 =@(x) (L/2/r)*b1*cos(rho1*2*x/L);
    V1D2 =@(x) (L/2/r)*b2*cos(rho2*2*x/L);
    V1D3 =@(x) (L/2/r)*b3*cos(rho3*2*x/L);
    V1D4 =@(x) (L/2/r)*b4*cos(rho4*2*x/L);

    % For identical boundary conditions at both ends of the cylinder the determinant reduces to 4x4 instead of the orginal 8x8.
    % The notation used in [3] for the determinant of the buckling problem is used here, where more details about the mathematical
    % nature of this procedure and its roots can be also found.
    DD2 = [ W1D1(L/2)  W1D2(L/2)  W1D3(L/2)  W1D4(L/2)
            dW1D1(L/2)  dW1D2(L/2) dW1D3(L/2) dW1D4(L/2)
            U1D1(L/2)   U1D2(L/2)  U1D3(L/2)  U1D4(L/2)
            V1D1(L/2)   V1D2(L/2)  V1D3(L/2)  V1D4(L/2)];

    % If any of the roots are complex then the determinant that would be obtained at this stage would be complex. Based on the fact that if there are any complex
    % roots then those come in conjugate pairs, the determinant is transformed to a real or imaginary one with the use of column operations, as presented in paragraph 4 of [4].
    if abs(real(rho1))>1e-10 && abs(imag(rho1))>1e-10
        DD2(:,1) = DD2(:,1) + DD2(:,2);
        DD2(:,2) = DD2(:,2) - 0.5 * DD2(:,1);
        DD2(:,2) = DD2(:,2)/sqrt(-1);
    end
    if abs(real(rho3))>1e-10 && abs(imag(rho3))>1e-10
        DD2(:,3) = DD2(:,3) + DD2(:,4);
        DD2(:,4) = DD2(:,4) - 0.5 * DD2(:,3);
        DD2(:,4) = DD2(:,4)/sqrt(-1);
    end

    detDelta = det(DD2);

end


function [duplicateRoots] = CheckForRootUniqueness(Z,N,p,r,t,nu,D,nearZeroTol,comparisonTol)

    duplicateRoots = false;

    L = sqrt(r*t*Z/sqrt(1-nu^2)); % Length of the cylinder corresponding to the current Batdorf parameter (eqs. 2.7.14 of [1])

    % The definition of the following parameters can be found in eqs. 2.7.2 of [1]
    lambda_p = p*(r^3)/D;
    kappa = (t^2)/12/(r^2);


    % The factors of the polynomial of eq. 2.7.12 of [1]
    factorPoly_r4 = 1;
    factorPoly_r3 = 4*(N^2);
    factorPoly_r2 = 6*(N^4) + (1-nu^2)/kappa - (N^2)*lambda_p;
    factorPoly_r1 = (N^4)*(4*(N^2)-2*lambda_p);
    factorPoly_r0 = (N^6)*((N^2)-lambda_p);

    polyn = [factorPoly_r4, factorPoly_r3, factorPoly_r2, factorPoly_r1, factorPoly_r0];
    rs2 = roots(polyn);

    for I=1:length(rs2)
        if abs(real(rs2(I)))<nearZeroTol && abs(imag(rs2(I)))>nearZeroTol; rs2(I) = sqrt(-1)*imag(rs2(I)); end
        if abs(imag(rs2(I)))<nearZeroTol && abs(real(rs2(I)))>nearZeroTol; rs2(I) = real(rs2(I)); end
        if abs(imag(rs2(I)))<nearZeroTol && abs(real(rs2(I)))<nearZeroTol; rs2(I) = 0; end
    end
    rs = sqrt(rs2);
    % Give priority to the pairs of complex conjugates roots in the arrangement of the roots in rs
    real_roots = []; imag_roots = []; complex_roots = [];
    for I =1:length(rs)
        if abs(real(rs(I)))>0 && abs(imag(rs(I)))>0
            complex_roots(end+1) = rs(I);
        elseif abs(real(rs(I)))>0 && abs(imag(rs(I)))<1e-4
            real_roots(end+1) = rs(I);
        elseif abs(real(rs(I)))<1e-4 && abs(imag(rs(I)))>0
            imag_roots(end+1) = rs(I);
        else
            real_roots(end+1) = rs(I); % a zero root goes to the real ones
        end
    end

    for I = 1:length(real_roots)
        for J = I:length(real_roots)
            if I~=J && abs((real_roots(J)-real_roots(I))/real_roots(I))<comparisonTol; duplicateRoots = true; return; end
            if I~=J && abs(real_roots(J)-real_roots(I))<comparisonTol^2; duplicateRoots = true; return; end
        end
    end

    for I = 1:length(imag_roots)
        for J = I+1:length(imag_roots)
            if I~=J && abs((imag_roots(J)-imag_roots(I))/imag_roots(I))<1e-3; duplicateRoots = true; return; end
            if I~=J && abs(imag_roots(J)-imag_roots(I))<1e-6; duplicateRoots = true; return; end
        end
    end

end


function [pcrit] = AQUINAS_CylinderLateralPressure_Donnell_perWavenumber(Z,N,r,t,nu,E,MatName)

    L = sqrt(r*t*Z/sqrt(1-nu^2)); % Length of the cylinder corresponding to the current Batdorf parameter (eqs. 2.7.14 of [1])

    M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)));

    SEG1 = AQUINAS_Segment_Object('type','Cone','material',MatName,'thickness',t,'rzbot',[r 0.0],'rztop',[r L],'els',40);

    % C1 Boundary condition at base of cylinder
    CON1 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'u'},'state','buckling');
    CON2 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'v'},'state','buckling');
    CON3 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'w'});
    CON4 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'beta'},'state','buckling');

    % C1 Boundary condition at top of cylinder
    CON5 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',L,'dofs',{'u'},'state','buckling');
    CON6 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',L,'dofs',{'v'},'state','buckling');
    CON7 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',L,'dofs',{'w'},'state','buckling');
    CON8 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',L,'dofs',{'beta'},'state','buckling');

    P1 = AQUINAS_Distributed_Pressure_Object('segments',{SEG1},'type','pn','functionHandle',@() -1);

    SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',6,'constrLagrMultF',1e-5);

    A1 = AQUINAS_Analysis_Object('type','LBA','circumferentialModes',N,'noEigenvalues',10);

    O1 = AQUINAS_Output_Mode_Object('mode','VarArgOut');

    LBA_C1_out = AQUINAS_Protocol(M1,SEG1,CON1,CON2,CON3,CON4,CON5,CON6,CON7,CON8,P1,SOL1,A1,O1);

    pcrit = LBA_C1_out.CircWave(N).EigenValue{1};

end


function [pcrit,Ncrit] = AQUINAS_CylinderLateralPressure_Donnell_SurrogateOptimisation(Z,r,t,nu,E,MatName)

    L = sqrt(r*t*Z/sqrt(1-nu^2)); % Length of the cylinder corresponding to the current Batdorf parameter (eqs. 2.7.14 of [1])

    M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)));

    SEG1 = AQUINAS_Segment_Object('type','Cone','material',MatName,'thickness',t,'rzbot',[r 0.0],'rztop',[r L],'els',40);

    % C1 Boundary condition at base of cylinder
    CON1 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'u'},'state','buckling');
    CON2 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'v'},'state','buckling');
    CON3 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'w'});
    CON4 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'beta'},'state','buckling');

    % C1 Boundary condition at top of cylinder
    CON5 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',L,'dofs',{'u'},'state','buckling');
    CON6 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',L,'dofs',{'v'},'state','buckling');
    CON7 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',L,'dofs',{'w'},'state','buckling');
    CON8 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',L,'dofs',{'beta'},'state','buckling');

    P1 = AQUINAS_Distributed_Pressure_Object('segments',{SEG1},'type','pn','functionHandle',@() -1);

    SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',6,'surrogate_optimisation',true,'constrLagrMultF',1e-5);

    A1 = AQUINAS_Analysis_Object('type','LBA','noEigenvalues',10);

    O1 = AQUINAS_Output_Mode_Object('mode','VarArgOut');

    LBA_C1_out = AQUINAS_Protocol(M1,SEG1,CON1,CON2,CON3,CON4,CON5,CON6,CON7,CON8,P1,SOL1,A1,O1);

    Nkeys = LBA_C1_out.CircWave.keys;

    pcrit = Inf; Ncrit = nan;
    for I = 1:length(Nkeys)
        if LBA_C1_out.CircWave(Nkeys{I}).EigenValue{1} < pcrit
            pcrit = LBA_C1_out.CircWave(Nkeys{I}).EigenValue{1};
            Ncrit = Nkeys{I};
        end
    end

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