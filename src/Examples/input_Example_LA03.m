%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example LA 03: LA of multi-strake cylindrical silos, with stepwise-varying
% thickness, under nonlinear loading. A BC1r / C1 condition is
% assumed at the base, and a BC2f / S3 condition at the top.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% [1] A. Boyez, A.J. Sadowski, B.A. Izzuddin - "A novel ‘boundary layer’
% finite element for the efficient analysis of thin cylindrical shells"
% Computers & Structures, Volume 182, 2017, Pages 573-587.
% https://doi.org/10.1016/j.compstruc.2016.10.016
% [2] A.J. Sadowski, J.M. Rotter - "Steel silos with different aspect
% ratios: I — Behaviour under concentric discharge"
% Journal of Constructional Steel Research,
% Volume 67, Issue 10, 2011, Pages 1537-1544.
% https://doi.org/10.1016/j.jcsr.2011.03.028
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
% DISCLAIMER : Any alteration to the following input variables (except from those used to define properties for plots)
% will result in poor comparison with the CSBL element from [1]. The CSBL element data for the comparison have already
% been generated for Example 2 of [1] and therefore any change to the following variables will have an impact only on AQUINAS' results.

% Silo names, see Fig. 7 of [1]
names = { 'VS' , 'S' , 'B' , 'I' , 'Q' };

% Geometrical properties for the definition of the silos, given in Table 2 of [1]
HDratio = [ 5.2 , 3 , 2.06 , 1.47 , 0.65 ]; % [-] - H/D : Height over Diameter ratio of the { VS , S , B , I , Q } silos
h = { [ 8800 , 3600 , 4400 , 5600 , 3600 ] ,...
      [ 8200 , 2800 , 3200 , 3800 ] ,...
      [ 8000 , 2400 , 2600 , 1000 ] ,...
      [ 8200 , 2200 , 800 ] ,...
      [ 3300 , 2700 , 500 ] }; % [mm] - Individual heights of the strakes of the { VS , S , B , I , Q } silos
t = { [ 3 , 4 , 5 , 6 , 7 ] ,...
      [ 3 , 4 , 5 , 6 ] ,...
      [ 3 , 4 , 5 , 6 ] ,...
      [ 3 , 4 , 5 ] ,...
      [ 1 , 2 , 3 ] }; % % [mm] - Individual thicknesses of the strakes of the { VS , S , B , I , Q } silos

% Material
MatName = 'Steel'; % Material name
E = 2.0e5; % [N/mm2] - Young's modulus of steel
nu = 0.3; % [-] - Poisson ratio of steel

% The notation of [2] has been adopted for the following variables (see Par. 2 of [2])
h0 = 0; % [mm] height between the equivalent surface of the solid and the highest solid to wall contact
mi = 0.4408; % [-] fully developed wall friction coefficient between the granular solid and the wall
K = 0.5994; % [-] lateral pressure ratio
unit_weight = 9e-6; % [N/mm3] unit weight of wheat (notated as 'gamma' in [2])
phir = 34; % [degrees]

% Plot controls
FS = 20; % [-] - Font Size (pts)
FSL = 18; % [-] - Font Size for Legends (pts)
MS = 3; % [-] - Markser size (pts)
LW = 3; % [-] - Line Width (pts)

% Subplot controls
sFS = 24; % [-] - Font Size (pts)
sFSL = 16; % [-] - Font Size for Legends (pts)
sMS = 3; % [-] - Markser size (pts)
sLW = 3; % [-] - Line Width (pts)

% Plot colors
colors = { 'red', 'green' , 'blue' , 'cyan' , 'magenta' }; % Colors for plotting the { VS , S , B , I , Q } silos

%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%
AQresults = cell(5,1);

for I=1:length(HDratio)
    disp(['Computing silo with height over diameter ratio H/d = ',num2str(HDratio(I))]);
    % Initializations for each silo case
    SEGMENTS = AQUINAS_Segment_Object.empty;
    H = sum(h{I}); %[mm] total height of silo
    R = H/HDratio(I)/2; %[mm] radius of silo

    z0 = R/2/mi/K; % Janssen reference depth, see eqs. (1), (2) of [2]
    ph0 = K*unit_weight*z0; % Wall pressure at infinite depth, see eqs. (1), (2) of [2]
    n = - (1 + tand(phir)) * (1 - h0/z0); % Parameter n, see eq. (2) of [2]
    % For the following categorization depending on the H/D ratio of the silos see Table 1 of [2]
    if HDratio(I) >= 2.0
        ph = @(z) ph0*(1-exp(-(H-z)/z0)); % [N/mm2] Radial Jansen pressure, eq. (1) of [2]
        pw = @(z) mi*ph0*(1-exp(-(H-z)/z0)); % [N/mm2] Tangential Jansen pressure (friction), Par. 2 of [2]
    elseif HDratio(I) > 0.4
        ph = @(z) ph0*(1-(((((H-z)-h0)/(z0-h0))+1)^n)); % [N/mm2] Radial Modified Reimbert pressure, eq. (2) of [2]
        pw = @(z) mi*ph0*(1-(((((H-z)-h0)/(z0-h0))+1)^n)); % [N/mm2] Tangential Modified Reimbert pressure (friction), Par. 2 of [2]
    end

    M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)));

    % Loop for definition of the different strakes of the silo
    hprev = H;
    for S = 1:length(h{I})
        SEGMENTS(end+1) = AQUINAS_Segment_Object('type','Cone','material',MatName,'thickness',t{I}(S),'rzbot',[R hprev-h{I}(S)],'rztop',[R hprev],'els',50,'meshtype','E','g',0.6);
        hprev = hprev - h{I}(S);
    end

    % BC1r condition at base of cylinder
    C1 = AQUINAS_Constraint_Object('rcoord',R,'zcoord',0.0,'dofs',{'u'});
    C2 = AQUINAS_Constraint_Object('rcoord',R,'zcoord',0.0,'dofs',{'w'});
    C3 = AQUINAS_Constraint_Object('rcoord',R,'zcoord',0.0,'dofs',{'b'});

    % BC2f condition at top of cylinder
    C4 = AQUINAS_Constraint_Object('rcoord',R,'zcoord',H,'dofs',{'u'});

    SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',6);

    A1 = AQUINAS_Analysis_Object('type','LA');

    O1 = AQUINAS_Output_Mode_Object('mode','VarArgOut');

    switch length(SEGMENTS)
    case 3
        P1 = AQUINAS_Distributed_Pressure_Object('segments',{SEGMENTS(1) SEGMENTS(2) SEGMENTS(3)},'type','pn','functionHandle',ph,'withRespectTo','z');
        P2 = AQUINAS_Distributed_Pressure_Object('segments',{SEGMENTS(1) SEGMENTS(2) SEGMENTS(3)},'type','pt','functionHandle',pw,'withRespectTo','z');
        AQresults{I} = AQUINAS_Protocol(M1,SEGMENTS(1),SEGMENTS(2),SEGMENTS(3),C1,C2,C3,C4,P1,P2,SOL1,A1,O1);
    case 4
        P1 = AQUINAS_Distributed_Pressure_Object('segments',{SEGMENTS(1) SEGMENTS(2) SEGMENTS(3) SEGMENTS(4)},'type','pn','functionHandle',ph,'withRespectTo','z');
        P2 = AQUINAS_Distributed_Pressure_Object('segments',{SEGMENTS(1) SEGMENTS(2) SEGMENTS(3) SEGMENTS(4)},'type','pt','functionHandle',pw,'withRespectTo','z');
        AQresults{I} = AQUINAS_Protocol(M1,SEGMENTS(1),SEGMENTS(2),SEGMENTS(3),SEGMENTS(4),C1,C2,C3,C4,P1,P2,SOL1,A1,O1);
    case 5
        P1 = AQUINAS_Distributed_Pressure_Object('segments',{SEGMENTS(1) SEGMENTS(2) SEGMENTS(3) SEGMENTS(4) SEGMENTS(5)},'type','pn','functionHandle',ph,'withRespectTo','z');
        P2 = AQUINAS_Distributed_Pressure_Object('segments',{SEGMENTS(1) SEGMENTS(2) SEGMENTS(3) SEGMENTS(4) SEGMENTS(5)},'type','pt','functionHandle',pw,'withRespectTo','z');
        AQresults{I} = AQUINAS_Protocol(M1,SEGMENTS(1),SEGMENTS(2),SEGMENTS(3),SEGMENTS(4),SEGMENTS(5),C1,C2,C3,C4,P1,P2,SOL1,A1,O1);
    end

    clear SEGMENTS M1 P1 P2 C1 C2 C3 C4 C5 C6 SOL1 A1 O1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import data for CSBL element
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading the results of Example 2 of [1] a .csv file
DataCSBL = readtable('input_Example_LA03.csv'); DataCSBL = table2array(DataCSBL);

CSBLind = cell(5,1);
CSBLstations = zeros(5,1);
CSBL_Z = cell(5,1); CSBL_u = cell(5,1); CSBL_w = cell(5,1); CSBL_kzz = cell(5,1); CSBL_kth = cell(5,1);
CSBL_ezz = cell(5,1); CSBL_eth = cell(5,1); CSBL_szz = cell(5,1); CSBL_sth = cell(5,1);
for I=1:5
    if any(~isnan(DataCSBL(:,9*(I-1)+1))); CSBLstations(I) = length(DataCSBL(:,9*(I-1)+1)); else; CSBLstations(I) = find(isnan(DataCSBL(:,9*(I-1)+1)),1) - 1; end
    CSBL_Z{I} = DataCSBL(1:CSBLstations(I),9*(I-1)+1);
    CSBL_u{I} = DataCSBL(1:CSBLstations(I),9*(I-1)+2);
    CSBL_w{I} = DataCSBL(1:CSBLstations(I),9*(I-1)+3);
    CSBL_szz{I} = [DataCSBL(1:CSBLstations(I),9*(I-1)+6), 0.5*(DataCSBL(1:CSBLstations(I),9*(I-1)+6) + DataCSBL(1:CSBLstations(I),9*(I-1)+7)), DataCSBL(1:CSBLstations(I),9*(I-1)+7)];
    CSBL_sth{I} = [DataCSBL(1:CSBLstations(I),9*(I-1)+8), 0.5*(DataCSBL(1:CSBLstations(I),9*(I-1)+8) + DataCSBL(1:CSBLstations(I),9*(I-1)+9)), DataCSBL(1:CSBLstations(I),9*(I-1)+9)];
    CSBL_ezz{I} = (CSBL_szz{I} - nu*CSBL_sth{I})/E;
    CSBL_eth{I} = (CSBL_sth{I} - nu*CSBL_szz{I})/E;
    CSBL_kzz{I} = zeros(CSBLstations(I),1);
    CSBL_kth{I} = zeros(CSBLstations(I),1);
    % Extraction of start-end indices of each segment of the silos, for plotting reasons
    % (in order not to force continuity of inherently discontinuous results where the segments connect, for strains for example)
    height_indices = zeros(length(h{I})+1,1);
    for S=1:length(h{I})
        height_indices(S+1) = find(CSBL_Z{I}==sum(h{I}(length(h{I})-S+1:end)),1);
    end
    CSBLind{I} = height_indices;
    for S=1:length(h{I})
        CSBL_kzz{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1)) = (CSBL_ezz{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),3) - CSBL_ezz{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),1)) / t{I}(length(t{I})+1-S);
        CSBL_kth{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1)) = (CSBL_eth{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),3) - CSBL_eth{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),1)) / t{I}(length(t{I})+1-S);
    end
end


%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%

%% Radial / normal displacement plot
figure('Name','Radial displacement plot','NumberTitle','off','WindowState','Maximized','color','w');
pH = cell(10,1);
lgnd = cell(10,1);

for I=1:5
    % CSBL solution
    pH{I} =  plot(CSBL_u{I}(1:CSBLstations(I)), CSBL_Z{I}(1:CSBLstations(I))./sum(h{I}),...
        'Color',colors{I},'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
    lgnd{I} = strcat('CSBL element solution for',32,names{I},32,'silo');

    % AQUINAS solution
    for S=1:length(h{I})
        pH{5+I} = plot(AQresults{I}.DOFs.u{S}, AQresults{I}.Shell_Geom.z{S}./sum(h{I}),...
            'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',colors{I},'MarkerSize',MS);
        plot(AQresults{I}.DOFs.u{S}, AQresults{I}.Shell_Geom.z{S}./sum(h{I}),...
            'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors{I},'MarkerSize',0.8*MS);
    end
    lgnd{5+I} = strcat('AQUINAS FE solution for',32,names{I},32,'silo');
end

legend([pH{1} pH{2} pH{3} pH{4} pH{5} pH{6} pH{7} pH{8} pH{9} pH{10}],lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Radial displacement of the wall [$mm$]','interpreter','latex','fontsize',FS);
ylabel('Normalized cylindrical axial coordinate z [-]','interpreter','latex','fontsize',FS);
set(gca,'ticklabelinterpreter','latex','tickdir','out','fontsize',FS);


%% Axial / meridional displacement plot
figure('Name','Axial displacement plot','NumberTitle','off','WindowState','Maximized','color','w');
pH = cell(10,1);
lgnd = cell(10,1);

for I=1:5
    % CSBL solution
    pH{I} =  plot(CSBL_w{I}(1:CSBLstations(I)), CSBL_Z{I}(1:CSBLstations(I))./sum(h{I}),...
        'Color',colors{I},'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
    lgnd{I} = strcat('CSBL element solution for',32,names{I},32,'silo');

    % AQUINAS solution
    for S=1:length(h{I})
        pH{5+I} = plot(AQresults{I}.DOFs.w{S}, AQresults{I}.Shell_Geom.z{S}./sum(h{I}),...
            'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',colors{I},'MarkerSize',MS);
        plot(AQresults{I}.DOFs.w{S}, AQresults{I}.Shell_Geom.z{S}./sum(h{I}),...
            'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors{I},'MarkerSize',0.8*MS);
    end
    lgnd{5+I} = strcat('AQUINAS FE solution for',32,names{I},32,'silo');
end

legend([pH{1} pH{2} pH{3} pH{4} pH{5} pH{6} pH{7} pH{8} pH{9} pH{10}],lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('Axial displacement of the wall [mm]','interpreter','latex','fontsize',FS);
ylabel('Normalized cylindrical axial coordinate z [-]','interpreter','latex','fontsize',FS);
set(gca,'ticklabelinterpreter','latex','tickdir','out','fontsize',FS);


%% Axial / meridional midsurface strain plot
figure('Name','Axial midsurface strain plot','NumberTitle','off','WindowState','Maximized','color','w');
pH = cell(2,1);
lgnd = cell(2,1);
lgnd{1} = "CSBL element";
lgnd{2} = "AQUINAS";

for I=1:5
    % CSBL solution
    ax = subplot(1,5,I); xlimMax = 0; xlimMin = 0;
    for S=1:length(h{I})
        pH{1} =  plot(CSBL_ezz{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),2), CSBL_Z{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1))./sum(h{I}),...
            'Color',colors{I},'LineStyle','-','LineWidth',sLW,'Marker','+','MarkerSize',sMS); grid on; hold on;
        if xlimMax < max(CSBL_ezz{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),2)); xlimMax = max(CSBL_ezz{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),2)); end
        if xlimMin > min(CSBL_ezz{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),2)); xlimMin = min(CSBL_ezz{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),2)); end
    end

    % AQUINAS solution
    for S=1:length(h{I})
        pH{2} = plot(AQresults{I}.Strains.Epsilon.Eps_phi.m{S}, AQresults{I}.Shell_Geom.z{S}./sum(h{I}),...
            'Color','k','LineStyle','-','LineWidth',0.6*sLW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',colors{I},'MarkerSize',sMS);
        plot(AQresults{I}.Strains.Epsilon.Eps_phi.m{S}, AQresults{I}.Shell_Geom.z{S}./sum(h{I}),...
            'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors{I},'MarkerSize',0.8*sMS);
        if xlimMax < max(AQresults{I}.Strains.Epsilon.Eps_phi.m{S}); xlimMax = max(AQresults{I}.Strains.Epsilon.Eps_phi.m{S}); end
        if xlimMin > min(AQresults{I}.Strains.Epsilon.Eps_phi.m{S}); xlimMin = min(AQresults{I}.Strains.Epsilon.Eps_phi.m{S}); end
    end

    % horizontal splitters plotting
    for S=2:length(h{I})
        plot([xlimMin xlimMax],[CSBL_Z{I}(CSBLind{I}(S)+1)./sum(h{I}) CSBL_Z{I}(CSBLind{I}(S)+1)./sum(h{I})],'Marker','none','Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',0.2*sLW);
    end

    legend(ax,[pH{1} pH{2}],lgnd,'interpreter','latex','fontsize',sFSL,'location','best');
    xlabel(strcat(names{I},32,'silo'),'interpreter','latex','fontsize',sFS);
    xlim([xlimMin xlimMax]);
    if I==1; ylabel('Normalized cylindrical axial coordinate z [-]','interpreter','latex','fontsize',FS); end
    set(gca,'ticklabelinterpreter','latex','tickdir','out','fontsize',FS);
end
sgtitle('Axial midsurface strain of the wall [-]','interpreter','latex','fontsize',FS);


%% Axial / meridional inner surface strain plot
figure('Name','Axial inner surface strain plot','NumberTitle','off','WindowState','Maximized','color','w');
pH = cell(2,1);
lgnd = cell(2,1);
lgnd{1} = "CSBL element";
lgnd{2} = "AQUINAS";

for I=1:5
    % CSBL solution
    ax = subplot(1,5,I); xlimMax = 0; xlimMin = 0;
    for S=1:length(h{I})
        pH{1} =  plot(CSBL_ezz{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),1), CSBL_Z{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1))./sum(h{I}),...
            'Color',colors{I},'LineStyle','-','LineWidth',sLW,'Marker','+','MarkerSize',sMS); grid on; hold on;
       if xlimMax < max(CSBL_ezz{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),1)); xlimMax = max(CSBL_ezz{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),1)); end
       if xlimMin > min(CSBL_ezz{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),1)); xlimMin = min(CSBL_ezz{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),1)); end
    end

    % AQUINAS solution
    for S=1:length(h{I})
        pH{2} = plot(AQresults{I}.Strains.Epsilon.Eps_phi.i{S}, AQresults{I}.Shell_Geom.z{S}./sum(h{I}),...
            'Color','k','LineStyle','-','LineWidth',0.6*sLW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',colors{I},'MarkerSize',sMS);
        plot(AQresults{I}.Strains.Epsilon.Eps_phi.i{S}, AQresults{I}.Shell_Geom.z{S}./sum(h{I}),...
            'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors{I},'MarkerSize',0.8*sMS);
        if xlimMax < max(AQresults{I}.Strains.Epsilon.Eps_phi.i{S}); xlimMax = max(AQresults{I}.Strains.Epsilon.Eps_phi.i{S}); end
        if xlimMin > min(AQresults{I}.Strains.Epsilon.Eps_phi.i{S}); xlimMin = min(AQresults{I}.Strains.Epsilon.Eps_phi.i{S}); end
    end

    % horizontal splitters plotting
    for S=2:length(h{I})
        plot([xlimMin xlimMax],[CSBL_Z{I}(CSBLind{I}(S)+1)./sum(h{I}) CSBL_Z{I}(CSBLind{I}(S)+1)./sum(h{I})],'Marker','none','Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',0.2*sLW);
    end

    legend(ax,[pH{1} pH{2}],lgnd,'interpreter','latex','fontsize',sFSL,'location','best');
    xlabel(strcat(names{I},32,'silo'),'interpreter','latex','fontsize',sFS);
    xlim([xlimMin xlimMax]);
    if I==1; ylabel('Normalized cylindrical axial coordinate z [-]','interpreter','latex','fontsize',FS); end
    set(gca,'ticklabelinterpreter','latex','tickdir','out','fontsize',FS);
end
sgtitle('Axial inner surface strain of the wall [-]','interpreter','latex','fontsize',FS);


%% Axial / meridional outer surface strain plot
figure('Name','Axial outer surface strain plot','NumberTitle','off','WindowState','Maximized','color','w');
pH = cell(2,1);
lgnd = cell(2,1);
lgnd{1} = "CSBL element";
lgnd{2} = "AQUINAS";

for I=1:5
    % CSBL solution
    ax = subplot(1,5,I); xlimMax = 0; xlimMin = 0;
    for S=1:length(h{I})
        pH{1} =  plot(CSBL_ezz{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),3), CSBL_Z{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1))./sum(h{I}),...
            'Color',colors{I},'LineStyle','-','LineWidth',sLW,'Marker','+','MarkerSize',sMS); grid on; hold on;
       if xlimMax < max(CSBL_ezz{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),3)); xlimMax = max(CSBL_ezz{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),3)); end
       if xlimMin > min(CSBL_ezz{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),3)); xlimMin = min(CSBL_ezz{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),3)); end
    end

    % AQUINAS solution
    for S=1:length(h{I})
        pH{2} = plot(AQresults{I}.Strains.Epsilon.Eps_phi.o{S}, AQresults{I}.Shell_Geom.z{S}./sum(h{I}),...
            'Color','k','LineStyle','-','LineWidth',0.6*sLW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',colors{I},'MarkerSize',sMS);
        plot(AQresults{I}.Strains.Epsilon.Eps_phi.o{S}, AQresults{I}.Shell_Geom.z{S}./sum(h{I}),...
            'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors{I},'MarkerSize',0.8*sMS);
        if xlimMax < max(AQresults{I}.Strains.Epsilon.Eps_phi.o{S}); xlimMax = max(AQresults{I}.Strains.Epsilon.Eps_phi.o{S}); end
        if xlimMin > min(AQresults{I}.Strains.Epsilon.Eps_phi.o{S}); xlimMin = min(AQresults{I}.Strains.Epsilon.Eps_phi.o{S}); end
    end

    % horizontal splitters plotting
    for S=2:length(h{I})
        plot([xlimMin xlimMax],[CSBL_Z{I}(CSBLind{I}(S)+1)./sum(h{I}) CSBL_Z{I}(CSBLind{I}(S)+1)./sum(h{I})],'Marker','none','Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',0.2*sLW);
    end

    legend(ax,[pH{1} pH{2}],lgnd,'interpreter','latex','fontsize',sFSL,'location','best');
    xlabel(strcat(names{I},32,'silo'),'interpreter','latex','fontsize',sFS);
    xlim([xlimMin xlimMax]);
    if I==1; ylabel('Normalized cylindrical axial coordinate z [-]','interpreter','latex','fontsize',FS); end
    set(gca,'ticklabelinterpreter','latex','tickdir','out','fontsize',FS);
end
sgtitle('Axial outer surface strain of the wall [-]','interpreter','latex','fontsize',FS);


%% Circumferential midsurface strain plot
figure('Name','Circumferential midsurface strain plot','NumberTitle','off','WindowState','Maximized','color','w');
pH = cell(2,1);
lgnd = cell(2,1);
lgnd{1} = "CSBL element";
lgnd{2} = "AQUINAS";

for I=1:5
    % CSBL solution
    ax = subplot(1,5,I); xlimMax = 0; xlimMin = 0;
    for S=1:length(h{I})
        pH{1} =  plot(CSBL_eth{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),2), CSBL_Z{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),1)./sum(h{I}),...
            'Color',colors{I},'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
       if xlimMax < max(CSBL_eth{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),2)); xlimMax = max(CSBL_eth{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),2)); end
       if xlimMin > min(CSBL_eth{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),2)); xlimMin = min(CSBL_eth{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),2)); end
    end

    % AQUINAS solution
    for S=1:length(h{I})
        pH{2} = plot(AQresults{I}.Strains.Epsilon.Eps_theta.m{S}, AQresults{I}.Shell_Geom.z{S}./sum(h{I}),...
            'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',colors{I},'MarkerSize',MS);
        plot(AQresults{I}.Strains.Epsilon.Eps_theta.m{S}, AQresults{I}.Shell_Geom.z{S}./sum(h{I}),...
            'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors{I},'MarkerSize',0.8*MS);
        if xlimMax < max(AQresults{I}.Strains.Epsilon.Eps_theta.m{S}); xlimMax = max(AQresults{I}.Strains.Epsilon.Eps_theta.m{S}); end
        if xlimMin > min(AQresults{I}.Strains.Epsilon.Eps_theta.m{S}); xlimMin = min(AQresults{I}.Strains.Epsilon.Eps_theta.m{S}); end
    end

    % horizontal splitters plotting
    for S=2:length(h{I})
        plot([xlimMin xlimMax],[CSBL_Z{I}(CSBLind{I}(S)+1)./sum(h{I}) CSBL_Z{I}(CSBLind{I}(S)+1)./sum(h{I})],'Marker','none','Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',0.2*sLW);
    end

    legend(ax,[pH{1} pH{2}],lgnd,'interpreter','latex','fontsize',sFSL,'location','best');
    xlabel(strcat(names{I},32,'silo'),'interpreter','latex','fontsize',sFS);
    xlim([xlimMin xlimMax]);
    if I==1; ylabel('Normalized cylindrical axial coordinate z [-]','interpreter','latex','fontsize',FS); end
    set(gca,'ticklabelinterpreter','latex','tickdir','out','fontsize',FS);
end
sgtitle('Circumferential midsurface strain of the wall [-]','interpreter','latex','fontsize',FS);


%% Circumferential inner surface strain plot
figure('Name','Circumferential inner surface strain plot','NumberTitle','off','WindowState','Maximized','color','w');
pH = cell(2,1);
lgnd = cell(2,1);
lgnd{1} = "CSBL element";
lgnd{2} = "AQUINAS";

for I=1:5
    % CSBL solution
    ax = subplot(1,5,I); xlimMax = 0; xlimMin = 0;
    for S=1:length(h{I})
        pH{1} =  plot(CSBL_eth{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),1), CSBL_Z{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),1)./sum(h{I}),...
            'Color',colors{I},'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
       if xlimMax < max(CSBL_eth{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),1)); xlimMax = max(CSBL_eth{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),1)); end
       if xlimMin > min(CSBL_eth{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),1)); xlimMin = min(CSBL_eth{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),1)); end
    end

    % AQUINAS solution
    for S=1:length(h{I})
        pH{2} = plot(AQresults{I}.Strains.Epsilon.Eps_theta.i{S}, AQresults{I}.Shell_Geom.z{S}./sum(h{I}),...
            'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',colors{I},'MarkerSize',MS);
        plot(AQresults{I}.Strains.Epsilon.Eps_theta.i{S}, AQresults{I}.Shell_Geom.z{S}./sum(h{I}),...
            'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors{I},'MarkerSize',0.8*MS);
        if xlimMax < max(AQresults{I}.Strains.Epsilon.Eps_theta.i{S}); xlimMax = max(AQresults{I}.Strains.Epsilon.Eps_theta.i{S}); end
        if xlimMin > min(AQresults{I}.Strains.Epsilon.Eps_theta.i{S}); xlimMin = min(AQresults{I}.Strains.Epsilon.Eps_theta.i{S}); end
    end

    % horizontal splitters plotting
    for S=2:length(h{I})
        plot([xlimMin xlimMax],[CSBL_Z{I}(CSBLind{I}(S)+1)./sum(h{I}) CSBL_Z{I}(CSBLind{I}(S)+1)./sum(h{I})],'Marker','none','Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',0.2*sLW);
    end

    legend(ax,[pH{1} pH{2}],lgnd,'interpreter','latex','fontsize',sFSL,'location','best');
    xlabel(strcat(names{I},32,'silo'),'interpreter','latex','fontsize',sFS);
    xlim([xlimMin xlimMax]);
    if I==1; ylabel('Normalized cylindrical axial coordinate z [-]','interpreter','latex','fontsize',FS); end
    set(gca,'ticklabelinterpreter','latex','tickdir','out','fontsize',FS);
end
sgtitle('Circumferential inner surface strain of the wall [$mm$]','interpreter','latex','fontsize',FS);


%% Circumferential outer surface strain plot
figure('Name','Circumferential outer surface strain plot','NumberTitle','off','WindowState','Maximized','color','w');
pH = cell(2,1);
lgnd = cell(2,1);
lgnd{1} = "CSBL element";
lgnd{2} = "AQUINAS";

for I=1:5
    % CSBL solution
    ax = subplot(1,5,I); xlimMax = 0; xlimMin = 0;
    for S=1:length(h{I})
        pH{1} =  plot(CSBL_eth{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),3), CSBL_Z{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),1)./sum(h{I}),...
            'Color',colors{I},'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
       if xlimMax < max(CSBL_eth{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),3)); xlimMax = max(CSBL_eth{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),3)); end
       if xlimMin > min(CSBL_eth{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),3)); xlimMin = min(CSBL_eth{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),3)); end
    end

    % AQUINAS solution
    for S=1:length(h{I})
        pH{2} = plot(AQresults{I}.Strains.Epsilon.Eps_theta.o{S}, AQresults{I}.Shell_Geom.z{S}./sum(h{I}),...
            'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',colors{I},'MarkerSize',MS);
        plot(AQresults{I}.Strains.Epsilon.Eps_theta.o{S}, AQresults{I}.Shell_Geom.z{S}./sum(h{I}),...
            'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors{I},'MarkerSize',0.8*MS);
        if xlimMax < max(AQresults{I}.Strains.Epsilon.Eps_theta.o{S}); xlimMax = max(AQresults{I}.Strains.Epsilon.Eps_theta.o{S}); end
        if xlimMin > min(AQresults{I}.Strains.Epsilon.Eps_theta.o{S}); xlimMin = min(AQresults{I}.Strains.Epsilon.Eps_theta.o{S}); end
    end

    % horizontal splitters plotting
    for S=2:length(h{I})
        plot([xlimMin xlimMax],[CSBL_Z{I}(CSBLind{I}(S)+1)./sum(h{I}) CSBL_Z{I}(CSBLind{I}(S)+1)./sum(h{I})],'Marker','none','Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',0.2*sLW);
    end

    legend(ax,[pH{1} pH{2}],lgnd,'interpreter','latex','fontsize',sFSL,'location','best');
    xlabel(strcat(names{I},32,'silo'),'interpreter','latex','fontsize',sFS);
    xlim([xlimMin xlimMax]);
    if I==1; ylabel('Normalized cylindrical axial coordinate z [-]','interpreter','latex','fontsize',FS); end
    set(gca,'ticklabelinterpreter','latex','tickdir','out','fontsize',FS);
end
sgtitle('Circumferential outer surface strain of the wall [-]','interpreter','latex','fontsize',FS);


%% Axial / meridional curvature plot
figure('Name','Axial curvature plot','NumberTitle','off','WindowState','Maximized','color','w');
pH = cell(2,1);
lgnd = cell(2,1);
lgnd{1} = "CSBL element";
lgnd{2} = "AQUINAS";

for I=1:5
    % CSBL solution
    ax = subplot(1,5,I); xlimMax = 0; xlimMin = 0;
    for S=1:length(h{I})
        pH{1} =  plot(CSBL_kzz{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1)), CSBL_Z{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),1)./sum(h{I}),...
            'Color',colors{I},'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
       if xlimMax < max(CSBL_kzz{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1))); xlimMax = max(CSBL_kzz{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1))); end
       if xlimMin > min(CSBL_kzz{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1))); xlimMin = min(CSBL_kzz{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1))); end
    end

    % AQUINAS solution
    for S=1:length(h{I})
        pH{2} = plot(AQresults{I}.Strains.Kappa.Kappa_phi{S}, AQresults{I}.Shell_Geom.z{S}./sum(h{I}),...
            'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',colors{I},'MarkerSize',MS);
        plot(AQresults{I}.Strains.Kappa.Kappa_phi{S}, AQresults{I}.Shell_Geom.z{S}./sum(h{I}),...
            'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors{I},'MarkerSize',0.8*MS);
        if xlimMax < max(AQresults{I}.Strains.Kappa.Kappa_phi{S}); xlimMax = max(AQresults{I}.Strains.Kappa.Kappa_phi{S}); end
        if xlimMin > min(AQresults{I}.Strains.Kappa.Kappa_phi{S}); xlimMin = min(AQresults{I}.Strains.Kappa.Kappa_phi{S}); end
    end

    % horizontal splitters plotting
    for S=2:length(h{I})
        plot([xlimMin xlimMax],[CSBL_Z{I}(CSBLind{I}(S)+1)./sum(h{I}) CSBL_Z{I}(CSBLind{I}(S)+1)./sum(h{I})],'Marker','none','Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',0.2*sLW);
    end

    legend(ax,[pH{1} pH{2}],lgnd,'interpreter','latex','fontsize',sFSL,'location','best');
    xlabel(strcat(names{I},32,'silo'),'interpreter','latex','fontsize',sFS);
    xlim([xlimMin xlimMax]);
    if I==1; ylabel('Normalized cylindrical axial coordinate z [-]','interpreter','latex','fontsize',FS); end
    set(gca,'ticklabelinterpreter','latex','tickdir','out','fontsize',FS);
end
sgtitle('Axial curvature of the wall [$mm^{-1}$]','interpreter','latex','fontsize',FS);


%% Axial / meridional midsurface stress plot
figure('Name','Axial midsurface stress plot','NumberTitle','off','WindowState','Maximized','color','w');
pH = cell(2,1);
lgnd = cell(2,1);
lgnd{1} = "CSBL element";
lgnd{2} = "AQUINAS";

for I=1:5
    % CSBL solution
    ax = subplot(1,5,I); xlimMax = 0; xlimMin = 0;
    for S=1:length(h{I})
        pH{1} =  plot(CSBL_szz{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),2), CSBL_Z{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),1)./sum(h{I}),...
            'Color',colors{I},'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
       if xlimMax < max(CSBL_szz{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),2)); xlimMax = max(CSBL_szz{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),2)); end
       if xlimMin > min(CSBL_szz{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),2)); xlimMin = min(CSBL_szz{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),2)); end
    end

    % AQUINAS solution
    for S=1:length(h{I})
        pH{2} = plot(AQresults{I}.Stresses.Sig_phi.m{S}, AQresults{I}.Shell_Geom.z{S}./sum(h{I}),...
            'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',colors{I},'MarkerSize',MS);
        plot(AQresults{I}.Stresses.Sig_phi.m{S}, AQresults{I}.Shell_Geom.z{S}./sum(h{I}),...
            'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors{I},'MarkerSize',0.8*MS);
        if xlimMax < max(AQresults{I}.Stresses.Sig_phi.m{S}); xlimMax = max(AQresults{I}.Stresses.Sig_phi.m{S}); end
        if xlimMin > min(AQresults{I}.Stresses.Sig_phi.m{S}); xlimMin = min(AQresults{I}.Stresses.Sig_phi.m{S}); end
    end

    % horizontal splitters plotting
    for S=2:length(h{I})
        plot([xlimMin xlimMax],[CSBL_Z{I}(CSBLind{I}(S)+1)./sum(h{I}) CSBL_Z{I}(CSBLind{I}(S)+1)./sum(h{I})],'Marker','none','Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',0.2*sLW);
    end

    legend(ax,[pH{1} pH{2}],lgnd,'interpreter','latex','fontsize',sFSL,'location','best');
    xlabel(strcat(names{I},32,'silo'),'interpreter','latex','fontsize',sFS);
    xlim([xlimMin xlimMax]);
    if I==1; ylabel('Normalized cylindrical axial coordinate z [-]','interpreter','latex','fontsize',FS); end
    set(gca,'ticklabelinterpreter','latex','tickdir','out','fontsize',FS);
end
sgtitle('Axial midsurface stress of the wall [$N$/$mm^{2}$]','interpreter','latex','fontsize',FS);


%% Circumferential midsurface stress plot
figure('Name','Circumferential midsurface stress plot','NumberTitle','off','WindowState','Maximized','color','w');
pH = cell(2,1);
lgnd = cell(2,1);
lgnd{1} = "CSBL element";
lgnd{2} = "AQUINAS";

for I=1:5
    % CSBL solution
    ax = subplot(1,5,I); xlimMax = 0; xlimMin = 0;
    for S=1:length(h{I})
        pH{1} =  plot(CSBL_sth{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),2), CSBL_Z{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),1)./sum(h{I}),...
            'Color',colors{I},'LineStyle','-','LineWidth',LW,'Marker','+','MarkerSize',MS); grid on; hold on;
       if xlimMax < max(CSBL_sth{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),2)); xlimMax = max(CSBL_sth{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),2)); end
       if xlimMin > min(CSBL_sth{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),2)); xlimMin = min(CSBL_sth{I}(CSBLind{I}(S)+1:CSBLind{I}(S+1),2)); end
    end

    % AQUINAS solution
    for S=1:length(h{I})
        pH{2} = plot(AQresults{I}.Stresses.Sig_theta.m{S}, AQresults{I}.Shell_Geom.z{S}./sum(h{I}),...
            'Color','k','LineStyle','-','LineWidth',0.6*LW,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',colors{I},'MarkerSize',MS);
        plot(AQresults{I}.Stresses.Sig_theta.m{S}, AQresults{I}.Shell_Geom.z{S}./sum(h{I}),...
            'LineStyle','none','Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors{I},'MarkerSize',0.8*MS);
        if xlimMax < max(AQresults{I}.Stresses.Sig_theta.m{S}); xlimMax = max(AQresults{I}.Stresses.Sig_theta.m{S}); end
        if xlimMin > min(AQresults{I}.Stresses.Sig_theta.m{S}); xlimMin = min(AQresults{I}.Stresses.Sig_theta.m{S}); end
    end

    % horizontal splitters plotting
    for S=2:length(h{I})
        plot([xlimMin xlimMax],[CSBL_Z{I}(CSBLind{I}(S)+1)./sum(h{I}) CSBL_Z{I}(CSBLind{I}(S)+1)./sum(h{I})],'Marker','none','Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',0.2*sLW);
    end

    legend(ax,[pH{1} pH{2}],lgnd,'interpreter','latex','fontsize',sFSL,'location','best');
    xlabel(strcat(names{I},32,'silo'),'interpreter','latex','fontsize',sFS);
    xlim([xlimMin xlimMax]);
    if I==1; ylabel('Normalized cylindrical axial coordinate z [-]','interpreter','latex','fontsize',FS); end
    set(gca,'ticklabelinterpreter','latex','tickdir','out','fontsize',FS);
end
sgtitle('Circumferential midsurface stress of the wall [$N$/$mm^{2}$]','interpreter','latex','fontsize',FS);


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
% THE USE OF THE MATLAB LANGUAGE DOES NOT IMPLY ENDORSEMENT BY MATHWORKS.s