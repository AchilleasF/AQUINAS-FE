function [ GMNAstep, OSTORAGE ] = AQUINAS_GMNA(MATERIALS,SEGMENTS,CONSTRAINTS,NODAL_FORCES,DISTRIBUTED_PRESSURES,DOF_TRACKERS,SOLVER,ANALYSIS,OUTPUT_MODES)
% AQUINAS_GMNA - Geometrically and/or Materially Nonlinear Analysis with AQUINAS
% AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.

if (strcmp(ANALYSIS.type,'GNA') || strcmp(ANALYSIS.type,'GMNA')); include_nonlinear_G = true; else; include_nonlinear_G = false; end
if (strcmp(ANALYSIS.type,'MNA') || strcmp(ANALYSIS.type,'GMNA')); include_nonlinear_M = true; else; include_nonlinear_M = false; end
nelems = 0;
for S = 1:length(SEGMENTS); nelems = nelems + SEGMENTS(S).els; end % Total number of elements used for the shell's meridian
nnodes = nelems + length(SEGMENTS); % Total node count
dof_count = nnodes*4; % Total DOF count, 4 DOFs per node of each segment
uw = setdiff(1:dof_count,2:2:dof_count); % vector of dof indices, excluding circumferential and rotational dofs, for computations relevant to the arc-length constraint, see par. 9.2 of Teng and Rotter (1989a)
freeuwb = setdiff(1:dof_count,2:4:dof_count); % vector of free active dof indices, excluding circumferential and rotational dofs, for checking for loss of stability of the tangent stiffness matrix
% Inclusion of contraints in the total DOF count, due to the Lagrangian multiplier method
for C = 1:length(CONSTRAINTS)
    if strcmp(CONSTRAINTS(C).type,'ecc'); error('AQUINAS Error: Eccentrically connected segments not supported for a GMNA / GNA / MNA'); end
    dof_count = dof_count + 1;
end
% Determination of output level for storage purposes
output_level = 'basic'; for I = 1:length(OUTPUT_MODES); if strcmp(OUTPUT_MODES(I).level,'extensive'); output_level = 'extensive';  break; end; end

% Initializations of auxiliary fields that will be needed during execution of the geometrically and/or materially nonlinear solver
KSTORAGE = zeros(12,12,nelems);
FDISTR = zeros(12,nelems); FINTER = zeros(12,nelems);
OSTORAGE = zeros(12,nelems);
DELTA = zeros(4*nnodes,1); DELTAconv = DELTA; dDelta_a_prev = DELTA;
elDOFs = zeros(12,nelems); elDOFsconv = elDOFs;
SIGMAS = zeros(6,nelems,SOLVER.No_Gauss_Stations); SIGMASconv = SIGMAS;
DTs = zeros(nelems,SOLVER.No_Gauss_Stations,6,6); DTsconv = DTs;
alphas = zeros(nelems,SOLVER.No_Gauss_Stations,SOLVER.No_Simpson_Stations); alphasconv = alphas;
subFlags = zeros(nelems,SOLVER.No_Gauss_Stations,SOLVER.No_Simpson_Stations);
epsilons = zeros(6,nelems,SOLVER.No_Gauss_Stations); epsilonsconv = epsilons;
sigmas = zeros(3,nelems,SOLVER.No_Gauss_Stations,SOLVER.No_Simpson_Stations); sigmasconv = sigmas;
epsilonsZ = zeros(3,nelems,SOLVER.No_Gauss_Stations,SOLVER.No_Simpson_Stations); epsilonsZconv = epsilonsZ;
epns = zeros(nelems,SOLVER.No_Gauss_Stations,SOLVER.No_Simpson_Stations);
sigma_ys = zeros(nelems,SOLVER.No_Gauss_Stations,SOLVER.No_Simpson_Stations);
if include_nonlinear_M
    el_counter = 0;
    for S = 1:length(SEGMENTS)
        sigma_ys((1+el_counter):(el_counter+SEGMENTS(S).els),:,:) = MATERIALS(SEGMENTS(S).material).sy(1); el_counter = el_counter + SEGMENTS(S).els;
    end
    sigma_ysconv = sigma_ys;
end
Kc_siz = 8; Kc_hsiz = 4; % Condensed stiffness matrix full size (8x8) and half size (4x4)
EIGVALUES = []; % Empty array of eigenvalues in order to avoid uninitialised variable error if bifurcation checks are not performed (necessary argument for SOLVER termination function)
GMNAstep = {}; % Cell array to store output per converged step of the nonlinear solver
for Dkey=DOF_TRACKERS.keys
    DOF_TRACKERS(Dkey{1}) = DOF_TRACKERS(Dkey{1}).extract_DOF_LPF_pair(DELTA,0.0);
end

% Evaluate necessary element properties once here to avoid re-computation
el_counter = 0; % Element counter
for S = 1:length(SEGMENTS)
    t = SEGMENTS(S).thickness; % Element thickness
    for E = 1:SEGMENTS(S).els
        el_counter = el_counter + 1;
        nodeB = E+1; nodeT = E; % Store bottom - top node IDs. The element's start node is the top one and its end the bottom one
        % Extract nodal geometric properties of elements end nodes and store them in a single array
        drds2 = SEGMENTS(S).drds(nodeB); dzds2 = SEGMENTS(S).dzds(nodeB); % dr/ds & dz/ds of element bottom node  (Point 2 in Eqs 1a-b in Rotter & Teng (Vol. 31))
        drds1 = SEGMENTS(S).drds(nodeT); dzds1 = SEGMENTS(S).dzds(nodeT); % dr/ds & dz/ds of element top node  (Point 1 in Eqs 1a-b in Rotter & Teng (Vol. 31))
        phi2 = SEGMENTS(S).phi(nodeB); phi1 = SEGMENTS(S).phi(nodeT); % phi angle at element bottom and top nodes (Points 2 and 1 accordingly in eq 1c in Rotter & Teng (Vol. 31))
        dphids2 = SEGMENTS(S).dphids(nodeB); d2phids2_2 = SEGMENTS(S).d2phids2(nodeB); % dphi / ds & d2phi / ds2 of element bottom node (Point 2 in eq 1c in Rotter & Teng (Vol. 31))
        dphids1 = SEGMENTS(S).dphids(nodeT); d2phids2_1 = SEGMENTS(S).d2phids2(nodeT); % dphi / ds & d2phi / ds2 of element top node (Point 1 in eq 1c in Rotter & Teng (Vol. 31))
        L = abs(SEGMENTS(S).s(nodeB)-SEGMENTS(S).s(nodeT))*0.5; % Element arc half-length, computed as the difference between the arc-length coordinates of bottom and top nodes
        % Segment psi factor regarding the importance of transverse shear strains on the overall response of the segment, when a thick formulation is employed
        if strcmp(SEGMENTS(S).formulation,'thin'); psi = 0; elseif strcmp(SEGMENTS(S).formulation,'thick'); psi = (MATERIALS(SEGMENTS(S).material).E/MATERIALS(SEGMENTS(S).material).G)*(t^2)*(6/5/4/((L)^2)); end
        OSTORAGE(:,el_counter) = [drds2 drds1 dzds2 dzds1 phi2 phi1 dphids2 dphids1 d2phids2_2 d2phids2_1 L psi]';
    end
end
Z =@(z) [1 0 0 z 0 0; 0 1 0 0 z 0; 0 0 1 0 0 z]; % Z matrix, according to eq. 12 of Teng and Rotter (1989a)

if SOLVER.visualise_nonlinear
    % Generation of figure window and subplots
    figure('Name','AQUINAS GMNA solver','NumberTitle','off','WindowState','Maximized','color','w');
    axRes = subplot(2,2,1); set(axRes,'YScale','log','Ylim',[SOLVER.zeta_t/100 100*SOLVER.zeta_d]);
    axDef = subplot(2,2,2);
    axDofT = subplot(2,2,3); hold(axDofT,'on');  title(axDofT,'Translational DOF Trackers'); xlabel(axDofT,'Displacement [L]'); ylabel(axDofT,'LPF [-]'); grid on;
    axDofR = subplot(2,2,4); hold(axDofR,'on');  title(axDofR,'Rotational DOF Trackers'); xlabel(axDofR,'Rotation [rad]'); ylabel(axDofR,'LPF [-]'); grid on;
    drawnow
    if include_nonlinear_M
        figure('Name','AQUINAS GMNA solver plasticity progression','NumberTitle','off','WindowState','Maximized','color','w');
        axAlp = subplot(4,1,1);
        axSvM = subplot(4,1,2);
        axSph = subplot(4,1,3);
        axSth = subplot(4,1,4);
        drawnow
    end
end

% Initialization of auxiliary variables for the entire nonlinear procedure
doOneMoreLoadStep = true; % boolean that determines whether to perform another load step or not
i = 0; % load step increment number
Jim1 = 0; % Number of iterations needed to reach convergence of previous load step
SOLVER.LPF = 0; % Load Proportionality Factor
increment_attempts = 0;
termination_cause = "Unexpected termination of analysis";

% main GMNA loop
while doOneMoreLoadStep

    % Necessary updates of fields relevant to the load incrementation procedure
    i = i + 1;
    if SOLVER.visualise_nonlinear; cla(axRes); hold(axRes,'on'); ylabel(axRes,'Change of nodal displacement magnitudes'); grid on; end
    % Necessary initializations of fields relevant to the Newton Raphson - Arc Length loop
    Dijm1 = zeros(nnodes,1); % Vector of nodal displacement magnitudes of previous iteration
    Depnsj = zeros(nelems,SOLVER.No_Gauss_Stations,SOLVER.No_Simpson_Stations); % Accumulated plastic strain of current iteration
    norm_zeta_prev = nan;
    step_converged = false;
    dDELTA_a = zeros(nnodes*4,1); % Vector of accumulated displacement increments for the current load step
    dDELTA_1 = zeros(nnodes*4,1); % Vector of displacement increment for the 1st iteration of the current load step

    % Newton Raphson - Arc Length loop
    for j = 1:SOLVER.Jmax
        % Necessary initialisations per iteration
        seg_counter = 0; % Segment counter
        el_counter = 0; % Element counter
        dof_sum = 0; % Dof counter

        if j == 1
            % Assemble the global tangent stiffness and the global load vector. The stiffness matrices and the distributed load vectors are only computed for the first modified Newtor Raphson or Arc-Length iteration per Load Step attempted
            MATRIX = cell(1,length(SEGMENTS));
            VECTOR = cell(1,length(SEGMENTS)+length(CONSTRAINTS)+length(NODAL_FORCES));

            for S = 1:length(SEGMENTS)

                seg_counter = seg_counter + 1;
                seg_els = SEGMENTS(S).els;
                el_range = (1:SEGMENTS(S).els) + el_counter;
                t = SEGMENTS(S).thickness; % Element thickness

                [MTV,VCV,offsets] = AQUINAS_initialise_offsets(seg_els,dof_sum,Kc_siz,Kc_hsiz);

                if strcmp(SOLVER.compiler,'C++') % Assembly via C++ MEX capability
                    if strcmp(SEGMENTS(S).type,'Cone'); SEG_type = 1; end
                    if strcmp(SEGMENTS(S).type,'Plate'); SEG_type = 2; end
                    if strcmp(SEGMENTS(S).type,'Ellipse'); SEG_type = 3; end

                    cDistrPressType = uint32(ones(max(length(SEGMENTS(S).distributed_pressures),1),1));
                    cDistrPressMagn = zeros(max(length(SEGMENTS(S).distributed_pressures),1),length(SEGMENTS(S).r));
                    for DP=1:length(SEGMENTS(S).distributed_pressures)
                        magnitudes = @(r,z,phi) DISTRIBUTED_PRESSURES(SEGMENTS(S).distributed_pressures(DP)).magnitude_arrayfun(r,z,phi);
                        cDistrPressType(DP) = uint32(DISTRIBUTED_PRESSURES(SEGMENTS(S).distributed_pressures(DP)).type_ToInt);
                        cDistrPressMagn(DP,:) = magnitudes(SEGMENTS(S).r,SEGMENTS(S).z,zeros(1,length(SEGMENTS(S).z)));
                    end
                    [MTV(:,3), VCV(:,2), KSTORAGE(:,:,el_range), FDISTR(:,el_range), DTs(el_range,:,:,:)] = AQUINAS_Cpp_interface_GMNA(...
                         uint32(SEGMENTS(S).els),...
                         SEGMENTS(S).r,...
                         SEGMENTS(S).z,...
                         uint32(length(SEGMENTS(S).distributed_pressures)),...
                         cDistrPressType,...
                         cDistrPressMagn,...
                         uint32(strcmp(SEGMENTS(S).formulation,'thick')),...
                         SEGMENTS(S).thickness,...
                         [MATERIALS(SEGMENTS(S).material).nu MATERIALS(SEGMENTS(S).material).E MATERIALS(SEGMENTS(S).material).G],...
                         uint32(length(MATERIALS(SEGMENTS(S).material).sy)),...
                         MATERIALS(SEGMENTS(S).material).ep,...
                         MATERIALS(SEGMENTS(S).material).sy,...
                         uint32(SOLVER.No_Gauss_Stations),...
                         uint32(SOLVER.No_Simpson_Stations),...
                         uint32(el_counter),...
                         elDOFs,...
                         SIGMAS,...
                         sigmas,...
                         sigma_ys,...
                         epns,...
                         uint32(include_nonlinear_G),...
                         uint32(include_nonlinear_M),...
                         MTV(:,3),...
                         VCV(:,2),...
                         OSTORAGE(:,el_range),...
                         uint32(sort(cell2mat(offsets.values))),...
                         uint32(SOLVER.No_Threads));
                    el_counter = el_range(end);

                else % Assembly within Matlab

                    if ~include_nonlinear_M; DT = MATERIALS(SEGMENTS(S).material).DT(t,[],[]); end
                    for E = 1:SEGMENTS(S).els
                        el_counter = el_counter + 1;
                        r2 = SEGMENTS(S).r(E+1); z2 = SEGMENTS(S).z(E+1); % Global r-z coordinates of element bottom node (Point 2 in Eqs 1a-b in Rotter & Teng (Vol. 31))
                        r1 = SEGMENTS(S).r(E); z1 = SEGMENTS(S).z(E); % Global r-z coordinates of element top node (Point 1 in Eqs 1a-b in Rotter & Teng (Vol. 31))

                        drds1 = OSTORAGE(2, el_counter); drds2 = OSTORAGE(1, el_counter);
                        dzds1 = OSTORAGE(4, el_counter); dzds2 = OSTORAGE(3, el_counter);
                        phi2 = OSTORAGE(5, el_counter); phi1 = OSTORAGE(6, el_counter);
                        dphids2 = OSTORAGE(7, el_counter); dphids1 = OSTORAGE(8, el_counter);
                        d2phids2_2 = OSTORAGE(9, el_counter); d2phids2_1 = OSTORAGE(10, el_counter);
                        L = OSTORAGE(11, el_counter); psi = OSTORAGE(12, el_counter);

                        dat = struct('elDOFs',elDOFs(:,el_counter),'z2',z2,'z1',z1,'r2',r2,'r1',r1,'L',L,'phi2',phi2,'phi1',phi1,'drds2',drds2,'drds1',drds1,'dzds2',dzds2,'dzds1',dzds1,...
                            'dphids2',dphids2,'dphids1',dphids1,'d2phids2_2',d2phids2_2,'d2phids2_1',d2phids2_1,'t',t,'psi',psi);

                        KBARE_T = zeros(12,12); FBARE = zeros(12,1);
                        for DP=SEGMENTS(S).distributed_pressures
                            FBARE = FBARE + DISTRIBUTED_PRESSURES(DP).eta_integral_equivalent_nodal_load_vector(dat);
                        end
                        if include_nonlinear_M
                            for I = 1:length(SOLVER.Gauss_Nodes)
                                DTs(el_counter,I,:,:) = MATERIALS(SEGMENTS(S).material).DT(t,squeeze(sigmas(:,el_counter,I,:)),squeeze(sigma_ys(el_counter,I,:)),squeeze(epns(el_counter,I,:)));
                            end
                        end
                        for I = 1:length(SOLVER.Gauss_Nodes)
                            eta = SOLVER.Gauss_Nodes(I);
                            if include_nonlinear_M; DT = squeeze(DTs(el_counter,I,:,:)); end
                            if include_nonlinear_G
                                KBARE_T = KBARE_T + SOLVER.Gauss_Weights(I) * L * AQUINAS_BDTB_GN0G_matrix(eta,dat,SIGMAS(:,el_counter,I),DT);
                            else
                                KBARE_T = KBARE_T + SOLVER.Gauss_Weights(I) * L * AQUINAS_B0DTB0_matrix(eta,0,0,dat,DT);
                            end
                        end
                        [MTV,VCV] = AQUINAS_condense_element_matrix('tangent',KBARE_T,[],FBARE,phi1,phi2,E,0,offsets,MTV,VCV); % This performs the transformation and condensation operations - Eqs 63
                        KSTORAGE(:,:,el_counter) = KBARE_T;
                        FDISTR(:,el_counter) = FBARE;
                    end

                end

                MATRIX(seg_counter) = {MTV};
                VECTOR(seg_counter) = {VCV};
                dof_sum = dof_sum + (seg_els + 1) * Kc_hsiz;

            end

            % Adding matrix entries for Lagrange multipliers
            for C = 1:length(CONSTRAINTS)
                if strcmp(CONSTRAINTS(C).origin,'user') && (length(CONSTRAINTS(C).local_dofs)==1) && strcmp(CONSTRAINTS(C).local_dofs{1},'v'); continue; end
                if ~(strcmp(CONSTRAINTS(C).state,'general')||strcmp(CONSTRAINTS(C).state,'prebuckling')); continue; end
                for I = 1:length(CONSTRAINTS(C).global_dofs)
                    if (strcmp(CONSTRAINTS(C).origin,'auto') || strcmp(CONSTRAINTS(C).type,'ecc')) && (all(mod(CONSTRAINTS(C).global_dofs{I},4)==2)); continue; end
                    dof_sum = dof_sum + 1; seg_counter = seg_counter + 1;
                    rowM = dof_sum*ones(length(CONSTRAINTS(C).global_dofs{I}),1);
                    colM = CONSTRAINTS(C).global_dofs{I}';
                    valM = CONSTRAINTS(C).coeffLHS{I}';
                    rowF = dof_sum; valF = CONSTRAINTS(C).coeffRHS(I);
                    MATRIX(seg_counter) = {[[rowM; colM] [colM; rowM] [valM; valM]]};
                    VECTOR(seg_counter) = {[rowF valF]};
                end
                if SOLVER.check_axisym_stability && (i==1 && j==1) && (length(CONSTRAINTS(C).global_dofs)==1) && (length(CONSTRAINTS(C).global_dofs{1})==1) && (mod(CONSTRAINTS(C).global_dofs{1},4)~=2); freeuwb = setdiff(freeuwb,CONSTRAINTS(C).global_dofs{1}); end % Get the free u-w-beta dofs in the system
            end

            % Adding vector entries corresponding to Nodal Force Objects
            for N = 1:length(NODAL_FORCES)
                for S = 1:length(NODAL_FORCES(N).Segment_Object_IDs)
                    seg_counter = seg_counter + 1;
                    VECTOR(seg_counter) = {[NODAL_FORCES(N).Segment_Object_IDs{S}{3}, NODAL_FORCES(N).magnitude]};
                end
            end

            MATRIX = cell2mat(MATRIX'); VECTOR = cell2mat(VECTOR');
            MATRIX = sparse(MATRIX(:,1),MATRIX(:,2),MATRIX(:,3),dof_sum,dof_sum);
            VECTOR = accumarray(VECTOR(:,1),VECTOR(:,2));

            % Exclude stiffness entries corresponding to v DOFs (u, w and beta dofs only)
            uwb = sort([1:4:SEGMENTS(end).bot_dofIDs(1) 3:4:SEGMENTS(end).bot_dofIDs(3) 4:4:SEGMENTS(end).bot_dofIDs(4) SEGMENTS(end).bot_dofIDs(4)+1:dof_sum]);
            DELTA_I = zeros(size(VECTOR)); dDELTA_R = zeros(size(VECTOR));
            DELTA_I(uwb) = MATRIX(uwb,uwb)\VECTOR(uwb);
            [elDOFs_I] = AQUINAS_post_process(DELTA_I,VECTOR,KSTORAGE,FDISTR,OSTORAGE,SEGMENTS,MATERIALS,ANALYSIS,SOLVER,[]); delDOFs_R = zeros(size(elDOFs_I));
        else
            % For iterations after the 1st, the nodal load vector is set to be equal to the out of balance forces acting on the nodes
            dDELTA_R = zeros(size(RESIDUAL));
            dDELTA_R(uwb) = MATRIX(uwb,uwb)\RESIDUAL(uwb);
            [delDOFs_R] = AQUINAS_post_process(dDELTA_R,VECTOR,KSTORAGE,(SOLVER.LPF*FDISTR)-FINTER,OSTORAGE,SEGMENTS,MATERIALS,ANALYSIS,SOLVER,[]);
        end

        % Update dof vector and compute load increment for the arc-length method
        [SOLVER,DELTA,dDELTA_a,dDELTA_1,Dij,zeta] = SOLVER.update_nonlinear_solution(i,j,nnodes,uw,Dijm1,Jim1,DELTA,DELTA_I,dDELTA_R,dDELTA_a,dDELTA_1,dDelta_a_prev);
        if SOLVER.flag < 0; break; end

        if strcmp(SOLVER.nonlinear_solver,'NewtonRaphson')
            if j==1; elDOFs = elDOFs + elDOFs_I*SOLVER.dLPF; else; elDOFs = elDOFs + delDOFs_R; end
        elseif strcmp(SOLVER.nonlinear_solver,'ArcLength')
            elDOFs = elDOFs + elDOFs_I*SOLVER.dLPF + delDOFs_R;
        end

        % Assemble the global tangent stiffness and residual load vector, without recomputing the element stiffness matrices or distributed load vectors
        MATRIX = cell(1,length(SEGMENTS));
        RESIDUAL = cell(1,length(SEGMENTS)+length(CONSTRAINTS)+length(NODAL_FORCES));
        seg_counter = 0; % Segment counter
        el_counter = 0; % Element counter
        dof_sum = 0; % Dof counter

        for S = 1:length(SEGMENTS)

            seg_counter = seg_counter + 1;
            seg_els = SEGMENTS(S).els;
            el_range = (1:SEGMENTS(S).els) + el_counter;

            [MTV,VCV,offsets] = AQUINAS_initialise_offsets(seg_els,dof_sum,Kc_siz,Kc_hsiz);

            if strcmp(SOLVER.compiler,'C++') % Internal forces calculation via C++ MEX capability

                if strcmp(SEGMENTS(S).type,'Cone'); SEG_type = 1; end
                if strcmp(SEGMENTS(S).type,'Plate'); SEG_type = 2; end
                if strcmp(SEGMENTS(S).type,'Ellipse'); SEG_type = 3; end
                [VCV(:,2), MTV(:,3), FINTER(:,el_range), alphas(el_range,:,:), SIGMAS(:,el_range,:), sigmas(:,el_range,:,:), sigma_ys(el_range,:,:), epsilons(:,el_range,:), Depnsj(el_range,:,:), subFlags(el_range,:,:)] = AQUINAS_Cpp_interface_internal_forces(...
                     uint32(SEGMENTS(S).els),...
                     SEGMENTS(S).r,...
                     SEGMENTS(S).z,...
                     SEGMENTS(S).thickness,...
                     [MATERIALS(SEGMENTS(S).material).nu MATERIALS(SEGMENTS(S).material).E MATERIALS(SEGMENTS(S).material).G],...
                     uint32(length(MATERIALS(SEGMENTS(S).material).sy)),...
                     MATERIALS(SEGMENTS(S).material).ep,...
                     MATERIALS(SEGMENTS(S).material).sy,...
                     uint32(SOLVER.No_Gauss_Stations),...
                     uint32(SOLVER.No_Simpson_Stations),...
                     uint32(el_counter),...
                     elDOFs,...
                     uint32(alphas(el_range,:,:)),...
                     sigmas(:,el_range,:,:),...
                     sigma_ys(el_range,:,:),...
                     epsilons(:,el_range,:,:),...
                     Depnsj(el_range,:,:),...
                     epns(el_range,:,:),...
                     [SOLVER.epsilon_s SOLVER.max_epsilon_bar],...
                     uint32(include_nonlinear_G),...
                     uint32(include_nonlinear_M),...
                     VCV(:,2),...
                     MTV(:,3),...
                     SOLVER.LPF*FDISTR(:,el_range),...
                     KSTORAGE(:,:,el_range),...
                     OSTORAGE(:,el_range),...
                     uint32(sort(cell2mat(offsets.values))),...
                     uint32(SOLVER.No_Threads));
                 el_counter = el_range(end);
                 if any(any(any(subFlags))); SOLVER.flag = -6; end % Too large effective strain increment at at least one material point

            else % Assembly within Matlab

                % Going through all of the elements in order to evaluate the internal forces that have developed, see eqs. 51 and 52 in Rotter and Teng(Vol. 31)
                if ~include_nonlinear_M; DT = MATERIALS(SEGMENTS(S).material).DT(t,[],[]); end
                for E = 1:SEGMENTS(S).els
                    el_counter = el_counter + 1;
                    r1 = SEGMENTS(S).r(E); z1 = SEGMENTS(S).z(E); % Global r-z coordinates of element top node (Point 1 in Eqs 1a-b in Rotter & Teng (Vol. 31))
                    r2 = SEGMENTS(S).r(E+1); z2 = SEGMENTS(S).z(E+1); % Global r-z coordinates of element bottom node (Point 2 in Eqs 1a-b in Rotter & Teng (Vol. 31))
                    drds1 = OSTORAGE(2, el_counter); drds2 = OSTORAGE(1, el_counter);
                    dzds1 = OSTORAGE(4, el_counter); dzds2 = OSTORAGE(3, el_counter);
                    phi1 = OSTORAGE(6, el_counter); phi2 = OSTORAGE(5, el_counter);
                    dphids1 = OSTORAGE(8, el_counter); dphids2 = OSTORAGE(7, el_counter);
                    d2phids2_1 = OSTORAGE(10, el_counter); d2phids2_2 = OSTORAGE(9, el_counter);
                    L = OSTORAGE(11, el_counter); psi = OSTORAGE(12, el_counter);

                    dat = struct('elDOFs',elDOFs(:,el_counter),'z2',z2,'z1',z1,'r2',r2,'r1',r1,'L',L,'phi2',phi2,'phi1',phi1,'drds2',drds2,'drds1',drds1,'dzds2',dzds2,'dzds1',dzds1,...
                                 'dphids2',dphids2,'dphids1',dphids1,'d2phids2_2',d2phids2_2,'d2phids2_1',d2phids2_1,'t',t,'psi',psi);

                    FINT = zeros(12,1);
                    for I = 1:length(SOLVER.Gauss_Nodes)
                        eta = SOLVER.Gauss_Nodes(I);
                        if include_nonlinear_M
                            % Rotter and Teng sub-incremental method for the evaluation of the nonlinear stresses
                            epsilons_prev = epsilons(:,el_counter,I);
                            epsilons(:,el_counter,I) = AQUINAS_strains_vector(eta,dat,include_nonlinear_G);
                            Delta_epsilon = epsilons(:,el_counter,I) - epsilons_prev;
                            % Simpson stations
                            for J = 1:SOLVER.No_Simpson_Stations
                                zJ = -t/2 + (J-1)*t/(SOLVER.No_Simpson_Stations-1);
                                epsilon0Z = Z(zJ) * epsilons_prev;
                                DeltaepsilonZ = Z(zJ) * Delta_epsilon;
                                [SOLVER.flag,sigmas(:,el_counter,I,J),Depnsj(el_counter,I,J),alphas(el_counter,I,J)] = MATERIALS(SEGMENTS(S).material).sub_incremental_computations(SOLVER.flag,epsilon0Z,DeltaepsilonZ,sigmas(:,el_counter,I,J),sigma_ys(el_counter,I,J),epns(el_counter,I,J),Depnsj(el_counter,I,J),alphas(el_counter,I,J),SOLVER.epsilon_s,SOLVER.max_epsilon_bar);
                                epsilonsZ(:,el_counter,I,J) = epsilon0Z + DeltaepsilonZ;
                            end
                            SIGMAintegral = zeros(6,1);
                             % Numerical integration according to Simpson's 1/3 rule
                            for J = 1:((SOLVER.No_Simpson_Stations-1)/2)
                                a = 2*J-1; m = 2*J; b = 2*J+1;
                                z_a = -t/2 + (a-1)*t/(SOLVER.No_Simpson_Stations-1); z_m = -t/2 + (m-1)*t/(SOLVER.No_Simpson_Stations-1); z_b = -t/2 + (b-1)*t/(SOLVER.No_Simpson_Stations-1);
                                SIGMAintegral = SIGMAintegral + (z_b-z_a)*(sigmas(:,el_counter,I,a)'*Z(z_a) + 4*sigmas(:,el_counter,I,m)'*Z(z_m) + sigmas(:,el_counter,I,b)'*Z(z_b))'/6;
                            end
                            SIGMAS(:,el_counter,I) = SIGMAintegral;
                        else
                            epsilons(:,el_counter,I) = AQUINAS_strains_vector(eta,dat,include_nonlinear_G);
                            SIGMAS(:,el_counter,I) = DT * epsilons(:,el_counter,I);
                        end
                        FINT = FINT + L * SOLVER.Gauss_Weights(I) * AQUINAS_Fint_vector(eta,dat,SIGMAS(:,el_counter,I),include_nonlinear_G);
                    end
                    FINTER(:,el_counter) = FINT;
                    [MTV,VCV] = AQUINAS_condense_element_matrix('tangent',KSTORAGE(:,:,el_counter),[],(SOLVER.LPF*FDISTR(:,el_counter))-FINT,phi1,phi2,E,0,offsets,MTV,VCV); % This performs the transformation and condensation operations - Eqs 63
                end

            end
            MATRIX(seg_counter) = {MTV};
            RESIDUAL(seg_counter) = {VCV};
            dof_sum = dof_sum + (seg_els + 1) * Kc_hsiz;
        end

        % Adding matrix entries for Lagrange multipliers
        for C = 1:length(CONSTRAINTS)
            if strcmp(CONSTRAINTS(C).origin,'user') && (length(CONSTRAINTS(C).local_dofs)==1) && strcmp(CONSTRAINTS(C).local_dofs{1},'v'); continue; end
            if ~(strcmp(CONSTRAINTS(C).state,'general')||strcmp(CONSTRAINTS(C).state,'prebuckling')); continue; end
            for I = 1:length(CONSTRAINTS(C).global_dofs)
                if (strcmp(CONSTRAINTS(C).origin,'auto') || strcmp(CONSTRAINTS(C).type,'ecc')) && (all(mod(CONSTRAINTS(C).global_dofs{I},4)==2)); continue; end
                dof_sum = dof_sum + 1; seg_counter = seg_counter + 1;
                rowM = dof_sum*ones(length(CONSTRAINTS(C).global_dofs{I}),1);
                colM = CONSTRAINTS(C).global_dofs{I}';
                valM = CONSTRAINTS(C).coeffLHS{I}';
                rowF = dof_sum; valF = CONSTRAINTS(C).coeffRHS(I);
                MATRIX(seg_counter) = {[[rowM; colM] [colM; rowM] [valM; valM]]};
                RESIDUAL(seg_counter) = {[rowF valF]};
            end
        end

        % Adding vector entries corresponding to Nodal Force Objects
        for N = 1:length(NODAL_FORCES)
            for S = 1:length(NODAL_FORCES(N).Segment_Object_IDs)
                seg_counter = seg_counter + 1;
                RESIDUAL(seg_counter) = {[NODAL_FORCES(N).Segment_Object_IDs{S}{3}, SOLVER.LPF*NODAL_FORCES(N).magnitude]};
            end
        end
        MATRIX = cell2mat(MATRIX'); RESIDUAL = cell2mat(RESIDUAL');
        MATRIX = sparse(MATRIX(:,1),MATRIX(:,2),MATRIX(:,3),dof_sum,dof_sum);
        RESIDUAL = accumarray(RESIDUAL(:,1),RESIDUAL(:,2));
        if SOLVER.flag < 0; break; end

        if SOLVER.visualise_nonlinear
            % Deformed shape plot of current iteration
            AQUINAS_Output_Mode_Object.nonlinear_visualiser(axDef,'deformedShape',i,j,DELTA,[],[],alphas,[],SOLVER,SEGMENTS,[],[],include_nonlinear_M);
            % Convergence criterion development
            if j > 1 && ~isempty(zeta); AQUINAS_Output_Mode_Object.nonlinear_visualiser(axRes,'convergenceDevelopment',i,j,[],norm_zeta_prev,zeta,[],[],SOLVER,[],[],[],[]); end
            if include_nonlinear_M
                AQUINAS_Output_Mode_Object.nonlinear_visualiser(axAlp,'alphas',i,j,[],[],[],alphas,[],SOLVER,[],[],OSTORAGE,[]);
                AQUINAS_Output_Mode_Object.nonlinear_visualiser(axSvM,'sigmaVonMises',i,j,[],[],[],[],sigmas,SOLVER,[],[],OSTORAGE,[]);
                AQUINAS_Output_Mode_Object.nonlinear_visualiser(axSph,'sigmaphis',i,j,[],[],[],[],sigmas,SOLVER,[],[],OSTORAGE,[]);
                AQUINAS_Output_Mode_Object.nonlinear_visualiser(axSth,'sigmathetas',i,j,[],[],[],[],sigmas,SOLVER,[],[],OSTORAGE,[]);
            end
        end

        % Check for loss of stability at the axisymmetric prebuckling path
        if SOLVER.check_axisym_stability && j==1
            for sigma_shift = [0,1e-4,1e-3,1e-2]
                try
                    SOLVER.axisym_path_eig = eigs(MATRIX(freeuwb,freeuwb),1,sigma_shift);
                    break
                catch
                    SOLVER.axisym_path_eig = nan;
                    continue
                end
            end
        end

        % Convergence check
        if (norm(zeta,Inf) <= SOLVER.zeta_t) && ~any(isnan(RESIDUAL))
            step_converged = true;
            % Bifurcation checks for current steps
            if any(ismember(SOLVER.termination_conditions,'C')) && (SOLVER.LPF >= SOLVER.check_bifurcation_after_LPF)
                prevEIGVALUES = EIGVALUES;
                if SOLVER.surrogate_optimisation
                    [ EIGVECTORS, EIGVALUES, ANALYSIS ] = AQUINAS_SO(SEGMENTS,CONSTRAINTS,SOLVER,MATERIALS,ANALYSIS,DTs,elDOFs,SIGMAS,OSTORAGE);
                else
                    [ EIGVECTORS, EIGVALUES ] = AQUINAS_LBA(SEGMENTS,CONSTRAINTS,SOLVER,MATERIALS,ANALYSIS,DTs,elDOFs,SIGMAS,OSTORAGE);
                end
                % Ignore the sign of the eigenvalues if the axisymmetric mode is the critical one.
                % There may be negative eigenvalues due to the global tangent stiffness matrix being almost singular towards a limit point of the fundamental path.
                if any(ANALYSIS.circumferential_modes==0) && min(abs(EIGVALUES))==abs(EIGVALUES(1,1)); EIGVALUES = abs(EIGVALUES); end
                % Check for eigenvalues close to unity
                if min(EIGVALUES) < (1 + SOLVER.ksi)
                    if min(EIGVALUES) > (1 - SOLVER.ksi) && sum(EIGVALUES < (1 + SOLVER.ksi)) == 1
                        % A single circumferential mode produces an eigenvalue close to unity, which means buckling into this specific mode has occurred
                        C = find(EIGVALUES < (1 + SOLVER.ksi));
                        GMNAstep{i}.BuckledIntoMode = ANALYSIS.circumferential_modes(C);
                    elseif strcmp(SOLVER.simultaneous_bifurcation_treatment,'pickSmallestWithinRange') && min(EIGVALUES) > (1 - SOLVER.ksi)
                        % More than one circumferential modes produce an eigenvalue close to unity. Buckling into a specific mode has not been strictly detected, but the mode with the smallest eigenvalue is picked as the candidate with the highest prospects for buckling.
                        C = find(EIGVALUES == min(EIGVALUES));
                        GMNAstep{i}.BuckledIntoMode = ANALYSIS.circumferential_modes(C);
                    elseif min(EIGVALUES) < (1 - SOLVER.ksi) || (strcmp(SOLVER.simultaneous_bifurcation_treatment,'divideKsiBy10') && min(EIGVALUES) > (1 - SOLVER.ksi))
                        if strcmp(SOLVER.simultaneous_bifurcation_treatment,'divideKsiBy10') && min(EIGVALUES) > (1 - SOLVER.ksi); SOLVER.ksi = SOLVER.ksi/10; end
                        % Undo the current increment and redo the load step with an appropriate LPF to investigate the possibility of buckling
                        step_converged = false;
                        DELTA = DELTAconv; elDOFs = elDOFsconv; SIGMAS = SIGMASconv;
                        if include_nonlinear_M; DTs = DTsconv; sigmas = sigmasconv; sigma_ys = sigma_ysconv; alphas = alphasconv; epsilons = epsilonsconv; end
                        if i==1
                            error('AQUINAS Error: The loading assigned causes buckling from the first step of the geometrically nonlinear analysis. Consider reducing the applied load or dLPF step factor, or perhaps perform an LBA instead.')
                        else
                            [~,minind] = min(EIGVALUES);
                            % Determine an appropriate load increment factor based on the magnitude of the critical eigenvalue. For analyses where the Gaussian Process Regression has been employed, the modes may be different between steps and therefore the minimum eigenvalue will be considered.
                            if SOLVER.surrogate_optimisation
                                if ~isempty(prevEIGVALUES); SOLVER.dLPF = (SOLVER.LPF - GMNAstep{i-1}.LPF)*(min(prevEIGVALUES) - 1)/(min(prevEIGVALUES)-EIGVALUES(minind))/2; else; SOLVER.dLPF = (SOLVER.LPF - GMNAstep{i-1}.LPF)/(2 + increment_attempts); end
                            else
                                if ~isempty(prevEIGVALUES); SOLVER.dLPF = (SOLVER.LPF - GMNAstep{i-1}.LPF)*(prevEIGVALUES(minind) - 1)/(prevEIGVALUES(minind)-EIGVALUES(minind))/2; else; SOLVER.dLPF = (SOLVER.LPF - GMNAstep{i-1}.LPF)/(2 + increment_attempts); end
                            end
                            if strcmp(SOLVER.nonlinear_solver,'NewtonRaphson')
                                if SOLVER.dLPF > 0; SOLVER.dLPF = sign(dDelta_a_prev'*DELTA_I(1:4*nnodes))*SOLVER.dLPF; else; SOLVER.dLPF = SOLVER.dLPF*SOLVER.reattempt_cutback; end
                            elseif strcmp(SOLVER.nonlinear_solver,'ArcLength')
                                if SOLVER.dLPF > 0; SOLVER.l = SOLVER.dLPF*sqrt(DELTA_I(uw)'*DELTA_I(uw)); else; SOLVER.l = SOLVER.l*SOLVER.reattempt_cutback; end
                            end
                        end
                        SOLVER.LPF = GMNAstep{i-1}.LPF;
                        increment_attempts = increment_attempts + 1;
                        i = i - 1; Jim1 = nan;
                        break
                    end
                end
                GMNAstep{i}.Bifurcation = containers.Map('KeyType','int32','ValueType','any');
                for C = 1:length(ANALYSIS.circumferential_modes)
                    GMNAstep{i}.Bifurcation(ANALYSIS.circumferential_modes(C)) = struct('EigenValue',EIGVALUES(C),'EigenVector',EIGVECTORS(:,:,C));
                end
                if ~isfield(GMNAstep{i},'BuckledIntoMode'); GMNAstep{i}.BuckledIntoMode = []; end
            else
                GMNAstep{i}.BuckledIntoMode = [];
                GMNAstep{i}.Bifurcation = containers.Map('KeyType','int32','ValueType','any');
            end
            % Storing non-bifurcation step results in cell array for output
            if strcmp(SOLVER.nonlinear_solver,'ArcLength'); GMNAstep{i}.l = SOLVER.l; end
            GMNAstep{i}.LPF = SOLVER.LPF;
            if i == 1; GMNAstep{i}.dLPF = SOLVER.LPF; else; GMNAstep{i}.dLPF = SOLVER.LPF - GMNAstep{i-1}.LPF; end
            GMNAstep{i}.Delta = DELTA;
            if strcmp(output_level,'extensive')
                GMNAstep{i}.elDOFs = elDOFs;
                GMNAstep{i}.alphas = alphas;
                GMNAstep{i}.epsilons = epsilons;
                GMNAstep{i}.sigmas = sigmas;
                GMNAstep{i}.SIGMAS = SIGMAS;
            end

            % Update equivalent plastic strains and yield stress of Simpson's stations. There cannot be plastic unloading between different load step (see par. 8.2 of Teng and Rotter (1989a))
            if include_nonlinear_M
                el_counter = 0;
                for S = 1:length(SEGMENTS)
                    for E = 1:SEGMENTS(S).els
                        el_counter = el_counter + 1;
                        for I = 1:SOLVER.No_Gauss_Stations
                            for J = 1:SOLVER.No_Simpson_Stations
                                sigma_ys(el_counter,I,J) = MATERIALS(SEGMENTS(S).material).sigmay(sigma_ys(el_counter,I,J),epns(el_counter,I,J),Depnsj(el_counter,I,J));
                            end
                        end
                    end
                end
            end
            Depnsj(Depnsj<0) = 0;
            epns = epns + Depnsj;

            % Store the converged results as a checkpoint in case a later load-step fails to converge. Also extract the DOF values corresponding to DOF trackers.
            DELTAconv = DELTA; elDOFsconv = elDOFs; SIGMASconv = SIGMAS; Jim1 = max(j - 1,1);
            if include_nonlinear_M; DTsconv = DTs; sigmasconv = sigmas; sigma_ysconv = sigma_ys; alphasconv = alphas; epsilonsconv = epsilons; epsilonsZconv = epsilonsZ; end
            for Dkey=DOF_TRACKERS.keys
                DOF_TRACKERS(Dkey{1}) = DOF_TRACKERS(Dkey{1}).extract_DOF_LPF_pair(DELTA,SOLVER.LPF);
            end

            if SOLVER.visualise_nonlinear % DOF Trackers plot
                AQUINAS_Output_Mode_Object.nonlinear_visualiser(axDofT,'translationalDOFtrackers',i,j,[],[],[],[],[],SOLVER,[],DOF_TRACKERS,[],[]);
                AQUINAS_Output_Mode_Object.nonlinear_visualiser(axDofR,'rotationalDOFtrackers',i,j,[],[],[],[],[],SOLVER,[],DOF_TRACKERS,[],[]);
            end
            break

        elseif j == SOLVER.Jmax || (j > 1 && (norm(zeta,Inf) > SOLVER.zeta_d)) || any(isnan(RESIDUAL))

            if i == 1; error('AQUINAS Error: GNA/MNA/GMNA fails to converge from the very first load step. Consider reducing the initial load increment (dLPF property) to be applied.'); end
            if increment_attempts < SOLVER.max_attempts
                % Undo the current increment
                DELTA = DELTAconv; elDOFs = elDOFsconv; SIGMAS = SIGMASconv;
                if include_nonlinear_M; DTs = DTsconv; sigmas = sigmasconv; sigma_ys = sigma_ysconv; alphas = alphasconv; epsilons = epsilonsconv; end
                SOLVER.LPF = GMNAstep{i-1}.LPF;
                if strcmp(SOLVER.nonlinear_solver,'NewtonRaphson')
                    % Redo the load step with half the change in the LPF
                    SOLVER.dLPF = SOLVER.dLPF*SOLVER.reattempt_cutback;
                elseif strcmp(SOLVER.nonlinear_solver,'ArcLength')
                    % Redo the load step with half the arc-length
                    SOLVER.l = SOLVER.l*SOLVER.reattempt_cutback;
                end
            else
                doOneMoreLoadStep = false;
            end
            i = i - 1; Jim1 = nan;
            break

        end

        % update fields for next iteration of nonlinear procedure
        Dijm1 = Dij;
        if j>1; norm_zeta_prev = norm(zeta,Inf); end
    end
    if SOLVER.flag < 0 && increment_attempts <= SOLVER.max_attempts
        % Undo the current increment
        DELTA = DELTAconv; elDOFs = elDOFsconv; SIGMAS = SIGMASconv;
        if include_nonlinear_M; DTs = DTsconv; sigmas = sigmasconv; epsilonsZ = epsilonsZconv; sigma_ys = sigma_ysconv; alphas = alphasconv; epsilons = epsilonsconv; end
        SOLVER.LPF = GMNAstep{i-1}.LPF;
        % Redo the load step with half the arc-length
        SOLVER.l = SOLVER.l*SOLVER.reattempt_cutback;
        i = i - 1; Jim1 = nan;
    end
    % Console output for GNA/MNA/GMNA step
    if SOLVER.console_output
        if step_converged
            disp([ANALYSIS.type,' increment : ',num2str(i),' converged after ',num2str(j),' iterations']);
            if i == 1; disp(['    LPF = ',num2str(SOLVER.LPF),' [-] | dLPF = 0 [-]',]); else; disp(['    LPF = ',num2str(SOLVER.LPF),'[-] | dLPF = ',num2str(GMNAstep{i}.dLPF),' [-]',]); end
            if ~isempty(EIGVALUES); fprintf('    Lowest Eigenvalue: %.6f for a circumferential mode n = %i\n\n',min(EIGVALUES(1,:)),ANALYSIS.circumferential_modes(EIGVALUES(1,:)==min(EIGVALUES(1,:)))); end
        else
            disp([ANALYSIS.type,' attempt : ',num2str(increment_attempts + 1),' at increment : ',num2str(i),' failed after ',num2str(j),' iterations']);
        end
    end
    if step_converged; increment_attempts = 0; dDelta_a_prev = dDELTA_a; else; increment_attempts = increment_attempts + 1; end
    % Criteria to stop execution of the Geometrically Nonlinear Analysis
    [SOLVER,termination_cause] = SOLVER.nonlinear_analysis_termination(i,j,step_converged,increment_attempts,zeta,RESIDUAL,EIGVALUES,DOF_TRACKERS,GMNAstep);
    if SOLVER.flag ~= 0; doOneMoreLoadStep = false; end
end
GMNAstep{end}.termination_cause = termination_cause;


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
