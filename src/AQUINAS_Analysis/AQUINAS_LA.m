function [ DELTA, VECTOR, KSTORAGE, FSTORAGE, OSTORAGE ] = AQUINAS_LA(SEGMENTS,CONSTRAINTS,NODAL_FORCES,DISTRIBUTED_PRESSURES,SOLVER,MATERIALS)
% LA - Linear Analysis with AQUINAS
% AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MATRIX-BUILDING LOOP %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 0; % only the axisymmetric mode (n = 0) to be considered
seg_counter = 0; % Segment counter
dof_sum = 0; % Dof counter
tot_els = 0; % Total element count
el_counter = 0; % Element counter
for S = 1:length(SEGMENTS); tot_els = tot_els + SEGMENTS(S).els; end
MATRIX = cell(1,length(SEGMENTS)); VECTOR = cell(1,length(SEGMENTS));
KSTORAGE = zeros(12,12,tot_els); FSTORAGE = zeros(12,tot_els); OSTORAGE = zeros(12,tot_els);

tloopstart = tic;
for S = 1:length(SEGMENTS)
    seg_counter = seg_counter + 1;
    el_range = (1:SEGMENTS(S).els) + el_counter;
    seg_els = SEGMENTS(S).els; Kc_siz = 8; Kc_hsiz = 4;
    t = SEGMENTS(S).thickness; % Element thickness

    [MTV,VCV,offsets] = AQUINAS_initialise_offsets(seg_els,dof_sum,Kc_siz,Kc_hsiz);

    DT = MATERIALS(SEGMENTS(S).material).DT(t,[],[]); % Element tangent modulus matrix DT according to eq. 40 of Teng and Rotter (1989a)

    if strcmp(SOLVER.compiler,'C++') % Assembly via C++ MEX capability
        SEG_type = uint32(SEGMENTS(S).type_ToInt());

        cDistrPressType = uint32(ones(max(length(SEGMENTS(S).distributed_pressures),1),1));
        cDistrPressMagn = zeros(max(length(SEGMENTS(S).distributed_pressures),1),length(SEGMENTS(S).r));
        for DP=1:length(SEGMENTS(S).distributed_pressures)
            magnitudes = @(r,z,phi) DISTRIBUTED_PRESSURES(SEGMENTS(S).distributed_pressures(DP)).magnitude_arrayfun(r,z,phi);
            cDistrPressType(DP) = uint32(DISTRIBUTED_PRESSURES(SEGMENTS(S).distributed_pressures(DP)).type_ToInt);
            cDistrPressMagn(DP,:) = magnitudes(SEGMENTS(S).r,SEGMENTS(S).z,zeros(1,length(SEGMENTS(S).z)));
        end

        [MTV(:,3), VCV(:,2), KSTORAGE(:,:,el_range), FSTORAGE(:,el_range), OSTORAGE(:,el_range)] = AQUINAS_Cpp_interface_LA(...
             uint32(SEGMENTS(S).els),...
             SEGMENTS(S).r,...
             SEGMENTS(S).z,...
             SEGMENTS(S).s,...
             SEGMENTS(S).drds,...
             SEGMENTS(S).dzds,...
             SEGMENTS(S).phi,...
             SEGMENTS(S).dphids,...
             SEGMENTS(S).d2phids2,...
             uint32(length(SEGMENTS(S).distributed_pressures)),...
             cDistrPressType,...
             cDistrPressMagn,...
             uint32(strcmp(SEGMENTS(S).formulation,'thick')),...
             SEGMENTS(S).thickness,...
             [MATERIALS(SEGMENTS(S).material).nu MATERIALS(SEGMENTS(S).material).E MATERIALS(SEGMENTS(S).material).G],...
             uint32(SOLVER.No_Gauss_Stations),...
             MTV(:,3),...
             VCV(:,2),...
             uint32(sort(cell2mat(offsets.values))),...
             uint32(SOLVER.No_Threads));

        el_counter = el_range(end);

    else % Assembly within Matlab

        for E = 1:SEGMENTS(S).els
            el_counter = el_counter + 1;
            nod1 = E; nod2 = E+1; % Store bottom - top node IDs. The element's start node is the top one and its end the bottom one
            r1 = SEGMENTS(S).r(nod1); z1 = SEGMENTS(S).z(nod1); % Global r-z coordinates of element top node (Point 1 in Eqs 1a-b in Rotter & Teng (Vol. 31))
            r2 = SEGMENTS(S).r(nod2); z2 = SEGMENTS(S).z(nod2); % Global r-z coordinates of element bottom node (Point 2 in Eqs 1a-b in Rotter & Teng (Vol. 31))
            drds2 = SEGMENTS(S).drds(nod2); dzds2 = SEGMENTS(S).dzds(nod2); % dr/ds & dz/ds of element bottom node  (Point 2 in Eqs 1a-b in Rotter & Teng (Vol. 31))
            drds1 = SEGMENTS(S).drds(nod1); dzds1 = SEGMENTS(S).dzds(nod1); % dr/ds & dz/ds of element top node  (Point 1 in Eqs 1a-b in Rotter & Teng (Vol. 31))
            phi2 = SEGMENTS(S).phi(nod2); phi1 = SEGMENTS(S).phi(nod1); % phi angle at element bottom and top nodes (Points 2 and 1 accordingly in eq 1c in Rotter & Teng (Vol. 31))
            dphids2 = SEGMENTS(S).dphids(nod2); d2phids2_2 = SEGMENTS(S).d2phids2(nod2); % dphi / ds & d2phi / ds2 of element bottom node (Point 2 in eq 1c in Rotter & Teng (Vol. 31))
            dphids1 = SEGMENTS(S).dphids(nod1); d2phids2_1 = SEGMENTS(S).d2phids2(nod1); % dphi / ds & d2phi / ds2 of element top node (Point 1 in eq 1c in Rotter & Teng (Vol. 31))
            L = abs(SEGMENTS(S).s(nod2)-SEGMENTS(S).s(nod1))*0.5; % Element arc half-length, computed as the difference between the arc-length coordinates of bottom and top nodes

            % Segment psi factor regarding the importance of transverse shear strains on the overall response of the segment, when a thick formulation is employed
            if strcmp(SEGMENTS(S).formulation,'thin')
                psi = 0;
            elseif strcmp(SEGMENTS(S).formulation,'thick')
                psi = (MATERIALS(SEGMENTS(S).material).E/MATERIALS(SEGMENTS(S).material).G)*(t^2)*(6/5/4/((L)^2));
            end

            dat = struct('z2',z2,'z1',z1,'r2',r2,'r1',r1,'L',L,'phi2',phi2,'phi1',phi1,'drds2',drds2,'drds1',drds1,'dzds2',dzds2,'dzds1',dzds1,...
            'dphids2',dphids2,'dphids1',dphids1,'d2phids2_2',d2phids2_2,'d2phids2_1',d2phids2_1,'t',t,'psi',psi);

            % In the following the integrands are evaluated trhough the Gaussian integration scheme
            % intA () dA = intS intTheta ()*r(eta) dTheta dS = intEta intTheta ()*L*r(eta) dTheta dEta
            KBARE = zeros(12); FBARE = zeros(12,1);
            for DP=SEGMENTS(S).distributed_pressures
                FBARE = FBARE + DISTRIBUTED_PRESSURES(DP).eta_integral_equivalent_nodal_load_vector(dat);
            end
            for I = 1:length(SOLVER.Gauss_Nodes)
                eta = SOLVER.Gauss_Nodes(I);
                KBARE = KBARE + SOLVER.Gauss_Weights(I)*L*AQUINAS_B0DTB0_matrix(eta,0,n,dat,DT); % material stiffness matrix, based on linear strains, for any circumferential mode (here only for n = 0 / theta = 0, axisymmetric mode)
            end

            [MTV,VCV] = AQUINAS_condense_element_matrix('material',KBARE,[],FBARE,phi1,phi2,E,n,offsets,MTV,VCV); % This performs the transformation and condensation operations - Eqs 63
            KSTORAGE(:,:,el_counter) = KBARE;
            FSTORAGE(:,el_counter) = FBARE;
            OSTORAGE(:,el_counter) = [drds2 drds1 dzds2 dzds1 phi2 phi1 dphids2 dphids1 d2phids2_2 d2phids2_1 L psi]';
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
VECTOR = sparse(VECTOR(:,1),ones(size(VECTOR(:,1))),VECTOR(:,2),dof_sum,1);

% Exclude stiffness entries corresponding to v DOFs
v_free_indices = sort([1:4:SEGMENTS(end).bot_dofIDs(1) 3:4:SEGMENTS(end).bot_dofIDs(3) 4:4:SEGMENTS(end).bot_dofIDs(4) SEGMENTS(end).bot_dofIDs(4)+1:dof_sum]);
DELTA = zeros(size(VECTOR));

tassembly = toc(tloopstart);

% Solution of the equilibrium problem
tsolvestart = tic;
DELTA(v_free_indices) = MATRIX(v_free_indices,v_free_indices)\VECTOR(v_free_indices);
tsolver = toc(tsolvestart);

if SOLVER.console_output
    disp('LA step (n = 0) completed.');
    disp(['    Assembly time: ',num2str(tassembly),' [s]']);
    disp(['    Solver time: ',num2str(tsolver),' [s]']);
    disp(' ');
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
