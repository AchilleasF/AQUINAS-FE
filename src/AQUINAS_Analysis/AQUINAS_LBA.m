function [ EIGVECTORS, EIGVALUES ] = AQUINAS_LBA(SEGMENTS,CONSTRAINTS,SOLVER,MATERIALS,ANALYSIS,DTs,elDOFs,Sigmas,OStorage)
% AQUINAS_LBA - Linear Bifurcation Analysis with AQUINAS
% AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.

if isempty(ANALYSIS.circumferential_modes); error('AQUINAS Error: Analysis Object definition must include circumferential modes argument for LBAs'); end
if isempty(ANALYSIS.no_eigenvalues); error('AQUINAS Error: Analysis Object definition must include number of eigenvalues argument for LBAs'); end


dof_count = 0;
for S = 1:length(SEGMENTS); dof_count = dof_count + (SEGMENTS(S).els + 1) * 4; end
for C = 1:length(CONSTRAINTS); if (strcmp(CONSTRAINTS(C).state,'general') || strcmp(CONSTRAINTS(C).state,'buckling')); dof_count = dof_count + 1; end; end
if (strcmp(ANALYSIS.type,'GNA') || strcmp(ANALYSIS.type,'GMNA')); include_nonlinear_G = true; else; include_nonlinear_G = false; end
if strcmp(ANALYSIS.type,'GMNA'); include_nonlinear_M = true; else; include_nonlinear_M = false; end
EIGVECTORS = zeros(dof_count, ANALYSIS.no_eigenvalues, length(ANALYSIS.circumferential_modes));
EIGVALUES = zeros(ANALYSIS.no_eigenvalues, length(ANALYSIS.circumferential_modes));

for CM = 1:length(ANALYSIS.circumferential_modes)
    tloopstart = tic;
    C_mode = ANALYSIS.circumferential_modes(CM);

    seg_counter = 0; % Segment counter
    dof_sum = 0; % Dof counter
    tot_els = 0; % Total element count
    el_counter = 0; % Element counter
    for S = 1:length(SEGMENTS); tot_els = tot_els + SEGMENTS(S).els; end
    MATRIX_M = cell(1,length(SEGMENTS));
    MATRIX_G = cell(1,length(SEGMENTS));
    KSTORAGE_M = zeros(12,12,tot_els);
    KSTORAGE_G = zeros(12,12,tot_els);

    for S = 1:length(SEGMENTS)
        seg_counter = seg_counter + 1;
        el_range = (1:SEGMENTS(S).els) + el_counter;
        seg_els = SEGMENTS(S).els; Kc_siz = 8; Kc_hsiz = 4;
        t = SEGMENTS(S).thickness; % Element thickness

        [MTV_M,~,offsets] = AQUINAS_initialise_offsets(seg_els,dof_sum,Kc_siz,Kc_hsiz); MTV_G = MTV_M;

        if strcmp(SOLVER.compiler,'C++') % Assembly via C++ MEX capability
            if strcmp(SEGMENTS(S).type,'Cone'); SEG_type = 1; end
            if strcmp(SEGMENTS(S).type,'Plate'); SEG_type = 2; end
            if strcmp(SEGMENTS(S).type,'Ellipse'); SEG_type = 3; end

            if ~include_nonlinear_M; DTs = MATERIALS(SEGMENTS(S).material).DT(t,[],[]); end
            [MTV_M(:,3), MTV_G(:,3), KSTORAGE_M(:,:,el_range), KSTORAGE_G(:,:,el_range)] = AQUINAS_Cpp_interface_LBA(...
                uint32(SEGMENTS(S).els),...
                SEGMENTS(S).r,...
                SEGMENTS(S).z,...
                SEGMENTS(S).thickness,...
                uint32(C_mode),...
                uint32(SOLVER.No_Gauss_Stations),...
                uint32(el_counter),...
                elDOFs,...
                Sigmas,...
                DTs,...
                uint32(sort(cell2mat(offsets.values))),...
                uint32(include_nonlinear_G),...
                uint32(include_nonlinear_M),...
                MTV_M(:,3),...
                MTV_G(:,3),...
                OStorage(:,el_range),...
                uint32(SOLVER.No_Threads));

            el_counter = el_range(end);

        else % Assembly within Matlab
            if ~include_nonlinear_M; DT = MATERIALS(SEGMENTS(S).material).DT(t,[],[]); end
            for E = 1:SEGMENTS(S).els
                el_counter = el_counter + 1;
                r1 = SEGMENTS(S).r(E); z1 = SEGMENTS(S).z(E); % Global r-z coordinates of element top node (Point 1 in Eqs 1a-b in Rotter & Teng (Vol. 31))
                r2 = SEGMENTS(S).r(E+1); z2 = SEGMENTS(S).z(E+1); % Global r-z coordinates of element bottom node (Point 2 in Eqs 1a-b in Rotter & Teng (Vol. 31))
                drds2 = OStorage(1, el_counter); drds1 = OStorage(2, el_counter);
                dzds2 = OStorage(3, el_counter); dzds1 = OStorage(4, el_counter);
                phi2 = OStorage(5, el_counter); phi1 = OStorage(6, el_counter);
                dphids2 = OStorage(7, el_counter); dphids1 = OStorage(8, el_counter);
                d2phids2_2 = OStorage(9, el_counter); d2phids2_1 = OStorage(10, el_counter);
                L = OStorage(11, el_counter); psi = OStorage(12, el_counter);
                dat = struct('elDOFs',elDOFs(:,el_counter),'z2',z2,'z1',z1,'r2',r2,'r1',r1,'L',L,'phi2',phi2,'phi1',phi1,'drds2',drds2,'drds1',drds1,'dzds2',dzds2,'dzds1',dzds1,...
                    'dphids2',dphids2,'dphids1',dphids1,'d2phids2_2',d2phids2_2,'d2phids2_1',d2phids2_1,'t',t,'psi',psi);

                % intA () dA = intS intTheta ()*r(eta) dTheta dS = intEta intTheta ()*L*r(eta) dTheta dEta
                KBARE_M = zeros(12,12); KBARE_G = zeros(12,12);
                for I = 1:length(SOLVER.Gauss_Nodes)
                    eta = SOLVER.Gauss_Nodes(I);
                    if include_nonlinear_M; DT = squeeze(DTs(el_counter,I,:,:)); end
                    % Two integrals are evaluated for non zero circumferential wavenumbers, to take account of the integration in
                    % the circumferential direction by observing that each of the non zero terms of equations 67
                    % of Rotter and Teng (1989b) (Vol 32) can be evaluated idependently for theta = 0 and theta = pi/2 accordingly
                    if C_mode==0
                        if include_nonlinear_G
                            KBARE_M = KBARE_M + SOLVER.Gauss_Weights(I) * L * AQUINAS_BDTB_matrix(eta,0,C_mode,dat,DT);
                        else
                            KBARE_M = KBARE_M + SOLVER.Gauss_Weights(I) * L * AQUINAS_B0DTB0_matrix(eta,0,C_mode,dat,DT);
                        end
                        KBARE_G = KBARE_G + SOLVER.Gauss_Weights(I) * L * AQUINAS_GN0G_matrix(eta,0,C_mode,Sigmas(:,el_counter,I),dat);
                    else
                        if include_nonlinear_G
                            KBARE_M = KBARE_M + SOLVER.Gauss_Weights(I) * L * (AQUINAS_BDTB_matrix(eta,0,C_mode,dat,DT) + ...
                                AQUINAS_BDTB_matrix(eta,pi/2/C_mode,C_mode,dat,DT));
                        else
                            KBARE_M = KBARE_M + SOLVER.Gauss_Weights(I) * L * (AQUINAS_B0DTB0_matrix(eta,0,C_mode,dat,DT) + ...
                                AQUINAS_B0DTB0_matrix(eta,pi/2/C_mode,C_mode,dat,DT));
                        end
                        KBARE_G = KBARE_G + SOLVER.Gauss_Weights(I) * L * (AQUINAS_GN0G_matrix(eta,0,C_mode,Sigmas(:,el_counter,I),dat) + ...
                            AQUINAS_GN0G_matrix(eta,pi/2/C_mode,C_mode,Sigmas(:,el_counter,I),dat));
                    end
                end
                [MTV_M,~] = AQUINAS_condense_element_matrix('material',KBARE_M,[],[],phi1,phi2,E,C_mode,offsets,MTV_M,[]); % This performs the transformation and condensation operations - Eqs 63
                [MTV_G,~] = AQUINAS_condense_element_matrix('geometric',KBARE_M,KBARE_G,[],phi1,phi2,E,C_mode,offsets,MTV_G,[]); % This performs the transformation and condensation operations - Eqs 63
                KSTORAGE_M(:,:,el_counter) = KBARE_M;
                KSTORAGE_G(:,:,el_counter) = KBARE_G;
            end
        end
        MATRIX_M(seg_counter) = {MTV_M};
        MATRIX_G(seg_counter) = {MTV_G};
        dof_sum = dof_sum + (seg_els + 1) * Kc_hsiz;
    end


    % Adding matrix entries for Lagrange multipliers
    for C = 1:length(CONSTRAINTS)
        if (C_mode==0) && strcmp(CONSTRAINTS(C).origin,'user') && (length(CONSTRAINTS(C).local_dofs)==1) && strcmp(CONSTRAINTS(C).local_dofs{1},'v'); continue; end
        if (C_mode==0) && strcmp(CONSTRAINTS(C).origin,'auto') && (length(CONSTRAINTS(C).global_dofs{:})==2) && (any(mod(CONSTRAINTS(C).global_dofs{:},4)==2)); continue; end
        if ~(strcmp(CONSTRAINTS(C).state,'general')||strcmp(CONSTRAINTS(C).state,'buckling')); continue; end
        for I = 1:length(CONSTRAINTS(C).global_dofs)
            dof_sum = dof_sum + 1; seg_counter = seg_counter + 1;
            rowM = dof_sum*ones(length(CONSTRAINTS(C).global_dofs{I}),1);
            colM = CONSTRAINTS(C).global_dofs{I}';
            valM = CONSTRAINTS(C).coeffLHS{I}';
            MATRIX_M(seg_counter) = {[[rowM; colM] [colM; rowM] [valM; valM]]};
            MATRIX_G(seg_counter) = {[[rowM; colM] [colM; rowM] [SOLVER.Constr_Langr_Multip_Factor*valM; SOLVER.Constr_Langr_Multip_Factor*valM]]};
        end
    end

    % The complete material stiffness matrix of the system for harmonic n
    MATRIX_M = cell2mat(MATRIX_M');
    MATRIX_M = sparse(MATRIX_M(:,1),MATRIX_M(:,2),MATRIX_M(:,3),dof_sum,dof_sum);

    % The complete eometric stiffness matrix of the system for harmonic n
    MATRIX_G = cell2mat(MATRIX_G');
    MATRIX_G = sparse(MATRIX_G(:,1),MATRIX_G(:,2),MATRIX_G(:,3),dof_sum,dof_sum);
    % Exclude stiffness entries corresponding to v DOFs
    v_free_indices = sort([1:4:SEGMENTS(end).bot_dofIDs(1) 3:4:SEGMENTS(end).bot_dofIDs(3) 4:4:SEGMENTS(end).bot_dofIDs(4) SEGMENTS(end).bot_dofIDs(4)+1:dof_sum]);
    tassembly = toc(tloopstart);

    % Solution of the linearised bifurcation problem
    tsolvestart = tic;
    try
        if C_mode == 0
            [EIGVECTORS(v_free_indices,:,CM), tmp_eigenvalues] = eigs(MATRIX_M(v_free_indices,v_free_indices), -MATRIX_G(v_free_indices,v_free_indices), ANALYSIS.no_eigenvalues, SOLVER.eigsSigma,...
                'Tolerance', SOLVER.eigsTolerance, 'SubspaceDimension', max(max(2*ANALYSIS.no_eigenvalues,20),min(SOLVER.eigsSubspaceDimensionMultiplier*ANALYSIS.no_eigenvalues,length(v_free_indices))));
        else
            [EIGVECTORS(:,:,CM), tmp_eigenvalues] = eigs(MATRIX_M, -MATRIX_G, ANALYSIS.no_eigenvalues, SOLVER.eigsSigma,...
                'Tolerance', SOLVER.eigsTolerance, 'SubspaceDimension', max(max(2*ANALYSIS.no_eigenvalues,20),min(SOLVER.eigsSubspaceDimensionMultiplier*ANALYSIS.no_eigenvalues,length(MATRIX_M))));
        end
        if strcmp(SOLVER.eigsNegativeTreatment,'Drop')
            tmp_eigenvalues(tmp_eigenvalues<0) = nan;
        end
        EIGVALUES(:,CM) = diag(tmp_eigenvalues);
        if SOLVER.eigsForceReal
            if sum(~isreal(EIGVALUES(:,CM)))>0 || sum(sum(~isreal(EIGVECTORS(:,:,CM))))>0
                warning(strcat("A complex eigenvalue has been found, for a circumferential wavenumber of n = ",char(32),num2str(C_mode),char(32)," with a maximum absolute ratio of imaginary/real parts of ",char(32),num2str(max(abs(imag(EIGVALUES(:,CM)))./abs(real(EIGVALUES(:,CM)))))));
                EIGVECTORS(:,:,CM) = real(EIGVECTORS(:,:,CM));
                EIGVALUES(:,CM) = real(EIGVALUES(:,CM));
            end
        end
    catch
        EIGVECTORS(:,:,CM) = NaN;
        EIGVALUES(:,CM) = NaN;
    end
    tsolver = toc(tsolvestart);

    % Normalize eigenvectors by dividing with the second norm of each eigenvector, taking into account only the terms up to the last DOF
    if strcmp(ANALYSIS.normalize_eigs,'SecondNorm')
        lastDOF = SEGMENTS(end).bot_dofIDs(4);
        for i=1:ANALYSIS.no_eigenvalues
            EIGVECTORS(1:lastDOF,i,CM)=EIGVECTORS(1:lastDOF,i,CM)/vecnorm(EIGVECTORS(1:lastDOF,i,CM));
        end
    elseif strcmp(ANALYSIS.normalize_eigs,'MaxValue')
        lastDOF = SEGMENTS(end).bot_dofIDs(4);
        for i=1:ANALYSIS.no_eigenvalues
            if ~all(isnan(EIGVECTORS(1:lastDOF,i,CM)))
                EIGVECTORS(1:lastDOF,i,CM)=EIGVECTORS(1:lastDOF,i,CM)/EIGVECTORS(find(abs(EIGVECTORS(1:lastDOF,i,CM))==max(abs(EIGVECTORS(1:lastDOF,i,CM))),1),i,CM);
            end
        end
    end

    if SOLVER.console_output
        if strcmp(ANALYSIS.type,'LBA'); disp(['LBA step (n = ',num2str(C_mode),') completed.']); else; disp(['Bifurcation step (n = ',num2str(C_mode),') completed.']); end
        disp(['    Assembly time: ',num2str(tassembly),' [s]']);
        disp(['    Solver time: ',num2str(tsolver),' [s]']);
        fprintf('    Lowest Eigenvalue: %.6f\n',EIGVALUES(1,CM))
        disp(' ');
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
