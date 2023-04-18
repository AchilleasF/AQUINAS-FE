function [ DataOut ] = AQUINAS_post_process(DELTA,VECTOR,KSTORAGE,FSTORAGE,OSTORAGE,SEGMENTS,MATERIALS,ANALYSIS,SOLVER,DataIn)
% AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.

switch ANALYSIS.type
    case 'LA'

        % Generation of strains, stresses and generalised stress resultants, evaluated at the nodes
        dofs = 0;
        el_counter = 0;
        DOF_Data = zeros(SEGMENTS(end).bot_dofIDs(4)*6/4, 1);
        midStrain_Data = zeros(SEGMENTS(end).bot_dofIDs(4)*6/4, 1);
        Strain_Data = zeros(SEGMENTS(end).bot_dofIDs(4)*9/4, 1);
        Stress_Reslt_Data = zeros(SEGMENTS(end).bot_dofIDs(4)*6/4, 1);
        Stress_Data = zeros(SEGMENTS(end).bot_dofIDs(4)*9/4, 1);
        vonMisesStress_Data = zeros(SEGMENTS(end).bot_dofIDs(4)*3/4, 1);

        for S = 1:length(SEGMENTS)
            if (S > 1); dofs = dofs + 4; end
            t = SEGMENTS(S).thickness;
            DT = MATERIALS(SEGMENTS(S).material).DT(t,[],[]); % Element tangent modulus matrix DT according to eq. 40 of Teng and Rotter (1989a)

            for E = 1:SEGMENTS(S).els
                el_counter = el_counter + 1;
                dat.t = t;
                dat.r1 = SEGMENTS(S).r(E); dat.z1 = SEGMENTS(S).z(E);
                dat.r2 = SEGMENTS(S).r(E+1); dat.z2 = SEGMENTS(S).z(E+1);
                dat.drds2 = OSTORAGE(1, el_counter); dat.drds1= OSTORAGE(2, el_counter);
                dat.dzds2 = OSTORAGE(3, el_counter); dat.dzds1= OSTORAGE(4, el_counter);
                dat.phi2 = OSTORAGE(5, el_counter); dat.phi1 = OSTORAGE(6, el_counter);
                dat.dphids2 = OSTORAGE(7, el_counter); dat.dphids1 = OSTORAGE(8, el_counter);
                dat.d2phids2_2 = OSTORAGE(9, el_counter); dat.d2phids2_1 = OSTORAGE(10, el_counter);
                dat.L = OSTORAGE(11, el_counter);
                dat.psi = OSTORAGE(12, el_counter);

                elDOFs = AQUINAS_vaporise_element_matrix(KSTORAGE(:,:,el_counter),FSTORAGE(:,el_counter),DELTA(dofs+(1:8)),dat.phi1,dat.phi2,0);
                dat.elDOFs = elDOFs;
                DOF_Data(dofs*1.5 + (1:12)) = elDOFs;

                if E==1
                    eta = -1;
                    if (dat.r1 < 0.05*dat.r2); eta = -0.95; end
                    % Evaluation of midsurface strains and curvatures
                    midStrain_Data(dofs*6/4 + (1:6)) = AQUINAS_strains_vector(eta,dat,false);

                    % Evaluation of membrane and bending stress resultants
                    Stress_Reslt_Data(dofs*6/4 + (1:6)) = DT * midStrain_Data(dofs*6/4 + (1:6));

                    % Evaluation  of inner/mid/outer surface stresses
                    Stress_Data(dofs*(9/4) + [1 2 3]) = Stress_Reslt_Data(dofs*6/4 + [1 2 3])/t - 6*Stress_Reslt_Data(dofs*6/4 + [4 5 6])/(t*t);
                    Stress_Data(dofs*(9/4) + [4 5 6]) = Stress_Reslt_Data(dofs*6/4 + [1 2 3])/t;
                    Stress_Data(dofs*(9/4) + [7 8 9]) = Stress_Reslt_Data(dofs*6/4 + [1 2 3])/t + 6*Stress_Reslt_Data(dofs*6/4 + [4 5 6])/(t*t);

                    % Evaluation  of inner/mid/outer surface Von Mises stresses
                    vonMisesStress_Data(dofs*(3/4) + 1) = sqrt((Stress_Data(dofs*(9/4)+1)^2) + (Stress_Data(dofs*(9/4)+2)^2) - (Stress_Data(dofs*(9/4)+1)*Stress_Data(dofs*(9/4)+2)) + 3*(Stress_Data(dofs*(9/4)+3)^2));
                    vonMisesStress_Data(dofs*(3/4) + 2) = sqrt((Stress_Data(dofs*(9/4)+4)^2) + (Stress_Data(dofs*(9/4)+5)^2) - (Stress_Data(dofs*(9/4)+4)*Stress_Data(dofs*(9/4)+5)) + 3*(Stress_Data(dofs*(9/4)+6)^2));
                    vonMisesStress_Data(dofs*(3/4) + 3) = sqrt((Stress_Data(dofs*(9/4)+7)^2) + (Stress_Data(dofs*(9/4)+8)^2) - (Stress_Data(dofs*(9/4)+7)*Stress_Data(dofs*(9/4)+8)) + 3*(Stress_Data(dofs*(9/4)+9)^2));

                    % Evaluation  of inner/mid/outer surface strains
                    Strain_Data(dofs*(9/4) + [1 2 3]) = midStrain_Data(dofs*6/4 + [1 2 3]) - midStrain_Data(dofs*6/4 + [4 5 6])*(t/2);
                    Strain_Data(dofs*(9/4) + [4 5 6]) = midStrain_Data(dofs*6/4 + [1 2 3]);
                    Strain_Data(dofs*(9/4) + [7 8 9]) = midStrain_Data(dofs*6/4 + [1 2 3]) + midStrain_Data(dofs*6/4 + [4 5 6])*(t/2);
                end

                eta = 1;
                if (dat.r2 < 0.05*dat.r1); eta = 0.95; end
                % Evaluation of midsurface strains and curvatures
                midStrain_Data(dofs*6/4 + (7:12)) = AQUINAS_strains_vector(eta,dat,false);

                % Evaluation of membrane and bending stress resultants
                Stress_Reslt_Data(dofs*6/4 + (7:12)) = DT * midStrain_Data(dofs*6/4 + (7:12));

                % Evaluation  of inner/mid/outer surface stresses
                Stress_Data(dofs*(9/4) + [10 11 12]) = Stress_Reslt_Data(dofs*6/4 + [7 8 9])/t - 6*Stress_Reslt_Data(dofs*6/4 + [10 11 12])/(t*t);
                Stress_Data(dofs*(9/4) + [13 14 15]) = Stress_Reslt_Data(dofs*6/4 + [7 8 9])/t;
                Stress_Data(dofs*(9/4) + [16 17 18]) = Stress_Reslt_Data(dofs*6/4 + [7 8 9])/t + 6*Stress_Reslt_Data(dofs*6/4 + [10 11 12])/(t*t);

                % Evaluation  of inner/mid/outer surface Von Mises stresses
                vonMisesStress_Data(dofs*(3/4) + 4) = sqrt((Stress_Data(dofs*(9/4)+10)^2) + (Stress_Data(dofs*(9/4)+11)^2) - (Stress_Data(dofs*(9/4)+10)*Stress_Data(dofs*(9/4)+11)) + 3*(Stress_Data(dofs*(9/4)+12)^2));
                vonMisesStress_Data(dofs*(3/4) + 5) = sqrt((Stress_Data(dofs*(9/4)+13)^2) + (Stress_Data(dofs*(9/4)+14)^2) - (Stress_Data(dofs*(9/4)+13)*Stress_Data(dofs*(9/4)+14)) + 3*(Stress_Data(dofs*(9/4)+15)^2));
                vonMisesStress_Data(dofs*(3/4) + 6) = sqrt((Stress_Data(dofs*(9/4)+16)^2) + (Stress_Data(dofs*(9/4)+17)^2) - (Stress_Data(dofs*(9/4)+16)*Stress_Data(dofs*(9/4)+17)) + 3*(Stress_Data(dofs*(9/4)+18)^2));

                % Evaluation  of inner/mid/outer surface strains
                Strain_Data(dofs*(9/4) + [10 11 12]) = midStrain_Data(dofs*6/4 + [7 8 9]) - midStrain_Data(dofs*6/4 + [10 11 12])*(t/2);
                Strain_Data(dofs*(9/4) + [13 14 15]) = midStrain_Data(dofs*6/4 + [7 8 9]);
                Strain_Data(dofs*(9/4) + [16 17 18]) = midStrain_Data(dofs*6/4 + [7 8 9]) + midStrain_Data(dofs*6/4 + [10 11 12])*(t/2);

                dofs = dofs + 4;
            end
        end

        DataOut = {DOF_Data midStrain_Data Strain_Data Stress_Reslt_Data Stress_Data vonMisesStress_Data VECTOR(1:(4/6)*length(DOF_Data)) DELTA(1:(4/6)*length(DOF_Data))};
        return

    case 'LBA'

        % The deformed shape of the LA solution needs to be
        % uncondensed in order to compute the resultants that lead to buckling

        numelems = 0; % Total element count
        for S = 1:length(SEGMENTS); numelems = numelems + SEGMENTS(S).els; end
        elDOFs = zeros(12,numelems);
        Stress_Reslt_Data = zeros(3,numelems,SOLVER.No_Gauss_Stations); % Nphi, Ntheta and Nphitheta per element per Gauss integration point

        seg_counter = 0; % Segment counter
        el_counter = 0; % Element counter
        dofs = 0;

        for S = 1:length(SEGMENTS)
            if S>1; dofs = dofs + 4; end
            seg_counter = seg_counter + 1;
            t = SEGMENTS(S).thickness;
            DT = MATERIALS(SEGMENTS(S).material).DT(t,[],[]); % Element tangent modulus matrix DT according to eq. 40 of Teng and Rotter (1989a) (the same for all elements when no material nonlirearities are employed, as is the case for an LBA)

            for E = 1:SEGMENTS(S).els
                el_counter = el_counter + 1;
                dat.t = t;
                dat.r1 = SEGMENTS(S).r(E); dat.z1 = SEGMENTS(S).z(E);
                dat.r2 = SEGMENTS(S).r(E+1); dat.z2 = SEGMENTS(S).z(E+1);
                dat.drds1 = OSTORAGE(2, el_counter); dat.drds2 = OSTORAGE(1, el_counter);
                dat.dzds1 = OSTORAGE(4, el_counter); dat.dzds2 = OSTORAGE(3, el_counter);
                dat.phi1 = OSTORAGE(6, el_counter); dat.phi2 = OSTORAGE(5, el_counter);
                dat.dphids1 = OSTORAGE(8, el_counter); dat.dphids2 = OSTORAGE(7, el_counter);
                dat.d2phids2_1 = OSTORAGE(10, el_counter); dat.d2phids2_2 = OSTORAGE(9, el_counter);
                dat.L = OSTORAGE(11, el_counter);
                dat.psi = OSTORAGE(12, el_counter);

                elDOFs(:,el_counter) = AQUINAS_vaporise_element_matrix(KSTORAGE(:,:,el_counter),FSTORAGE(:,el_counter),DELTA(dofs+(1:8)),dat.phi1,dat.phi2,0);
                dat.elDOFs = elDOFs(:,el_counter);

                % Evaluate membrane stress resultants per Gauss point
                for I = 1:length(SOLVER.Gauss_Nodes)
                    eta = SOLVER.Gauss_Nodes(I);
                    epsilon = AQUINAS_strains_vector(eta,dat,false);
                    Stress_Reslt_Data(:,el_counter,I) = DT(1:3,:) * epsilon;
                end
                dofs = dofs + 4;
            end

        end

        DataOut = {elDOFs Stress_Reslt_Data};
        return

    case {'GNA','MNA','GMNA'}

        if isempty(DataIn)

            numelems = 0; % Total element count
            for S = 1:length(SEGMENTS); numelems = numelems + SEGMENTS(S).els; end
            elDOFs = zeros(12,numelems);
            seg_counter = 0; % Segment counter
            el_counter = 0; % Element counter
            dofs = 0;
            for S=1:length(SEGMENTS)
                if S>1; dofs = dofs + 4; end
                seg_counter = seg_counter + 1;
                for E = 1:SEGMENTS(S).els
                    el_counter = el_counter + 1;
                    PHI2 = OSTORAGE(5, el_counter); % PHI2 = PHIB
                    PHI1 = OSTORAGE(6, el_counter); % PHI1 = PHIT
                    elDOFs(:,el_counter) = AQUINAS_vaporise_element_matrix(KSTORAGE(:,:,el_counter),FSTORAGE(:,el_counter),DELTA(dofs+(1:8)),PHI1,PHI2,0);
                    dofs = dofs + 4;
                end
            end
            DataOut = elDOFs;

        else

            DataOut = DataIn;

        end

        return

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