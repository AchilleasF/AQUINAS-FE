classdef AQUINAS_Material_Object
    % Class definition for a Material Object - used to build the stiffness matrix of
    % a single axisymmetric shell segment
    % At least one Material Object must exist
    % AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.

    properties
        class_type = 'AQUINAS_Material_Object';
        % General isotropic-orthotropic elasticity
        name % Material Object name
        type % type of material object, can be either isotropic or orthotrpic
        E % Elastic modulus in both phi and theta directions [F/L^2] for isotropic material
        Ephi % Elastic modulus in phi direction [F/L^2] for othotropic material
        Etheta % Elastic modulus in theta direction [F/L^2] for othotropic material
        nu % Poisson ratio [-] for isotropic material
        nuphi_theta % Poisson ratio with phi-theta subscript [-] for othotropic material
        nutheta_phi % Poisson ratio in theta-phi subscript [-] for othotropic material
        G % Shear modulus [F/L^2]
        De % Elastic material rigidity matrix (P1: Eq. 33; P3: Eq. 80)

        % Isotropic plasticity
        curve_def % System that the stress-strain curve provided is defined in. Can be either 'True' or 'Engineering'.
        sy = [] % Engineering yield stress [F/L^2]
        ep = [] % Engineering plastic strain [-]


        % Auxiliary fields for material nonlinearities
        yield_tol = 1e-8; % Tolerance for testing yielding
    end

    methods

        % Initiator method for a Material Object
        function obj = AQUINAS_Material_Object(options)
            arguments
                options.name { mustBeNonzeroLengthText }
                options.E { mustBeReal }
                options.Ephi { mustBeReal }
                options.Etheta { mustBeReal }
                options.nu { mustBeReal }
                options.nuphitheta { mustBeReal }
                options.nuthetaphi { mustBeReal }
                options.G { mustBeReal }
                options.sy { mustBeReal } = []
                options.ep { mustBeReal } = []
                options.curveDef { mustBeNonzeroLengthText } = 'Undefined'
            end

            if ~isfield(options,'name'); error('AQUINAS Error: Material Object definition must include a name for the material.'); end
            if ~(isfield(options,'E') || (isfield(options,'Ephi') && isfield(options,'Etheta'))); error('AQUINAS Error: Material Object must either include elastic modulus E for isotropic material or both Ephi and Etheta for orthrotropic material.'); end
            if (isfield(options,'E') && (isfield(options,'Ephi') || isfield(options,'Etheta'))); error('AQUINAS Error: Material Object must either include elastic modulus E for isotropic material or both Ephi and Etheta for orthrotropic material. The material object cannot be both isotropic and orthotrpic.'); end
            if ~(isfield(options,'nu') || (isfield(options,'nuphitheta') && isfield(options,'nuthetaphi'))); error('AQUINAS Error: Material Object must either include Poisson ratio nu for isotropic material or both nuphitheta and nuthetaphi for orthrotropic material.'); end
            if (isfield(options,'nu') && (isfield(options,'nuphitheta') || isfield(options,'nuthetaphi'))); error('AQUINAS Error: Material Object must either include Poisson ratio nu for isotropic material or both nuphitheta and nuthetaphi for orthrotropic material. The material object cannot be both isotropic and orthotrpic.'); end
            if ~(isfield(options,'E') && isfield(options,'nu')) && ~(isfield(options,'Ephi') && isfield(options,'Etheta') && isfield(options,'nuphitheta') && isfield(options,'nuthetaphi')); error('AQUINAS Error: Material Object must either include E and nu for isotropic material or Ephi, Etheta, nuphitheta and nuthetaphi for orthotropic material.'); end
            if ~isfield(options,'G'); error('AQUINAS Error: Material Object must include shear modulus'); end

            % Elastic material properties
            obj.name = options.name;
            if isfield(options,'E') && isfield(options,'nu')
                obj.type = 'isotropic';
                obj.E = options.E; obj.nu = options.nu;
            elseif isfield(options,'Ephi') && isfield(options,'Etheta') && isfield(options,'nuphitheta') && isfield(options,'nuthetaphi')
                obj.type = 'orthotropic';
                obj.Ephi = options.Ephi; obj.Etheta = options.Etheta; obj.nuphi_theta = options.nuphitheta; obj.nutheta_phi = options.nuthetaphi;
            end
            obj.G = options.G;
            % Construction of elastic material rigidity matrix [De]
            % This matrix will be invariant at every point during the analysis
            if strcmp(obj.type,'isotropic')
                EPP = obj.E/(1 - obj.nu*obj.nu); EPT = obj.nu*EPP; ETP = EPT; ETT = obj.E/(1 - obj.nu*obj.nu);
            elseif strcmp(obj.type,'orthotropic')
                EPP = obj.Ephi/(1 - obj.nuphi_theta*obj.nutheta_phi); EPT = obj.nuphi_theta*EPP; ETP = EPT; ETT = obj.Etheta/(1 - obj.nuphi_theta*obj.nutheta_phi);
            end
            obj.De = [EPP EPT 0     ;
                      ETP ETT 0     ;
                      0   0   obj.G];
            % Elastoplastic material properties
            if ~isempty(options.ep)
                if options.ep(1) ~= 0; error('AQUINAS Error : The plastic strain at the point of first yield must be zero.'); end
                if length(options.ep) ~= length(options.sy); error('AQUINAS Error : The yield stress - plastic strain vectors provided for the definition of the elastoplastic stress-strain curve must be of equal length.'); end
            end
            % Construction of post-yield true stress-true strain curve, in case the input stresses-strains are engineering
            if strcmp(options.curveDef,'Engineering') && (~isempty(options.sy) || ~isempty(options.ep))
                obj.sy = options.sy.*(1 + options.ep); obj.ep = log(1 + options.ep); % Conversion to true stress & strain
            elseif strcmp(options.curveDef,'True') && (~isempty(options.sy) || ~isempty(options.ep))
                obj.sy = options.sy; obj.ep = options.ep; % No need for conversion, the input curves are already true stresses - strains
            end
            obj.curve_def = options.curveDef;

        end


        % Compute strain hardening parameter H, according to eq. 27 of Teng and Rotter (1989a)
        function H = H(obj,epn)
            if length(obj.ep) <= 1; H = 0; return; end
            % Find the range that of the stres-strain curve that the current equivalnent plastic strain belongs to
            for I = 1:length(obj.ep)
                if I == length(obj.ep); break; end
                if obj.ep(I) < epn && obj.ep(I+1) > epn; break; end
            end
            if I == length(obj.ep); I = I - 1; end
            ep_p = obj.ep(I); ep_n = obj.ep(I+1);
            sy_p = obj.sy(I); sy_n = obj.sy(I+1);
            % Compute strain hardening parameter
            H = (sy_n - sy_p)/(ep_n - ep_p);
        end

        % Compute yield stress of material point, according to par. 10-(8)-(v)-(b) of Teng and Rotter (1989a)
        function sigy = sigmay(obj,sigy,epn,depn)
            if epn + depn < 1e-16; return; end % if there is no equivalnet plastic strain developing, do not alter the value of the yield stress
            sigy = sigy + obj.H(epn)*depn;
        end

        % Compute Dep matrix at a through-thickness station of coordinate z, according to eq. 34 of Teng and Rotter (1989a)
        function Dep = Dep(obj,sigma,sigmay,epn)
            % von Mises equivalent stress, according to eq. 24 of Teng and Rotter (1989a)
            sigmabar = sqrt(sigma(1)^2 + sigma(2)^2 - sigma(1)*sigma(2) + 3*sigma(3)^2);
            if sigmabar > sigmay*(1-obj.yield_tol)
                % Deviatoric sigmas, according to eqs. 30 of Teng and Rotter (1989a)
                sphi = (2*sigma(1) - sigma(2))/3; stheta = (2*sigma(2) - sigma(1))/3;
                % S1, S2, S4 and S5 parameters, according to eqs. 35 of Teng and Rotter (1989a)
                S1 = sphi + obj.nu*stheta; S2 = stheta + obj.nu*sphi; S4 = sphi^2 + stheta^2 + 2*obj.nu*sphi*stheta + 2*(1-obj.nu)*sigma(3)^2; S5 = 2*(sigmabar^2)*(1-obj.nu)*obj.H(epn)/(9*obj.G) + S4;
                % Dp matrix, according to eqs. 34 of Teng and Rotter (1989a)
                Dp = obj.E/((1-obj.nu^2)*S5)*[S1*S1 , S1*S2 , S1*(1-obj.nu)*sigma(3); S1*S2 , S2*S2 , S2*(1-obj.nu)*sigma(3); S1*(1-obj.nu)*sigma(3) , S2*(1-obj.nu)*sigma(3) , (sigma(3)^2)*((1-obj.nu)^2)];
                % Dep matrix, according to eqs. 34 of Teng and Rotter (1989a), for a point where the material is plastic-hardening (alpha = 1)
                Dep = obj.De - Dp;
            else
                % Dep matrix, according to eqs. 34 of Teng and Rotter (1989a), for a point where the material is elastic (alpha = 0)
                Dep = obj.De;
            end
        end

        % Compute tangent modulus matrix DT according to eq. 40 of Teng and Rotter (1989a)
        function DT = DT(obj,t,sigmas,sigma_ys,epns)

            Z =@(z) [1 0 0 z 0 0; 0 1 0 0 z 0; 0 0 1 0 0 z;]; % Z matrix, according to eq. 12 of Teng and Rotter (1989a)
            zmin = -t/2; zmax = t/2;
            if isempty(sigmas) || isempty(sigma_ys)
                % Only elastic portion of modulus matrix if no sigmas or number of integration stations provided
                zm = zmax - zmin; zm2 = (zmax*zmax - zmin*zmin)*0.5; zm3 = (zmax*zmax*zmax - zmin*zmin*zmin)/3;
                DT = [zm*obj.De zm2*obj.De; zm2*obj.De zm3*obj.De];
            else
                noStations = length(sigma_ys);
                if strcmp(obj.type,'orthotropic'); error('AQUINAS Error: Orthotropic material currently not supported for materially nonlinear type of analysis.'); end
                DT = zeros(6,6);
                % Numerical integration according to Simpson's 1/3 rule
                noLayers = noStations - 1;
                z_stations = zmin:(t/noLayers):zmax; % through thickness coordinate z of integration stations for the tangent modulus matrix DT, starting from the bottom fiber of the shell section and moving to the top
                for I = 1:(noLayers/2)
                    a = 2*I-1; m = 2*I; b = 2*I+1;
                    Z_a = Z(z_stations(a)); Z_m = Z(z_stations(m)); Z_b = Z(z_stations(b));
                    DT_a = Z_a'*obj.Dep(sigmas(:,a),sigma_ys(a),epns(a))*Z_a; % Tangent modulus matrix of bottom fiber a of current integration range
                    DT_m = Z_m'*obj.Dep(sigmas(:,m),sigma_ys(m),epns(m))*Z_m; % Tangent modulus matrix of middle fiber m of current integration range
                    DT_b = Z_b'*obj.Dep(sigmas(:,b),sigma_ys(b),epns(b))*Z_b; % Tangent modulus matrix of top fiber j of current integration range
                    % Numerical integration according to Simpson's 1/3 rule
                    DT = DT + (z_stations(b) - z_stations(a))*(DT_a + 4*DT_m + DT_b)/6;
                end
            end

        end

        % Compute (1-ksi) portion of plastic strain out of the current strain increment, according to eq. 85 of Teng and Rotter (1989a)
        function ksi = elastic_strain_portion(obj,sigma0,Deltasigma,sigmay)
            A = Deltasigma(1)^2 + Deltasigma(2)^2 - Deltasigma(1)*Deltasigma(2) + 3*Deltasigma(3)^2;
            B = Deltasigma(1)*(2*sigma0(1)-sigma0(2)) + Deltasigma(2)*(2*sigma0(2)-sigma0(1)) + 6*sigma0(3)*Deltasigma(3);
            C = sigma0(1)^2 + sigma0(2)^2 - sigma0(1)*sigma0(2) + 3*sigma0(3)^2 - sigmay^2;
            if C > 0; C = 0; end
            ksi = (-B + sqrt(B^2 - 4*A*C))/2/A;
            if ksi < 0; ksi = 0; end
        end

        % Compute plastic strains - sigmas using the sub incremental technique
        function [flag,sigman,Depn,alpha] = sub_incremental_computations(obj,flag,epsilon0,DeltaepsilonZ,sigma0,sigmay,epn,Depn,alpha,epsilon_s,max_epsilon_bar)

            if alpha == 0 || (alpha == 1 && Depn < 0) % Point that is assumed to behave elastically
                Deltasigma = obj.De * DeltaepsilonZ;
                sigmabar = sqrt((sigma0(1)+Deltasigma(1))^2 + (sigma0(2)+Deltasigma(2))^2 - (sigma0(1)+Deltasigma(1))*(sigma0(2)+Deltasigma(2)) + 3*(sigma0(3)+Deltasigma(3))^2); % von Mises stress
                if sigmabar > sigmay*(1-obj.yield_tol) % The point has yielded
                   % (1-ksi) is the portion of the current strain increment that causes plasticity to occur for the material
                    ksi = obj.elastic_strain_portion(sigma0,Deltasigma,sigmay);
                    alpha = 1; % This point is now plastic hardening
                else % The point is still behaving elastically
                    sigman = sigma0 + Deltasigma;
                    return
                end
            elseif alpha == 1 % Point which is already plastic hardening
                ksi = 0; % Since the point is already plastic hardening, (1-ksi) is the portion of it (all of it) that causes plasticity to occur for the material
            end

            sigman = sigma0 + ksi*obj.De*DeltaepsilonZ;
            epsilonn = epsilon0 + ksi*DeltaepsilonZ;
            % Effective strain increment, according to par. 10-(8)-(iv) of Teng and Rotter (1989a)
            Delta_epsilon_bar = (2*(1-ksi)/sqrt(3))*sqrt(DeltaepsilonZ(1)^2 + DeltaepsilonZ(2)^2 + DeltaepsilonZ(1)*DeltaepsilonZ(2) + (DeltaepsilonZ(3)^2)/4);
            % An effective strain increment that is higher than the max_epsilon_bar will not be attempted
            if Delta_epsilon_bar > max_epsilon_bar; flag = -6; end
            % Number of required sub increments, according to par. 10-(8)-(iv) of Teng and Rotter (1989a)
            Nsb = ceil(Delta_epsilon_bar/epsilon_s);
            % Strain increment to be applied per sub-increment, according to par. 10-(8)-(iv) of Teng and Rotter (1989a)
            depsilon_Z = ((1-ksi)/Nsb)*DeltaepsilonZ;
            for n = 1:Nsb
                if alpha == 0 % The point is behaving elastically (unloaded during the n-1 step)
                    [flag,sigman,Depn,alpha] = obj.sub_incremental_computations(flag,epsilonn,depsilon_Z,sigman,obj.sigmay(sigmay,epn,max(0,Depn)),epn,Depn,alpha,epsilon_s,max_epsilon_bar); % recursive call to the function it self, in order to return to step (ii) of methodology of paragraph 10 - (8) of Teng and Rotter 1989a
                    if flag < 0; return; end
                    epsilonn = epsilonn + depsilon_Z;
                elseif alpha == 1 % The point is plastic hardening
                    % von Mises stress, according to eq. 24 of Teng and Rotter (1989a)
                    sigmabar = sqrt(sigman(1)^2 + sigman(2)^2 - sigman(1)*sigman(2) + 3*sigman(3)^2);
                    % Deviatoric sigmas, according to eqs. 30 of Teng and Rotter (1989a)
                    sphi = (2*sigman(1) - sigman(2))/3; stheta = (2*sigman(2) - sigman(1))/3; sphitheta = 2*sigman(3);
                    % S1, S2, S3, S4 and S5 parameters, according to eqs. 35 of Teng and Rotter (1989a)
                    S1 = sphi + obj.nu*stheta; S2 = stheta + obj.nu*sphi; S4 = sphi^2 + stheta^2 + 2*obj.nu*sphi*stheta + 2*(1-obj.nu)*sigman(3)^2; S5 = 2*(sigmabar^2)*(1-obj.nu)*obj.H(epn)/(9*obj.G) + S4;
                    S3 = S1*depsilon_Z(1) + S2*depsilon_Z(2) + (1-obj.nu)*sigman(3)*depsilon_Z(3);
                    % Sub-incremental equivalent plastic strain, according to eq. 36 of of Teng and Rotter (1989a)
                    depn = (2*sigmabar/3)*S3/S5;
                    % Accumulated sub-incremental equivalent plastic strain, according to par. 10-(8)-(v)-(b) of Teng and Rotter (1989a)
                    Depn = Depn + depn;
                    if Depn > 0
                        % Update yield stress for current sub-increment
                        sigmayn = obj.sigmay(sigmay,epn,Depn);
                        % Stresses at half the step, according to eq. 72 of Teng and Rotter (1989a)
                        sigma_n_halfstep = sigman + obj.Dep(sigman,sigmayn,epn)*depsilon_Z;
                        % Current stress increment, according to eq. 73 of Teng and Rotter (1989a)
                        dsigma_n = 0.5*(obj.Dep(sigman,sigmayn,epn) + obj.Dep(sigma_n_halfstep,sigmayn,epn))*depsilon_Z;
                    else
                        alpha = 0;
                        depn_Z = Depn*(3/2/sigmabar)*[sphi stheta sphitheta]';
                        deen_Z = (depsilon_Z - depn_Z);
                        dsigma_n = obj.De*deen_Z;
                    end
                    sigman = sigman + dsigma_n;
                    epsilonn = epsilonn + depsilon_Z;
                end
            end
            % von Mises stress, according to eq. 24 of Teng and Rotter (1989a)
            sigmabar = sqrt(sigman(1)^2 + sigman(2)^2 - sigman(1)*sigman(2) + 3*sigman(3)^2);
            % Current yield stress
            sigmayn = obj.sigmay(sigmay,epn,max(0,Depn));
            % Yield surface function value
            F1 = sigmabar - sigmayn;
            if F1 < -obj.yield_tol*sigmayn; alpha = 0; else; alpha = 1; end
            while F1 > obj.yield_tol*sigmayn
                % Deviatoric sigmas, according to eqs. 30 of Teng and Rotter (1989a)
                sphi = (2*sigman(1) - sigman(2))/3; stheta = (2*sigman(2) - sigman(1))/3; sphitheta = 2*sigman(3);
                % Compute stress corrections, according to eq. 79 of Teng and Rotter (1989a)
                deltasigma = -2*sigmabar*F1/(3*(sphi^2+stheta^2+4*sigman(3)^2))*[sphi stheta sphitheta]';
                % Update sigmas of current sub increment n with the stress corrections
                sigman = sigman + deltasigma;
                % von Mises stress, according to eq. 24 of Teng and Rotter (1989a)
                sigmabar = sqrt(sigman(1)^2 + sigman(2)^2 - sigman(1)*sigman(2) + 3*sigman(3)^2);
                % Yield surface function value
                F1 = sigmabar - sigmayn;
            end

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