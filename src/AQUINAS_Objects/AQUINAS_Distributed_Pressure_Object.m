classdef AQUINAS_Distributed_Pressure_Object
    % AQUINAS_Distributed_Pressure_Object Class definition for a Distributed Pressure Object - used to apply a distributed pressuse or traction along the meridian of a segment
    % The Distributed Pressure Object must be associated with a specific
    % Segment. In case the user wishes to create a normal pressure on the shell's surface
    % with type='pn', the pressure is considered positive if pointing outwards
    % (and negative inwards). In case the user wishes to create a tangential
    % traction using type='pt', the traction is considered positive if following
    % the eta direction of Fig.1 in the article by Rotter and Teng (1989a) Vol. 31
    % AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.

    properties
        class_type = 'AQUINAS_Distributed_Pressure_Object';
        type % Distributed Pressure Object type string defining direction with respect to the surface ('pn' for normal or 'pt' for tangential)
        functionHandle % The function handle that defines the distribution of the pressure/traction
        withRespectTo % The variable that the distributed pressure function is expressed with respect to (possible values for this are '-', 'r', 'z', and 'phi')
        integrationStations % The number of integration stations that should be used for the numerical integration of the distributed pressure along the meridian of the segment.
                            % The number of integration stations should be one more than those used for the numerical integration of the stiffness matrix.
        ID % The unique ID of the Distributed Pressure object, as defined during the analysis
        segments = {} % A struct holding the radial and axial coordinates of the top and bottom edges of all the segmetns this pressure is applied to
        segIDs = [] % A vector containing the IDs of all the Segment objects that this Distributed Pressure Object is applied to
    end

    methods
        % Initiator method for a Distributed Pressure Object
        function obj = AQUINAS_Distributed_Pressure_Object(options)
            arguments
                options.segments { mustBeA(options.segments,"cell") }
                options.type { mustBeNonzeroLengthText }
                options.functionHandle % must be a function handle
                options.withRespectTo { mustBeNonzeroLengthText } = '-'
            end
            if ~isfield(options,'segments'); error('AQUINAS Error: Distributed Pressure Object definition must include the segments that the pressure is applied on, provided in the form of a cell array (even with only one entry)'); end
            if ~isfield(options,'type'); error('AQUINAS Error: Distributed Pressure Object definition must include the type (pn or pt) of the distributed pressure'); end
            if ~isfield(options,'functionHandle'); error('AQUINAS Error: Distributed Pressure Object definition must include function handle of the distributed pressure'); end
            if ~isfield(options,'withRespectTo'); error('AQUINAS Error: Distributed Pressure Object definition must include the variable that the distributed pressure function is expressed with respect to'); end
            if ~isa(options.functionHandle,'function_handle'); error('AQUINAS Error: The function argument for a Distributed Pressure Object must be a function handle'); end
            for segment=options.segments
                if ~isa(segment{1},'AQUINAS_Segment_Object'); error('AQUINAS Error: Distributed Pressure Object segments property must only include already defined AQUINAS_Segment_Objects'); end
            end

            for S=1:length(options.segments)
                obj.segments{end+1} = struct('r',options.segments{S}.r,'z',options.segments{S}.z);
            end
            obj.type = options.type;
            obj.functionHandle = options.functionHandle;
            obj.withRespectTo = options.withRespectTo;
        end

        % Method to assign a unique global ID to the Distributed Pressure Object
        function obj = assign_global_ID(obj,ID)
            obj.ID = ID;
        end

        % Method to associate the Distributed Pressure Object with the Segment object under pressure
        function obj = associate_Segment_Object(obj,SID)
            obj.segIDs(end+1) = SID;
        end

        % Method to update properties of Distributed Pressure Object based on the Solver Object provided for the analysis
        function obj = update_on_solver(obj,solver)
            obj.integrationStations = solver.No_Gauss_Stations + 1;
        end

        % Method to generate equivalnet nodal loads for an axisymmetric shell element due to the distributed pressure
        function F = eta_integral_equivalent_nodal_load_vector(obj,dat)

            F = zeros(12,1);

            L = dat.L;
            Gauss = AQUINAS_Solver_Object.Gauss_Legendre_nodes_weights(obj.integrationStations);

            for i = 1:length(Gauss.nodes)

                eta = Gauss.nodes(i); eta2 = eta*eta; eta3 = eta2*eta;
                N01 = (2-3*eta+eta3)*0.25; N11 = L*(1-eta-eta2+eta3)*0.25; N02 = (2+3*eta-eta3)*0.25; N12 = L*(-1-eta+eta2+eta3)*0.25; % Eq. 2 in Rotter & Teng (Vol. 31)
                r = N01*dat.r1 + N02*dat.r2 + N11*dat.drds1 + N12*dat.drds2; % Eqs. 1a & 2 in Rotter & Teng (Vol. 31)
                phi = N01*dat.phi1 + N02*dat.phi2 + N11*dat.dphids1 + N12*dat.dphids2;
                c = cos(phi); s = sin(phi);

                magnitude = @(eta) obj.magnitude_arrayfun([dat.r1,dat.r2],[dat.z1,dat.z2],[dat.phi1,dat.phi2])*[(1-eta)/2 (1+eta)/2]';

                switch obj.type
                case 'pn'
                    F = F + Gauss.weights(i)*L*magnitude(eta)*r*[N01*s N11*s 0 0 N01*c N11*c N02*s N12*s 0 0 N02*c N12*c]';
                case 'pt'
                    F = F + Gauss.weights(i)*L*magnitude(eta)*r*[N01*c N11*c 0 0 -N01*s -N11*s N02*c N12*c 0 0 -N02*s -N12*s]';
                end

            end

        end

        % Method to assign an integer ID for the Distributed Pressure Object type, useful for passing data into C++
        function int = type_ToInt(obj)

            switch obj.type
            case 'pn'
                int = 1;
            case 'pt'
                int = 2;
            end

        end

        % Method to generate an array function for computing the magnitude of the distributed pressure at the nodes of a segment
        function out = magnitude_arrayfun(obj,r,z,phi)

            switch obj.withRespectTo
            case {'-'}
                magnitude = @(r,z,phi) obj.functionHandle();
            case {'R','r'}
                magnitude = @(r,z,phi) obj.functionHandle(r);
            case {'Z','z'}
                magnitude = @(r,z,phi) obj.functionHandle(z);
            case {'PHI','Phi','phi'}
                magnitude = @(r,z,phi) obj.functionHandle(phi);
            end

            out = arrayfun(magnitude,r,z,phi);

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