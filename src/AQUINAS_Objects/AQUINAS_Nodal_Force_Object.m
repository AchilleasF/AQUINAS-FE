classdef AQUINAS_Nodal_Force_Object
    % Class definition for a Nodal Force Object - used to apply 'point' forces and moments
    % The Nodal Force Object can only refer to a node at a Segment Object boundary
    % AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.

    properties
        class_type = 'AQUINAS_Nodal_Force_Object';
        % Properties of all Nodal Force Objects
        type % Nodal Force Object type string defining direction ('Fu','Fw' or 'M')
        ID % Integer global Nodal Force Object identifier
        rcoord % Global r coordinates identifying each node relevant to the Nodal Force Object
        zcoord % Global z coordinates identifying each node relevant to the Nodal Force Object
        magnitude % Magnitude of force [F] or moment [F.L], treated as positive if acting in the positive direction of the condensed dof
        origin % the origin of the generation of the Nodal Force Object ('user','auto')
        Segment_Object_IDs = {} % A cell vector of Segment Object IDs relevant to this Nodal Force Object
        % Note that the magnitude will be split equally between all Segment Objects which share this node
        tol_zero_radial_coord = 1e-4 % Tolerance to check whether a nodal force is applied on the axis of revolution
    end

    methods
        % Initiator method for a Nodal Forcel Object
        function obj = AQUINAS_Nodal_Force_Object(options)
            arguments
                options.type { mustBeNonzeroLengthText }
                options.rcoord { mustBeNonnegative }
                options.zcoord { mustBeReal }
                options.magnitude { mustBeReal }
                options.origin { mustBeNonzeroLengthText } = 'user'
            end
            if ~isfield(options,'type'); error('AQUINAS Error: Nodal Force Object definition must include the type of the load to ba applied'); end
            if ~isfield(options,'rcoord'); error('AQUINAS Error: Nodal Force Object definition must include the radial coordinate r of the node where the load is applied'); end
            if ~isfield(options,'zcoord'); error('AQUINAS Error: Nodal Force Object definition must include the axial coordinate z of the node where the load is applied'); end
            if ~isfield(options,'magnitude'); error('AQUINAS Error: Nodal Force Object definition must include the magnitude of the load to be applied'); end

            obj.type = options.type; obj.rcoord = options.rcoord; obj.zcoord = options.zcoord;
            if (strcmp(options.origin,'user'))
                if options.rcoord < obj.tol_zero_radial_coord
                    obj.magnitude = options.magnitude;
                else
                    obj.magnitude = options.rcoord*options.magnitude;
                end
            elseif (strcmp(options.origin,'auto'))
                obj.magnitude = options.magnitude;
            end
        end

        % Method to assign a unique global ID to Nodal Force Object
        function obj = assign_global_ID(obj,ID)
            obj.ID = ID;
        end

        % Identify Segment Objects which are affected by this Nodal Force Object,
        % as well as the global condensed dof ID relevant to the force
        function obj = associate_Segment_Object(obj,SID)
            switch obj.type
            case 'Fu'
                dof = SID{3}(1);
            case 'Fw'
                dof = SID{3}(3);
            case 'M'
                dof = SID{3}(4);
            otherwise
                error(['AQUINAS Error: Unsupported nodal force type:',char(32),obj.type,char(32),'applied at coordinates: r =',char(32),num2str(obj.rcoord),char(32),'| z =',char(32),num2str(obj.zcoord)]);
            end
            obj.Segment_Object_IDs{end+1} = {SID(1),SID(2),dof};
        end

        % Reduce magnitude of Nodal Force Object to account for it being shared across N Segment Objects
        function obj = split_magnitude(obj,N)
            obj.magnitude = obj.magnitude/N;
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