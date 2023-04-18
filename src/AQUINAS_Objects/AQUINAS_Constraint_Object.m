classdef AQUINAS_Constraint_Object
    % Class definition for a Constraint Object - used to apply boundary
    % conditions on nodal variables and other constraints through the use
    % of Lagrange multipliers.
    % At least one Constraint Object must be user-created to prevent rigid body motion
    % The Constraint Object can only refer to a node at a Segment Object boundary
    % AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.

    properties
        class_type = 'AQUINAS_Constraint_Object';
        % Properties of all Constraint Objects
        type % Type of the constraint object, can be either 'bc' for boundary condition (default) value or 'ecc' for eccentricity constraint
        origin % Constraint Object origin ('user' - created by user; 'auto' - created automatically)
        state % The state of the solution process that this constraint should be considered ('prebuckling' - only during the LA - GMNA steps, 'buckling' - only during the bifurcation checks, 'general' - both of the previous states)
        ID % Integer global Constraint Object identifier
        rcoord % The global r coordinate of the boundary constraint, identifying the node where the Constraint Object is applied. Only relevant if 'type' == 'bc'
        zcoord % The global z coordinate of the boundary constraint, identifying the node where the Constraint Object is applied. Only relevant if 'type' == 'bc'
        local_dofs = {} % A cell vector of condensed global nodal dof types (up to 4 per node: u, v, w, b)
        global_dofs = {} % A cell vector, each entry of which is a row vector of condensed global nodal dof IDs
        coeffLHS = {} % A vector of factoring coefficients on each of the dofs in dofIDs
        coeffRHS = [] % The RHS side of the equation built from dofIDs and coeffLHS
        eccEnds = {} % A two entry cell vector holding the ends of the shell segments that are eccentrically connected. Only relevant if type == 'ecc'
        eccSegIDs = [] % A two entry cell that holds the IDs of the segments to be eccentrically connected. Appropriate values are 'top' or 'bot', where each entry refers to the ends of the segments provided in the segments input cell. Only relevant if type == 'ecc'
        eccSegEndCoords = {} % A two entry cell that holds structs with the top and bottom coordinates of the shell segments to be eccentrically connected. Only relevant if type == 'ecc'
        % e.g. to specify a restrained global vertical displacement dof:
        % dofIDs = [1]; coeffLHS = [1]; coeffRHS = 0;
        % e.g. to specify a constraint such that u + v + w = 0 at a node:
        % dofIDs = [1 2 3]; coeffLHS = [1 1 1]; coeffRHS = 0;
        Segment_Object_Props = {} % A cell vector of Segment Object Properties relevant to this Constraint Object. The properties stored here are different depending on whether this is a boundary or eccentricity constraint.
    end

    methods
        % Initiator method for a Constraint Object
        function obj = AQUINAS_Constraint_Object(options)
            arguments
                options.type { mustBeNonzeroLengthText } = 'bc'
                options.origin { mustBeNonzeroLengthText } = 'user'
                options.state { mustBeNonzeroLengthText } = 'general'
                options.rcoord { mustBeNonnegative }
                options.zcoord { mustBeReal }
                options.dofs
                options.coeffLHS { mustBeReal } = 1
                options.coeffRHS { mustBeReal } = 0
                options.segments { mustBeA(options.segments,"cell") }
                options.ends { mustBeA(options.ends,"cell") }
            end
            if strcmp(options.type,'bc')
                if ~isfield(options,'rcoord'); error('AQUINAS Error: Boundary constraint object definition must include the radial coordinate r of the node to be constrained'); end
                if ~isfield(options,'zcoord'); error('AQUINAS Error: Boundary constraint object definition must include the axial coordinate z of the node to be constrained'); end
                if ~isfield(options,'dofs'); error('AQUINAS Error: Boundary constraint object definition must include the cell array of dofs to be constrained'); end
            elseif strcmp(options.type,'ecc')
                if ~isfield(options,'segments'); error('AQUINAS Error: Eccentricity constraint object definition must include which the segments to be connected eccentrically'); end
                if ~isfield(options,'ends'); error('AQUINAS Error: Eccentricity constraint object definition must include which ends of the provided segments are to be connected eccentrically'); end
                if (length(options.segments)~=2); error('AQUINAS Error: Eccentricity constraint object definition must only include the two segments to be connected eccentrically'); end
                if (length(options.ends)~=2); error('AQUINAS Error: Eccentricity constraint object definition must only include the two ends of the segments of which the ends are to be connected eccentrically'); end
                if ~(strcmp(options.ends{1},'bot')||strcmp(options.ends{1},'top')); error('AQUINAS Error: Eccentricity constraint object only supports ''bot'' and ''top'' as options for the ends to be connected'); end
                if ~(strcmp(options.ends{2},'bot')||strcmp(options.ends{2},'top')); error('AQUINAS Error: Eccentricity constraint object only supports ''bot'' and ''top'' as options for the ends to be connected'); end
            else
                error('AQUINAS Error: Constraint object type can be either bc (for boundary condition) or ecc (for eccentricity)')
            end
            obj.type = options.type; obj.origin = options.origin; obj.state = options.state;
            if strcmp(options.type,'bc')
                obj.rcoord = options.rcoord; obj.zcoord = options.zcoord;
                if strcmp(obj.origin,'user'); obj.local_dofs = options.dofs; end % In this case dofs is a cell vector of chars
                if strcmp(obj.origin,'auto'); obj.global_dofs = {options.dofs}; end % In this case dofs is a vector of global nodal dof IDs
                obj.coeffLHS = {options.coeffLHS}; obj.coeffRHS = [options.coeffRHS];
            elseif strcmp(options.type,'ecc')
                obj.eccSegEndCoords = cell(2,1);
                obj.eccSegEndCoords{1} = struct('rbot',options.segments{1}.rbot,'zbot',options.segments{1}.zbot,'rtop',options.segments{1}.rtop,'ztop',options.segments{1}.ztop);
                obj.eccSegEndCoords{2} = struct('rbot',options.segments{2}.rbot,'zbot',options.segments{2}.zbot,'rtop',options.segments{2}.rtop,'ztop',options.segments{2}.ztop);
                obj.eccEnds = options.ends;
                obj.Segment_Object_Props = cell(2,1);
            end
        end

        % Method to assign a unique global ID to Constraint Object
        function obj = assign_global_ID(obj,ID)
            obj.ID = ID;
        end

        % Identify Segment Objects which are affected by this Constraint Object
        function obj = bc_associate_Segment_Object(obj,segProps)
            obj.Segment_Object_Props{end+1} = segProps;
            % Identify global nodal dof IDs relevant to this Constraint Object
            if length(segProps) > 2
                IDs = zeros(size(obj.local_dofs)); segEndDOFs = cell2mat(segProps(3));
                for I = 1:length(obj.local_dofs)
                    switch obj.local_dofs{I}
                    case 'u'  % Global radial nodal displacement
                        IDs(I) = segEndDOFs(1);
                    case 'v'  % Global circumferential nodal displacement
                        IDs(I) = segEndDOFs(2);
                    case 'w'  % Global axial nodal displacement
                        IDs(I) = segEndDOFs(3);
                    case { 'b' , 'beta' }  % Global meridional nodal rotation
                        IDs(I) = segEndDOFs(4);
                    end
                end
                obj.global_dofs = {IDs};
            end
        end

        % Identify segment objects which are affected by this segment object
        function obj = ecc_check_and_associate_Segment_Object(obj,segment)
            tol = segment.geom_tolerance;
            if (abs(segment.rbot-obj.eccSegEndCoords{1}.rbot)<tol) && (abs(segment.zbot-obj.eccSegEndCoords{1}.zbot)<tol) && (abs(segment.rtop-obj.eccSegEndCoords{1}.rtop)<tol) && (abs(segment.ztop-obj.eccSegEndCoords{1}.ztop)<tol)
                obj.eccSegIDs(1) = segment.ID;
            elseif (abs(segment.rbot-obj.eccSegEndCoords{2}.rbot)<tol) && (abs(segment.zbot-obj.eccSegEndCoords{2}.zbot)<tol) && (abs(segment.rtop-obj.eccSegEndCoords{2}.rtop)<tol) && (abs(segment.ztop-obj.eccSegEndCoords{2}.ztop)<tol)
                obj.eccSegIDs(2) = segment.ID;
            end
        end

        % Set Left Hand Side coefficients of eccentricity constraints to be inlcuded in the equilibrium equation, according to the Lagrange multipliers method
        function obj = ecc_set_coeffLHS(obj,segmentI,segmentJ)
            obj.global_dofs = cell(4,1); obj.coeffLHS = cell(4,1); obj.coeffRHS = zeros(4,1);
            if strcmp(obj.eccEnds{1},'bot') && strcmp(obj.eccEnds{2},'bot')
                Sx = segmentJ.rbot- segmentI.rbot; Sy = segmentJ.zbot- segmentI.zbot; dofsI = segmentI.bot_dofIDs; dofsJ = segmentJ.bot_dofIDs;
            elseif strcmp(obj.eccEnds{1},'bot') && strcmp(obj.eccEnds{2},'top')
                Sx = segmentJ.rtop- segmentI.rbot; Sy = segmentJ.ztop- segmentI.zbot; dofsI = segmentI.bot_dofIDs; dofsJ = segmentJ.top_dofIDs;
            elseif strcmp(obj.eccEnds{1},'top') && strcmp(obj.eccEnds{2},'bot')
                Sx = segmentJ.rbot- segmentI.rtop; Sy = segmentJ.zbot- segmentI.ztop; dofsI = segmentI.top_dofIDs; dofsJ = segmentJ.bot_dofIDs;
            elseif strcmp(obj.eccEnds{1},'top') && strcmp(obj.eccEnds{2},'top')
                Sx = segmentJ.rtop- segmentI.rtop; Sy = segmentJ.ztop- segmentI.ztop; dofsI = segmentI.top_dofIDs; dofsJ = segmentJ.top_dofIDs;
            end
            obj.global_dofs{1} = [dofsI(1) dofsI(4) dofsJ(1)];  obj.coeffLHS{1} = [1, Sy, -1];
            obj.global_dofs{2} = [dofsI(2) dofsJ(2)];           obj.coeffLHS{2} = [1, -1];
            obj.global_dofs{3} = [dofsI(3) dofsI(4) dofsJ(3)];  obj.coeffLHS{3} = [1, -Sx, -1];
            obj.global_dofs{4} = [dofsI(4) dofsJ(4)];           obj.coeffLHS{4} = [1, -1];
            for I=1:length(obj.global_dofs)
                [obj.global_dofs{I},sorted_indeces] = sort(obj.global_dofs{I});
                obj.coeffLHS{I} = obj.coeffLHS{I}(sorted_indeces);
            end
        end

        function equal = eq(obj1,obj2)
            equal = true;
            if strcmp(obj1.type,obj2.type) && strcmp(obj1.type,'bc')
                if obj1.rcoord~=obj2.rcoord; equal=false; return; end
                if obj1.zcoord~=obj2.zcoord; equal=false; return; end
                if obj1.coeffLHS{1}~=obj2.coeffLHS{1}; equal=false; return; end
                if obj1.coeffRHS~=obj2.coeffRHS; equal=false; return; end
                if length(obj1.local_dofs)~=length(obj2.local_dofs); equal=false; return; end
                for I=1:length(obj1.local_dofs)
                    if obj1.local_dofs{I}~=obj2.local_dofs{I}; equal=false; return; end
                end
            elseif strcmp(obj1.type,obj2.type) && strcmp(obj1.type,'ecc')
                if strcmp(obj1.eccEnds{1},obj2.eccEnds{1}); equal=false; end
                if strcmp(obj1.eccEnds{2},obj2.eccEnds{2}); equal=false; end
                if obj1.eccSegEndCoords{1}.rbot~=obj2.eccSegEndCoords{1}.rbot; equal=false; end
                if obj1.eccSegEndCoords{1}.zbot~=obj2.eccSegEndCoords{1}.zbot; equal=false; end
                if obj1.eccSegEndCoords{1}.rtop~=obj2.eccSegEndCoords{1}.rtop; equal=false; end
                if obj1.eccSegEndCoords{1}.ztop~=obj2.eccSegEndCoords{1}.ztop; equal=false; end
                if obj1.eccSegEndCoords{2}.rbot~=obj2.eccSegEndCoords{2}.rbot; equal=false; end
                if obj1.eccSegEndCoords{2}.zbot~=obj2.eccSegEndCoords{2}.zbot; equal=false; end
                if obj1.eccSegEndCoords{2}.rtop~=obj2.eccSegEndCoords{2}.rtop; equal=false; end
                if obj1.eccSegEndCoords{2}.ztop~=obj2.eccSegEndCoords{2}.ztop; equal=false; end
            else
                equal = false;
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