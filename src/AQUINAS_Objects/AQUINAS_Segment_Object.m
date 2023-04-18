classdef AQUINAS_Segment_Object
% Class definition for a Segment Object - used to build the geometry of
% a single axisymmetric shell segment
% At least one Segment Object must exist
% AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.

    properties
        class_type = 'AQUINAS_Segment_Object';
        % Properties of all Segment Objects
        ID % Integer global Segment Object identifier
        type % 'Plate', 'Cone', 'Ellipse', '3pointArc', 'Interp'
        material % Name of Material Object that is to be associated with this Segment Object
        formulation % Thin or thick shell axisymmetric shell formulation ('thin' and 'thick' accordingly)
        thickness % Thickness of the Segment Object

        % A cylinder is a special case of a cone where rtop = rbot
        % A sphere is a special case of an ellipse
        rbot % Global cylindrical radial reference surface r coordinate of base node (LHS node for ring)
        zbot % Global cylindrical vertical reference surface z coordinate of base node (LHS node for ring)
        rtop % Global cylindrical radial reference surface r coordinate of top node (RHS node for ring)
        ztop % Global cylindrical vertical reference surface z coordinate of top node (RHS node for ring)
        bot_connectivity = [] % Vector of Segment Object IDs for base node connectivity (RHS node for ring)
        top_connectivity = [] % Vector of Segment Object IDs for top node connectivity (LHS node for ring)
                              % Empty connectivity vectors mean no connection to another Segment Object (either free node or one where BCs will be enforced)
        geom = [] % Vector of auxiliary data if needed to generate the geometry

        % Equation handles for segment reference surface geometry in terms of parametric coordinate on the Segment Object
        botT % Local parametric coordinate of bottom node (LHS node for ring)
        topT % Local parametric coordinate of top node (LHS node for ring)
        rT % Global r coordinate of reference meridian i.t.o. local parametric T coordinate
        zT % Global z coordinate of reference meridian i.t.o. local parametric T coordinate

        % Properties of Segment Object of type Ellipse
        cenr % Global position of centre of generating ellipse i.t.o. global r coordinate
        cenz % Global position of centre of generating ellipse i.t.o. global z coordinate
        rhor % Semi-radius along global r axis, beginning from cenr in positive global r direction
        rhoz % Semi-radius along global z axis, beginning from cenz in positive global z direction

        % Meshing information
        els % No. of element edges along the arc length of the Segment Object
        g % Grid stretching ratio (g = 0: uniform spacing; 0 < g < 1: towards top of the Segment Object; -1 < g < 0: towards base of Segment Object)
        meshtype % Mesh type: S - single-direction grading; M - graded towards middle; E - graded towards both ends;
        r = [] % Vector of nodal global r coordinates within the Segment Object
        z = [] % Vector of nodal global z coordinates within the Segment Object
        s = [] % Vector of nodal cummulative arc-length coordinate S from the start of the Segment Object
        bot_dofIDs = [] % 1x4 vector of condensed global dof IDs at the base node
        top_dofIDs = [] % 1x4 vector of condensed global dof IDs at the top node

        % Arrays of pre-computed geometric properties at the nodes of a segment
        drdz = [] % dr/dz of reference meridian i.t.o. global z coordinate
        d2rdz2 = [] % d2r/dz2 of reference meridian i.t.o. global z coordinate
        d3rdz3 = [] % d3r/dz3 of reference meridian i.t.o. global z coordinate
        phi = [] % Meridional angle of normal to reference meridian i.t.o. r(z) i.t.o. global z coordinate
        dphids = [] % dphi/ds curvature of reference meridian i.t.o. r(z) i.t.o. global z coordinate
        d2phids2 = [] % d2phi/ds2 aberrancy of reference meridian i.t.o. r(z) i.t.o. global z coordinate
        drds = [] % dr/ds i.t.o. global z coordinate
        dzds = [] % dz/ds i.t.o. global z coordinate
        pp % Enriched piecewise polynomial structure for the definition of the arcs between segment nodes. The breaks and values refer to the input r - z coordinate array. Only relevant for 'Interp' type of segments

        % Constraint Object association
        bot_constraints = [] % A vector of Constraint Object IDs at the base node of the Segment Object
        top_constraints = [] % A vector of Constraint Object IDs at the top node of the Segment Object

        % Nodal Force Object association
        bot_nodal_forces = [] % A vector of Nodal Force Object IDs at the base node of the Segment Object
        mid_nodal_forces = [] % A vector of Nodal Force Object IDs across the surface of the Segment Object
        top_nodal_forces = [] % A vector of Nodal Force Object IDs at the top node of the Segment Object

        % Distributed  Pressure Object association
        distributed_pressures = [] % A vector of Distributed Pressure Objects IDs associated with the Segment Object

        % Miscellaneous
        geom_tolerance = 1e-6 % Tolerance for geometric operations
    end

    methods
        % Initiator method for a Segment Object
        function obj = AQUINAS_Segment_Object(options)
            arguments
                options.type { mustBeNonzeroLengthText }
                options.material { mustBeNonzeroLengthText }
                options.formulation { mustBeNonzeroLengthText } = 'thin'
                options.thickness { mustBePositive }
                options.rzbot { mustBeReal }
                options.rztop { mustBeReal }
                options.rz3rd { mustBeReal }
                options.rzpts { mustBeReal }
                options.geom { mustBeReal } = []
                options.els { mustBeInteger, mustBeNonnegative }
                options.g { mustBeReal } = 0.0
                options.meshtype { mustBeNonzeroLengthText } = 'S'
            end

            if ~isfield(options,'type'); error('AQUINAS Error: Segment Object definition must include the type of the segment.'); end
            if ~isfield(options,'material'); error('AQUINAS Error: Segment Object definition must include a material name for the segment.'); end
            if ~isfield(options,'thickness'); error('AQUINAS Error: Segment Object definition must include the thickness of the segment.'); end
            if ~strcmp(options.type,'Interp') && ~isfield(options,'rzbot'); error('AQUINAS Error: Segment Object definition must include the r - z coordinates of the bottom end of the segment.'); end
            if ~strcmp(options.type,'Interp') && ~isfield(options,'rztop'); error('AQUINAS Error: Segment Object definition must include the r - z coordinates of the top end of the segment.'); end
            if ~isfield(options,'els'); error('AQUINAS Error: Segment Object definition must include number of elements for discretization of segment.'); end
            obj.type = options.type;
            obj.material = options.material;
            obj.formulation = options.formulation;
            obj.thickness = options.thickness;
            obj.geom = options.geom;
            if strcmp(options.type,'Interp')
                if ~isfield(options,'rzpts'); error('AQUINAS Error: Segment Object definition must include the array of r - z coordinates for a meridian to be determined through interpolation.'); end
                if options.rzpts(1,2) < options.rzpts(end,2); warning('AQUINAS Warning : rzpts provided need to be in top - down order. The rzpts array will be flipped'); options.rzpts = flipud(options.rzpts); end
                obj.rbot = options.rzpts(end,1); obj.zbot = options.rzpts(end,2); obj.rtop = options.rzpts(1,1); obj.ztop = options.rzpts(1,2);
            else
                obj.rbot = options.rzbot(1); obj.zbot = options.rzbot(2); obj.rtop = options.rztop(1); obj.ztop = options.rztop(2);
                if obj.ztop < obj.zbot; error('AQUINAS Error: Segment Object''s rzbot point must be below the rztop point.'); end
            end

            % Initial geometry of reference surface meridian
            switch options.type
                case 'Plate'
                    obj.botT = 0;
                    obj.topT = 1;
                    obj.rT = @(T) obj.rbot + (obj.rtop - obj.rbot)*T;
                    obj.zT = @(T) obj.zbot + (obj.ztop - obj.zbot)*T;
                case 'Cone'
                    obj.botT = 0;
                    obj.topT = 1;
                    obj.rT = @(T) obj.rbot + (obj.rtop - obj.rbot)*T;
                    obj.zT = @(T) obj.zbot + (obj.ztop - obj.zbot)*T;
                case 'Ellipse'
                    if isempty(obj.geom); error('AQUINAS Error: Segment Object definition must include geometrical properties of the elliptical segment [cenr cenz rhor rhoz].'); end
                    obj.cenr = obj.geom(1); obj.cenz = obj.geom(2); obj.rhor = obj.geom(3); obj.rhoz = obj.geom(4);
                    obj.botT = atan2(obj.rhor*(obj.zbot - obj.cenz),obj.rhoz*(obj.rbot - obj.cenr));
                    obj.topT = atan2(obj.rhor*(obj.ztop - obj.cenz),obj.rhoz*(obj.rtop - obj.cenr));
                    obj.rT = @(T) obj.cenr + obj.rhor*cos(T);
                    obj.zT = @(T) obj.cenz + obj.rhoz*sin(T);
                case '3pointArc'
                    if ~isfield(options,'rz3rd'); error('AQUINAS Error: Segment Object definition must include the r - z coordinates of 3rd point for 3pointArc segments.'); end
                    % A 3 point arc is basically a circular segment part, which is under the 'Ellipse' segment category.
                    % The '3pointArc' definition only helps to simplify the process of defining such an arc for the analyst.
                    slope1 = - 1/((options.rz3rd(2) - obj.zbot)/(options.rz3rd(1) - obj.rbot));      slope2 = - 1/((obj.ztop - options.rz3rd(2))/(obj.rtop - options.rz3rd(1)));
                    Rm1 = (options.rz3rd(1) + obj.rbot)/2; Zm1 = (options.rz3rd(2) + obj.zbot)/2;  Rm2 = (options.rz3rd(1) + obj.rtop)/2; Zm2 = (options.rz3rd(2) + obj.ztop)/2;
                    obj.geom = zeros(4,1);
                    obj.geom(1) = (Zm2 - Zm1 + slope1*Rm1 - slope2*Rm2)/(slope1 - slope2); % r coordinate of circular arcs centre
                    obj.geom(2) = Zm1 + slope1*(obj.geom(1) - Rm1); % z coordinate of circular arcs centre
                    obj.geom(3) = sqrt(((obj.rbot - obj.geom(1))^2) + ((obj.zbot - obj.geom(2))^2)); % radius of circular arc
                    obj.geom(4) = obj.geom(3);
                    obj.cenr = obj.geom(1); obj.cenz = obj.geom(2); obj.rhor = obj.geom(3); obj.rhoz = obj.geom(4);
                    % Evaluation of top, mid and bottom nodal angles with respect to the centre of the circle that the arc belongs to, in the [-pi,pi] range of angles
                    obj.botT = atan2(obj.rhor*(obj.zbot - obj.cenz),obj.rhoz*(obj.rbot - obj.cenr));
                    midT = atan2(obj.rhor*(options.rz3rd(2) - obj.cenz),obj.rhoz*(options.rz3rd(1) - obj.cenr));
                    obj.topT = atan2(obj.rhor*(obj.ztop - obj.cenz),obj.rhoz*(obj.rtop - obj.cenr));
                    obj.rT = @(T) obj.cenr + obj.rhor*cos(T);
                    obj.zT = @(T) obj.cenz + obj.rhoz*sin(T);
                case 'Interp'
                    if ~isfield(options,'rzpts'); error('AQUINAS Error: Segment Object definition must include array of r - z coordinates for Interp segments.'); end
                    obj.pp = AQUINAS_Segment_Object.akima_interp(options.rzpts(:,2),options.rzpts(:,1));
                    Svector = zeros(length(options.rzpts(:,1))-1,1); obj.pp.arcls = zeros(length(options.rzpts(:,1))-1,1); % vector of arc-lengths for the set of r - z points provided
                    for i = 1:length(Svector)
                        if isnan(obj.pp.coefs(i,1))
                            Svector(i) = abs(obj.pp.values(i+1) - obj.pp.values(i));
                        else
                            c1 = obj.pp.coefs(i,3) - 2*obj.pp.coefs(i,2)*obj.pp.breaks(i) + 3*obj.pp.coefs(i,1)*obj.pp.breaks(i)*obj.pp.breaks(i);
                            c2 = obj.pp.coefs(i,2) - 3*obj.pp.coefs(i,1)*obj.pp.breaks(i);
                            c3 = obj.pp.coefs(i,1);
                            Svector(i) = integral(@(z) (1+(3*c3*z.*z+2*c2*z+c1).^2).^(0.5),obj.pp.breaks(i+1),obj.pp.breaks(i));
                        end
                        if i == 1; obj.pp.arcls(i) = Svector(i); else; obj.pp.arcls(i) = obj.pp.arcls(i-1) + Svector(i); end
                    end
            end

            % Meshing along the arc length of the Segment Object
            obj.els = options.els; obj.g = options.g; obj.meshtype = options.meshtype;
            switch options.meshtype
                case 'S' % Mesh graded in a single direction (or ungraded)
                    t_factor = (1:-1/obj.els:0) + (obj.g/pi)*sin(pi*(1:-1/obj.els:0));
                case 'M' % Mesh graded towards the middle of the Segment Object
                    B = (1:-2/obj.els:0) + (obj.g/pi)*sin(pi*(1:-2/obj.els:0));
                    T = (1:-2/obj.els:0) - (obj.g/pi)*sin(pi*(1:-2/obj.els:0));
                    t_factor = [0.5*(1+T) B(2:end)*0.5];
                case 'E' % Mesh graded towards the ends of the Segment Object
                    B = (1:-2/obj.els:0) - (obj.g/pi)*sin(pi*(1:-2/obj.els:0));
                    T = (1:-2/obj.els:0) + (obj.g/pi)*sin(pi*(1:-2/obj.els:0));
                    t_factor = [0.5*(1+T) B(2:end)*0.5];
            end
            if ~strcmp(obj.meshtype,'S') && mod(obj.els,2); t_factor = [t_factor 0]; end
            if strcmp(obj.type,'3pointArc')
                anglespan = obj.topT - obj.botT; if anglespan < 0; anglespan = anglespan + 2*pi; end % counter-clockwise angle span from bottom to top end
                anglemidspan = midT - obj.botT; if anglemidspan < 0; anglemidspan = anglemidspan + 2*pi; end % counter-clockwise angle span from bottom end to internal 3rd point
                if anglespan < anglemidspan; anglespan = anglespan - 2*pi; end
            else
                anglespan = obj.topT - obj.botT;
            end
            % Re-define 3pointArc segments as elliptical ones, and be treated as such from this point on
            if strcmp(obj.type,'3pointArc'); obj.type = 'Ellipse'; end
            % Compute geometric properties at nodes, based on the type of segment used
            obj.r = zeros(1,obj.els+1); obj.z = zeros(1,obj.els+1); obj.phi = zeros(1,obj.els+1); obj.dphids = zeros(1,obj.els+1); obj.d2phids2 = zeros(1,obj.els+1);
            obj.drdz = zeros(1,obj.els+1); obj.drds = zeros(1,obj.els+1); obj.dzds = zeros(1,obj.els+1); obj.d2rdz2 = zeros(1,obj.els+1); obj.d3rdz3 = zeros(1,obj.els+1);
            obj = obj.precompute_nodal_geometric_properties(t_factor,anglespan);
            % Check for nodes defined beyond the axis of revolution, with some numerical tolerance
            if any(obj.r<-obj.geom_tolerance)
                error(['AQUINAS Error: Violation of geometric axisymmetry by segment with bottom node :  r = ',num2str(obj.rbot), ' | z = ',num2str(obj.zbot),' and top node :  r = ',num2str(obj.rtop), ' | z = ',num2str(obj.ztop)]);
            else
                obj.r = abs(obj.r);
            end
        end


        % Method to compute nodal values of geometric properties, necessary for defining the geometry along the meridian within an element, according to eqs. 1 of Teng and Rotter (1989a)
        function obj = precompute_nodal_geometric_properties(obj,t_factor,anglespan)
            switch obj.type
            case 'Plate'
                obj.r = obj.rT(obj.botT + anglespan*t_factor);
                obj.z = obj.zT(obj.botT + anglespan*t_factor);
                obj.s = abs(obj.rbot-obj.rtop)*(1-t_factor);
                obj.drdz = Inf*ones(size(t_factor));
                obj.d2rdz2 = zeros(size(t_factor));
                obj.d3rdz3 = zeros(size(t_factor));
                if (obj.rbot>obj.rtop); obj.drds = ones(size(t_factor)); else; obj.drds = -ones(size(t_factor)); end
                obj.dzds = zeros(size(t_factor));
                obj.phi = zeros(size(t_factor));
                obj.dphids = zeros(size(t_factor));
                obj.d2phids2 = zeros(size(t_factor));
            case 'Cone'
                obj.r = obj.rT(obj.botT + anglespan*t_factor);
                obj.z = obj.zT(obj.botT + anglespan*t_factor);
                obj.s = sqrt((obj.rtop-obj.rbot)^2 + (obj.ztop-obj.zbot)^2)*(1-t_factor);
                obj.drdz = (obj.rbot - obj.rtop)/(obj.zbot - obj.ztop)*ones(size(obj.r));
                obj.d2rdz2 = zeros(size(t_factor));
                obj.d3rdz3 = zeros(size(t_factor));
                if any(isinf(obj.drdz)); obj.drds = ones(size(t_factor)); else; obj.drds = -obj.drdz./sqrt(1+(obj.drdz.^2)); end
                if any(isinf(obj.drdz)); obj.dzds = zeros(size(t_factor)); else; obj.dzds = -1./sqrt(1+(obj.drdz.^2)); end
                obj.phi = pi/2-atan2(-obj.drdz,ones(size(t_factor)));
                obj.dphids = zeros(size(t_factor));
                obj.d2phids2 = zeros(size(t_factor));
            case 'Ellipse'
                obj.r = obj.rT(obj.botT + anglespan*t_factor);
                obj.z = obj.zT(obj.botT + anglespan*t_factor);
                obj.s = zeros(size(t_factor));
                % Identify the sign of the meridional curvature of the segment, for elliptical segments.
                if all((obj.r-obj.cenr) > -obj.geom_tolerance)
                    mer_curv_sign = 1;
                elseif all((obj.r-obj.cenr) < obj.geom_tolerance)
                    mer_curv_sign = -1;
                else
                    error(['AQUINAS Error: Regions of opposite meridional curvature are not allowed in the same segment. Consider dividing the segment into different parts.',...
                        'Violation found in segment with bottom node :  r = ',num2str(obj.rbot), ' | z = ',num2str(obj.zbot),' and top node :  r = ',num2str(obj.rtop), ' | z = ',num2str(obj.ztop)]);
                end
                obj.drdz = mer_curv_sign*(obj.rhor/obj.rhoz)*(obj.cenz - obj.z).*((obj.rhoz + obj.cenz - obj.z).*(obj.rhoz - obj.cenz + obj.z)).^(-1/2);
                obj.d2rdz2 = -mer_curv_sign*obj.rhor*obj.rhoz*((obj.rhoz + obj.cenz - obj.z).*(obj.rhoz - obj.cenz + obj.z)).^(-3/2);
                obj.d3rdz3 = 3*mer_curv_sign*obj.rhor*obj.rhoz*(obj.cenz - obj.z).*((obj.rhoz + obj.cenz - obj.z).*(obj.rhoz - obj.cenz + obj.z)).^(-5/2);
                for i = 1:length(t_factor)
                    fdRdZ = @(z) mer_curv_sign*(obj.rhor/obj.rhoz)*(obj.cenz - z).*((obj.rhoz + obj.cenz - z).*(obj.rhoz - obj.cenz + z)).^(-1/2);
                    if i > 1; obj.s(i) = integral(@(z) (1+fdRdZ(z).^2).^(0.5),obj.z(i),obj.ztop); end
                    if isinf(obj.drdz(i)); obj.dzds(i) = 0.0; else; obj.dzds(i) = -1/sqrt(1+(obj.drdz(i)^2)); end
                    if isinf(obj.drdz(i)); obj.drds(i) = 1.0; else; obj.drds(i) = -obj.drdz(i)/sqrt(1+(obj.drdz(i)^2)); end
                    obj.phi(i) = pi/2-atan2(-obj.drdz(i),1);
                    if isinf(obj.drdz(i)); obj.dphids(i) = mer_curv_sign*obj.rhoz/obj.rhor/obj.rhor; else; obj.dphids(i) = obj.dzds(i)*(obj.d2rdz2(i)/(1 + obj.drdz(i)^2)); end
                    if isinf(obj.drdz(i)); obj.d2phids2(i) = 0.0; else; obj.d2phids2(i) =(-obj.dzds(i))*((obj.d3rdz3(i)/(sqrt(1+obj.drdz(i)^2)^3)) - 3*(obj.drdz(i)*(obj.d2rdz2(i)^2)/(sqrt(1+obj.drdz(i)^2)^5))); end
                end
            case 'Interp'
                obj.s = (1-t_factor)*obj.pp.arcls(end);
                for i = 1:length(obj.s)
                    % Perform interpolation to find where along the meridian the current node lives
                    for j = 1:length(obj.pp.arcls)
                        if obj.pp.arcls(j) >= obj.s(i); ind = j; break; end
                    end
                    nodei = ind; nodej = ind+1;
                    if any(isnan(obj.pp.coefs(ind,:))) || any(isinf(obj.pp.coefs(ind,:)))
                        obj.r(i) = obj.pp.values(nodei) + (obj.pp.values(nodej) - obj.pp.values(nodei))*(obj.s(i) - obj.pp.arcls(nodei))/(obj.pp.arcls(nodej) - obj.pp.arcls(nodei));
                        obj.z(i) = obj.pp.breaks(nodei) + (obj.pp.breaks(nodej) - obj.pp.breaks(nodei))*(obj.s(i) - obj.pp.arcls(nodei))/(obj.pp.arcls(nodej) - obj.pp.arcls(nodei));
                        obj.phi(i) = 0; obj.dphids(i) = 0; obj.d2phids2(i) = 0;
                        obj.drdz(i) = Inf; obj.d2rdz2(i) = 0; obj.d3rdz3(i) = 0; obj.dzds(i) = 0;
                        if obj.pp.values(nodei) <= obj.pp.values(nodej); obj.drds(i) = 1; else; obj.drds(i) = -1; end
                    else
                        c0 = obj.pp.coefs(ind,4) - obj.pp.coefs(ind,3)*obj.pp.breaks(ind) + obj.pp.coefs(ind,2)*obj.pp.breaks(ind)*obj.pp.breaks(ind) - obj.pp.coefs(ind,1)*obj.pp.breaks(ind)*obj.pp.breaks(ind)*obj.pp.breaks(ind);
                        c1 = obj.pp.coefs(ind,3) - 2*obj.pp.coefs(ind,2)*obj.pp.breaks(ind) + 3*obj.pp.coefs(ind,1)*obj.pp.breaks(ind)*obj.pp.breaks(ind);
                        c2 = obj.pp.coefs(ind,2) - 3*obj.pp.coefs(ind,1)*obj.pp.breaks(ind);
                        c3 = obj.pp.coefs(ind,1);
                        obj.z(i) = fzero(@(x) integral(@(z) (1+(3*c3*z.*z+2*c2*z+c1).^2).^(0.5),obj.pp.breaks(nodej),x) - (obj.pp.arcls(ind) - obj.s(i)),obj.pp.breaks(nodej));
                        obj.r(i) = c3*obj.z(i)*obj.z(i)*obj.z(i) + c2*obj.z(i)*obj.z(i) + c1*obj.z(i) + c0;
                        obj.drdz(i) = 3*c3*obj.z(i)*obj.z(i) + 2*c2*obj.z(i) + c1; obj.d2rdz2(i) = 6*c3*obj.z(i) + 2*c2; obj.d3rdz3(i) = 6*c3;
                        obj.drds(i) = - obj.drdz(i)/sqrt(1+(obj.drdz(i)^2));
                        obj.dzds(i) = - 1/sqrt(1+(obj.drdz(i)^2));
                        obj.phi(i) = pi/2-atan2(-obj.drdz(i),1);
                        obj.dphids(i) = obj.dzds(i)*(obj.d2rdz2(i)/(1 + obj.drdz(i)^2));
                        obj.d2phids2(i) =(-obj.dzds(i))*((obj.d3rdz3(i)/(sqrt(1+obj.drdz(i)^2)^3)) - 3*(obj.drdz(i)*(obj.d2rdz2(i)^2)/(sqrt(1+obj.drdz(i)^2)^5)));
                    end
                end
            end
        end

        % Method to assign a unique global ID to Segment Object
        function obj = assign_global_ID(obj,ID)
            obj.ID = ID;
        end

        % Method to calculate the global condensed dof IDs at the base and top Segment Object nodes
        function [dof_sum,obj] = assign_global_dof_IDs(obj,dof_sum)
            obj.top_dofIDs = dof_sum + (1:4);
            obj.bot_dofIDs = dof_sum + (1:4) + 4*obj.els;
            dof_sum = dof_sum + 4*(obj.els+1);
        end

        % Methods to register bottom / top node connectivity to other Segment Objects
        function obj = assign_bot_connectivity(obj,ID)
            obj.bot_connectivity(end+1) = ID;
        end
        function obj = assign_top_connectivity(obj,ID)
            obj.top_connectivity(end+1) = ID;
        end

        % Methods to register bottom / top node connecitivity to Constraint Objects
        function obj = assign_constraint(obj,CID,node)
            switch node
            case 'top'
                obj.top_constraints(end+1) = CID;
            case 'bot'
                obj.bot_constraints(end+1) = CID;
            end
        end

        % Methods to register bottom / intermediate / top node connecitivity to Nodal Force Objects
        function obj = assign_bot_nodal_force(obj,FID)
            obj.bot_nodal_forces(end+1) = FID;
        end
        function obj = assign_intermediate_nodal_force(obj,FID)
            obj.mid_nodal_forces(end+1) = FID;
        end
        function obj = assign_top_nodal_force(obj,FID)
            obj.top_nodal_forces(end+1) = FID;
        end

        % Method to register connectivity to Distributed Pressure Objects
        function obj = assign_distributed_pressure(obj,PID)
            obj.distributed_pressures(end+1)=PID;
        end

        % Method to assign an integer ID for the Segment Object type, useful for passing data into C++
        function int = type_ToInt(obj)

            switch obj.type
            case 'Cone'
                int = 1;
            case 'Plate'
                int = 2;
            case 'Ellipse'
                int = 3;
            case 'Interp'
                int = 4;
            end

        end

        % Equal operator overloading for Segment Objects
        function equal = eq(obj1,obj2)
            equal = true;
            if ~strcmp(obj1.type,obj2.type); equal=false; return; end
            if ~strcmp(obj1.formulation,obj2.formulation); equal=false; return; end
            if obj1.thickness~=obj2.thickness; equal=false; return; end
            if obj1.rbot~=obj2.rbot; equal=false; return; end
            if obj1.zbot~=obj2.zbot; equal=false; return; end
            if obj1.rtop~=obj2.rtop; equal=false; return; end
            if obj1.ztop~=obj2.ztop; equal=false; return; end
            if length(obj1.geom)~=length(obj2.geom); equal=false; return; end
            if ~all(obj1.geom==obj2.geom); equal=false; return; end
            if length(obj1.r)~=length(obj2.r); equal=false; return; end
            if ~all(abs(obj1.r-obj2.r)<obj1.geom_tolerance); equal=false; return; end
            if ~all(abs(obj1.z-obj2.z)<obj1.geom_tolerance); equal=false; return; end
        end

    end

    methods (Static)

        % Method to perform Akima interpolation and compute cubic polynomial coefficients line segments between nodes, according to Akima H. (1970)
        function pp = akima_interp(x,y)
            if length(x) ~= length(y)
                error('AQUINAS Error : Unequal length of r - z coordinates.')
            end
            if length(x) == 1; error('AQUINAS Error : Cannot define a meridian with a single point.'); end
            if length(x) == 2; error('AQUINAS Error : For straight meridians (interpolation applied on two end-points) use the ''Cone'' option.'); end
            if length(x) == 3; error('AQUINAS Error : For meridians defined through 3 points use the ''3pointArc'' option.'); end
            pp.form = 'pp';
            pp.breaks = nan(1,length(x));
            pp.coefs = nan(length(x)-1,4);
            pp.pieces = length(x);
            pp.order = 4;
            pp.dim = 1;
            pp.values = nan(1,length(y));
            t = nan(length(x),1);
            % Compute slope at each of the nodes
            for i = 1:length(x)
                if i == 1
                    x3 = x(i); x4 = x(i+1); x5 = x(i+2);
                    y3 = y(i); y4 = y(i+1); y5 = y(i+2);
                    x1 = 2*x3 - x5; x2 = x4 + x3 - x5;
                    y2 = (x3-x2)*((y5-y4)/(x5-x4)-2*(y4-y3)/(x4-x3)) + y3; if isnan(y2); y2 = y4 + y3 - y5; end
                    y1 = (x2-x1)*((y4-y3)/(x4-x3)-2*(y3-y2)/(x3-x2)) + y2; if isnan(y1); y1 = y3 + y2 - y4; end
                elseif i == 2
                    x2 = x(i-1); x3 = x(i); x4 = x(i+1); x5 = x(i+2);
                    y2 = y(i-1); y3 = y(i); y4 = y(i+1); y5 = y(i+2);
                    x1 = 2*x2 - x4;
                    y1 = (x2-x1)*((y4-y3)/(x4-x3)-2*(y3-y2)/(x3-x2)) + y2; if isnan(y1); y1 = y3 + y2 - y4; end
                elseif i == length(x)-1
                    x1 = x(i-2); x2 = x(i-1); x3 = x(i); x4 = x(i+1);
                    y1 = y(i-2); y2 = y(i-1); y3 = y(i); y4 = y(i+1);
                    x5 = x4 + x3 - x2;
                    y5 = (x5-x4)*(2*(y4-y3)/(x4-x3)-(y3-y2)/(x3-x2)) + y4; if isnan(y5); y5 = y4 + y3 - y2; end
                elseif i == length(x)
                    x1 = x(i-2); x2 = x(i-1); x3 = x(i);
                    y1 = y(i-2); y2 = y(i-1); y3 = y(i);
                    x4 = x3 + x2 - x1; x5 = x4 + x3 - x2;
                    y4 = (x4-x3)*(2*(y3-y2)/(x3-x2)-(y2-y1)/(x2-x1)) + y3; if isnan(y4); y4 = y3 + y2 - y1; end
                    y5 = (x5-x4)*(2*(y4-y3)/(x4-x3)-(y3-y2)/(x3-x2)) + y4; if isnan(y5); y5 = y4 + y3 - y2; end
                else
                    x1 = x(i-2); x2 = x(i-1); x3 = x(i); x4 = x(i+1); x5 = x(i+2);
                    y1 = y(i-2); y2 = y(i-1); y3 = y(i); y4 = y(i+1); y5 = y(i+2);
                end
                m1 = (y2-y1)/(x2-x1); m2 = (y3-y2)/(x3-x2); m3 = (y4-y3)/(x4-x3); m4 = (y5-y4)/(x5-x4);
                t(i) = (m2*abs(m4-m3)+m3*abs(m2-m1))/(abs(m4-m3)+abs(m2-m1));
                if isinf(t(i)) || isnan(t(i))
                    if m1 == m2 && m3 == m4
                        t(i) = (y4-y2)/(x4-x2);
                    else
                        error('AQUINAS Error : Unhandled exception in segment object during execution of Akima interpolation.');
                    end
                end
            end
            % Compute polynomial coefficients for each of the line segments between nodes
            for i = 1:(length(x)-1)
                p0 = y(i);
                p1 = t(i);
                p2 = (3*(y(i+1)-y(i))/(x(i+1)-x(i))-2*t(i)-t(i+1))/(x(i+1)-x(i));
                p3 = (t(i)+t(i+1)-2*(y(i+1)-y(i))/(x(i+1)-x(i)))/(x(i+1)-x(i))/(x(i+1)-x(i));
                pp.coefs(i,:) = [p3 p2 p1 p0];
                if i == 1; pp.breaks(i) = x(i); end
                pp.breaks(i+1) = x(i+1);
                if i == 1; pp.values(i) = y(i); end
                pp.values(i+1) = y(i+1);
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
% THE USE OF THE MATLAB LANGUAGE DOES NOT IMPLY ENDORSEMENT BY MATHWORKS.s