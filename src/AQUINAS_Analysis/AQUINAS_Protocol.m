function [varargout] = AQUINAS_Protocol(varargin)
% AQUINAS - FE for shells of revolution
% Master pre-processing function
% AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PRE-PROCESSOR %%%%
%%%%%%%%%%%%%%%%%%%%%%%%
% Gather input objects into vector or dictionary depositories
SEGMENTS = AQUINAS_Segment_Object.empty;
CONSTRAINTS = AQUINAS_Constraint_Object.empty;
NODAL_FORCES = AQUINAS_Nodal_Force_Object.empty;
DISTRIBUTED_PRESSURES = AQUINAS_Distributed_Pressure_Object.empty;
DOF_TRACKERS = containers.Map;
OUTPUT_MODES = AQUINAS_Output_Mode_Object.empty;
MATERIALS = containers.Map;
for V = 1:length(varargin)
    if strcmp(varargin{V}.class_type,'AQUINAS_Analysis_Object'); ANALYSIS = varargin{V};    end
    if strcmp(varargin{V}.class_type,'AQUINAS_Constraint_Object'); CONSTRAINTS(end+1) = varargin{V}; end
    if strcmp(varargin{V}.class_type,'AQUINAS_Distributed_Pressure_Object'); DISTRIBUTED_PRESSURES(end+1) = varargin{V}; end
    if strcmp(varargin{V}.class_type,'AQUINAS_Material_Object'); MATERIALS(varargin{V}.name) = varargin{V}; end
    if strcmp(varargin{V}.class_type,'AQUINAS_Nodal_Force_Object'); NODAL_FORCES(end+1) = varargin{V}; end
    if strcmp(varargin{V}.class_type,'AQUINAS_Output_Mode_Object'); OUTPUT_MODES(end+1) = varargin{V}; end
    if strcmp(varargin{V}.class_type,'AQUINAS_Segment_Object'); SEGMENTS(end+1) = varargin{V}; end
    if strcmp(varargin{V}.class_type,'AQUINAS_Solver_Object'); SOLVER = varargin{V}; end
    if strcmp(varargin{V}.class_type,'AQUINAS_DOF_Tracker_Object')
        if DOF_TRACKERS.isKey(varargin{V}.name)
            error(['AQUINAS Error: Duplicate DOF Tracker Object : ',varargin{V}.name]);
        else
            DOF_TRACKERS(varargin{V}.name) = varargin{V};
        end
    end
end
clear varargin;

% Order the Segment Objects by decreasing global z position of its z-midpoint
zmid = zeros(1,length(SEGMENTS));
for I = 1:length(SEGMENTS)
    zmid(I) = 0.5*(SEGMENTS(I).zbot + SEGMENTS(I).ztop);
end
[~,I] = sort(-zmid); SEGMENTS = SEGMENTS(I); tol = SEGMENTS(1).geom_tolerance;
% Enumerate dofs at the ends of each segment
dof_sum = 0;
for I = 1:length(SEGMENTS)
    SEGMENTS(I) = SEGMENTS(I).assign_global_ID(I); % Assign unique integer global ID to each Segment Object
    [dof_sum,SEGMENTS(I)] = SEGMENTS(I).assign_global_dof_IDs(dof_sum); % Assign unique integer global IDs to base / top condensed dofs
end
% Detect Segment Object connectivity
seg_connectivity = {};
for I = 1:length(SEGMENTS)
    for J = I+1:length(SEGMENTS)
        if (abs(SEGMENTS(I).rtop - SEGMENTS(J).rtop) < tol) && (abs(SEGMENTS(I).ztop - SEGMENTS(J).ztop) < tol) % top to top connection
            SEGMENTS(I) = SEGMENTS(I).assign_top_connectivity(J); SEGMENTS(J) = SEGMENTS(J).assign_top_connectivity(I);
            seg_connectivity{end+1} = struct('segments',[I J],'endI','top','endJ','top','rcoord',SEGMENTS(I).rtop,'zcoord',SEGMENTS(I).ztop,'dofsI',SEGMENTS(I).top_dofIDs,'dofsJ',SEGMENTS(J).top_dofIDs);
        end
        if (abs(SEGMENTS(I).rtop - SEGMENTS(J).rbot) < tol) && (abs(SEGMENTS(I).ztop - SEGMENTS(J).zbot) < tol) % top to bottom connection
            SEGMENTS(I) = SEGMENTS(I).assign_top_connectivity(J); SEGMENTS(J) = SEGMENTS(J).assign_bot_connectivity(I);
            seg_connectivity{end+1} = struct('segments',[I J],'endI','top','endJ','bot','rcoord',SEGMENTS(I).rtop,'zcoord',SEGMENTS(I).ztop,'dofsI',SEGMENTS(I).top_dofIDs,'dofsJ',SEGMENTS(J).bot_dofIDs);
        end
        if (abs(SEGMENTS(I).rbot - SEGMENTS(J).rtop) < tol) && (abs(SEGMENTS(I).zbot - SEGMENTS(J).ztop) < tol) % bottom to top connection
            SEGMENTS(I) = SEGMENTS(I).assign_bot_connectivity(J); SEGMENTS(J) = SEGMENTS(J).assign_top_connectivity(I);
            seg_connectivity{end+1} = struct('segments',[I J],'endI','bot','endJ','top','rcoord',SEGMENTS(J).rbot,'zcoord',SEGMENTS(J).zbot,'dofsI',SEGMENTS(I).bot_dofIDs,'dofsJ',SEGMENTS(J).top_dofIDs);
        end
        if (abs(SEGMENTS(I).rbot - SEGMENTS(J).rbot) < tol) && (abs(SEGMENTS(I).zbot - SEGMENTS(J).zbot) < tol) % bottom to bottom connection
            SEGMENTS(I) = SEGMENTS(I).assign_bot_connectivity(J); SEGMENTS(J) = SEGMENTS(J).assign_bot_connectivity(I);
            seg_connectivity{end+1} = struct('segments',[I J],'endI','bot','endJ','bot','rcoord',SEGMENTS(J).rbot,'zcoord',SEGMENTS(J).zbot,'dofsI',SEGMENTS(I).bot_dofIDs,'dofsJ',SEGMENTS(J).bot_dofIDs);
        end
    end
end
%Check for identical Segment Objects
for S1=1:length(SEGMENTS)
    for S2=1:length(SEGMENTS)
        if S1==S2; continue; end
        if SEGMENTS(S1)==SEGMENTS(S2); error(['AQUINAS Error: Identical Segments: ',num2str(SEGMENTS(S1).ID),' - ',num2str(SEGMENTS(S2).ID)]); end
    end
end
%Check for identical Constraint Objects
for C1=1:length(CONSTRAINTS)
    for C2=1:length(CONSTRAINTS)
        if C1==C2; continue; end
        if CONSTRAINTS(C1)==CONSTRAINTS(C2)
            if strcmp(CONSTRAINTS(C1).type,'bc')
                error(['AQUINAS Error: Identical Boundary Constraints at coordinates: r = ',num2str(CONSTRAINTS(C1).rcoord),' | z  = ',num2str(CONSTRAINTS(C1).zcoord),' and r = ',num2str(CONSTRAINTS(C2).rcoord),' | z  = ',num2str(CONSTRAINTS(C2).zcoord)]);
            elseif strcmp(CONSTRAINTS(C1).type,'ecc')
                error('AQUINAS Error: Identical Eccentricity Constraints found');
            end
        end
    end
end
% Gaussian Process Regression with pre-defined range of circumferential wave numbers to be checked
if (strcmp(ANALYSIS.type,'LBA')||strcmp(ANALYSIS.type,'GNA')) && ~isempty(ANALYSIS.circumferential_modes) && SOLVER.surrogate_optimisation && SOLVER.so_Auto_Bounds
    warning('Circumferential wavenumbers provided through AQUINAS_Analysis_Object will be neglected and overwritten due to the automatic detection of the critical circumferential wavenumber. Turn so_Auto_Bounds in the AQUINAS_Solver_Object to turn the automatic detection off.');
end
% Uniquely number existing objects
for I = 1:length(CONSTRAINTS)
    CONSTRAINTS(I) = CONSTRAINTS(I).assign_global_ID(I); % Assign unique integer global ID to each user-created Constraint Object
end
for I = 1:length(NODAL_FORCES)
    NODAL_FORCES(I) = NODAL_FORCES(I).assign_global_ID(I); % Assign unique integer global ID to each user-created Nodal Force Object
end
for I = 1:length(DISTRIBUTED_PRESSURES)
    DISTRIBUTED_PRESSURES(I) = DISTRIBUTED_PRESSURES(I).update_on_solver(SOLVER);
    DISTRIBUTED_PRESSURES(I) = DISTRIBUTED_PRESSURES(I).assign_global_ID(I); % Assign unique integer global ID to each user-created Distributed Pressure Object
end
DOF_ID = 0;
for Ikey = DOF_TRACKERS.keys
    DOF_ID = DOF_ID + 1;
    DOF_TRACKERS(Ikey{1}) = DOF_TRACKERS(Ikey{1}).assign_global_ID(DOF_ID); % Assign unique integer global ID to each user-created DOF Tracker Object
end
% Associate various objects with Segment Objects
for I = 1:length(SEGMENTS)
    % Associate Constraint Objects
    for J = 1:length(CONSTRAINTS)
        if strcmp(CONSTRAINTS(J).type,'bc')
            if (abs(SEGMENTS(I).rbot - CONSTRAINTS(J).rcoord) < tol) && (abs(SEGMENTS(I).zbot - CONSTRAINTS(J).zcoord) < tol)
                SEGMENTS(I) = SEGMENTS(I).assign_constraint(J,'bot');
                CONSTRAINTS(J) = CONSTRAINTS(J).bc_associate_Segment_Object({SEGMENTS(I).ID,'B',SEGMENTS(I).bot_dofIDs});
            end
            if (abs(SEGMENTS(I).rtop - CONSTRAINTS(J).rcoord) < tol) && (abs(SEGMENTS(I).ztop - CONSTRAINTS(J).zcoord) < tol)
                SEGMENTS(I) = SEGMENTS(I).assign_constraint(J,'top');
                CONSTRAINTS(J) = CONSTRAINTS(J).bc_associate_Segment_Object({SEGMENTS(I).ID,'T',SEGMENTS(I).top_dofIDs});
            end
        end
    end
    % Associate Nodal Force Objects
    for J = 1:length(NODAL_FORCES)
        if (abs(SEGMENTS(I).rbot - NODAL_FORCES(J).rcoord) < tol) && (abs(SEGMENTS(I).zbot - NODAL_FORCES(J).zcoord) < tol)
            SEGMENTS(I) = SEGMENTS(I).assign_bot_nodal_force(J);
            NODAL_FORCES(J) = NODAL_FORCES(J).associate_Segment_Object({SEGMENTS(I).ID,'B',SEGMENTS(I).bot_dofIDs});
        end
        if (abs(SEGMENTS(I).rtop - NODAL_FORCES(J).rcoord) < tol) && (abs(SEGMENTS(I).ztop - NODAL_FORCES(J).zcoord) < tol)
            SEGMENTS(I) = SEGMENTS(I).assign_top_nodal_force(J);
            NODAL_FORCES(J) = NODAL_FORCES(J).associate_Segment_Object({SEGMENTS(I).ID,'T',SEGMENTS(I).top_dofIDs});
        end
    end
    % Associate DOF_Tracker Objects
    for Jkey = DOF_TRACKERS.keys
        if strcmp(DOF_TRACKERS(Jkey{1}).active_end,'top')
             if (abs(SEGMENTS(I).rtop - DOF_TRACKERS(Jkey{1}).r) < tol) && (abs(SEGMENTS(I).ztop - DOF_TRACKERS(Jkey{1}).z) < tol)
                DOF_TRACKERS(Jkey{1}) = DOF_TRACKERS(Jkey{1}).associate_Segment_Object(SEGMENTS(I).ID);
                DOF_TRACKERS(Jkey{1}) = DOF_TRACKERS(Jkey{1}).global_Tracker(SEGMENTS(I).top_dofIDs);
            end
        elseif strcmp(DOF_TRACKERS(Jkey{1}).active_end,'bottom')||strcmp(DOF_TRACKERS(Jkey{1}).active_end,'bot')
             if (abs(SEGMENTS(I).rbot - DOF_TRACKERS(Jkey{1}).r) < tol) && (abs(SEGMENTS(I).zbot - DOF_TRACKERS(Jkey{1}).z) < tol)
                DOF_TRACKERS(Jkey{1}) = DOF_TRACKERS(Jkey{1}).associate_Segment_Object(SEGMENTS(I).ID);
                DOF_TRACKERS(Jkey{1}) = DOF_TRACKERS(Jkey{1}).global_Tracker(SEGMENTS(I).bot_dofIDs);
             end
        end
    end
    % Associate Distributed Pressure Objects
    for J = 1:length(DISTRIBUTED_PRESSURES)
        for K = 1:length(DISTRIBUTED_PRESSURES(J).segments)
            if all(abs(SEGMENTS(I).r - DISTRIBUTED_PRESSURES(J).segments{K}.r) < tol) && all(abs(SEGMENTS(I).z - DISTRIBUTED_PRESSURES(J).segments{K}.z) < tol)
                SEGMENTS(I) = SEGMENTS(I).assign_distributed_pressure(J);
                DISTRIBUTED_PRESSURES(J) = DISTRIBUTED_PRESSURES(J).associate_Segment_Object(SEGMENTS(I).ID);
            end
        end
    end
end
% Account for sharing of Nodal Force Object among multiple Segment Objects, and also check for unconnected Nodal Force objects
for I = 1:length(NODAL_FORCES)
    if isempty(NODAL_FORCES(I).Segment_Object_IDs); error(['AQUINAS Error: Unconnected Nodal Force at coordinates: r = ',num2str(NODAL_FORCES(I).rcoord),' | z  = ',num2str(NODAL_FORCES(I).zcoord)]); end
    NODAL_FORCES(I) = NODAL_FORCES(I).split_magnitude(length(NODAL_FORCES(I).Segment_Object_IDs));
end
% Check for unconnected Constraint objects
for I = 1:length(CONSTRAINTS)
    if strcmp(CONSTRAINTS(I).type,'bc') && isempty(CONSTRAINTS(I).Segment_Object_Props); error(['AQUINAS Error: Unconnected Constraint at coordinates: r = ',num2str(CONSTRAINTS(I).rcoord),' | z  = ',num2str(CONSTRAINTS(I).zcoord)]); end
end
% Create new Constraint Objects to ensure Segment Object nodal connectivity
for C = 1:length(seg_connectivity)
    for K = 1:4
        CONSTRAINTS(end+1) = AQUINAS_Constraint_Object('rcoord',seg_connectivity{C}.rcoord,'zcoord',seg_connectivity{C}.zcoord,'dofs',[seg_connectivity{C}.dofsI(K) seg_connectivity{C}.dofsJ(K)],'coeffLHS',[1 -1],'coeffRHS',0,'origin','auto');
        CONSTRAINTS(end) = CONSTRAINTS(end).assign_global_ID(length(CONSTRAINTS));
        CONSTRAINTS(end) = CONSTRAINTS(end).bc_associate_Segment_Object({SEGMENTS(seg_connectivity{C}.segments(1)).ID,seg_connectivity{C}.endI});
        CONSTRAINTS(end) = CONSTRAINTS(end).bc_associate_Segment_Object({SEGMENTS(seg_connectivity{C}.segments(2)).ID,seg_connectivity{C}.endJ});
        SEGMENTS(seg_connectivity{C}.segments(1)) = SEGMENTS(seg_connectivity{C}.segments(1)).assign_constraint(length(CONSTRAINTS),seg_connectivity{C}.endI);
        SEGMENTS(seg_connectivity{C}.segments(2)) = SEGMENTS(seg_connectivity{C}.segments(2)).assign_constraint(length(CONSTRAINTS),seg_connectivity{C}.endJ);
    end
end
% Process eccentricity constraint objects
for I = 1:length(CONSTRAINTS)
    if strcmp(CONSTRAINTS(I).type,'ecc')
        for S = 1:length(SEGMENTS)
            CONSTRAINTS(I) = CONSTRAINTS(I).ecc_check_and_associate_Segment_Object(SEGMENTS(S));
        end
        CONSTRAINTS(I) = CONSTRAINTS(I).ecc_set_coeffLHS(SEGMENTS(CONSTRAINTS(I).eccSegIDs(1)),SEGMENTS(CONSTRAINTS(I).eccSegIDs(2)));
    end
end
% Unconnected segments error check
unconnected = []; % Check for unconnected Segment Objects
for S = 1:length(SEGMENTS)
    if isempty(SEGMENTS(S).top_connectivity) && isempty(SEGMENTS(S).bot_connectivity) && length(SEGMENTS) > 1;
        eccentrically_connected = false;
        for C = 1:length(CONSTRAINTS)
            if strcmp(CONSTRAINTS(C).type,'ecc') && (CONSTRAINTS(C).eccSegIDs(1)==SEGMENTS(S).ID || CONSTRAINTS(C).eccSegIDs(2)==SEGMENTS(S).ID)
                eccentrically_connected = true;
            end
        end
        if ~eccentrically_connected; unconnected(end+1) = SEGMENTS(S).ID; end
    end
end
if ~isempty(unconnected); error(['AQUINAS Error: Unconnected segments: ',num2str(unconnected)]); end
% Check that material plastic properties have been provided as input for materially nonlinear types of analysis
if strcmp(ANALYSIS.type,'MNA') || strcmp(ANALYSIS.type,'GMNA')
    for S = 1:length(SEGMENTS)
        if ~strcmp(MATERIALS(SEGMENTS(S).material).curve_def,'True') && ~strcmp(MATERIALS(SEGMENTS(S).material).curve_def,'Engineering')
            error('AQUINAS Error: The stress-strain curve provided for a material used in an MNA or GMNA needs to be explicitly defined as either ''True'' or ''Engineering''');
        end
        if isempty(MATERIALS(SEGMENTS(S).material).sy) || isempty(MATERIALS(SEGMENTS(S).material).ep)
            error('AQUINAS Error: The values of the stress-strain yield curve must be provided as input for materials participating in a materially nonlinear type of analysis.');
        end
    end
end
% Check that a valid DOF Tracker name has been provided in the AQUINAS_Solver_Object in order to check for termination according to the CIP methodology
if (strcmp(ANALYSIS.type,'GNA') || strcmp(ANALYSIS.type,'MNA') || strcmp(ANALYSIS.type,'GMNA')) && any(ismember(SOLVER.termination_conditions,'E'))
    if strcmp(SOLVER.cip.dof_tracker,'-') || ~isKey(DOF_TRACKERS,SOLVER.cip.dof_tracker)
        error('AQUINAS Error: A valid DOF Tracker name needs to be provided as input in the Solver object in order to check for termination of a nonlinear analysis according to the CIP algorithm.');
    end
    SOLVER.cip = struct('dof_tracker',SOLVER.cip.dof_tracker,'tol_cv',0.1,'tol_cip_rpl',1e-5,'counter',0,'prev_cip_rpl',0,'all_omega_bar',[],'all_pms',[]);
end

%%%%%%%%%%%%%%%%%%%
%%%%% ANALYSIS %%%%
%%%%%%%%%%%%%%%%%%%
GUI_call = false;
% LA or pre-LBA/GNA LA step
if strcmp(ANALYSIS.type,'LA') || strcmp(ANALYSIS.type,'LBA')
    [ LA_Delta, LA_Vector, LA_KStorage, LA_FStorage, LA_OStorage ] = AQUINAS_LA(SEGMENTS,CONSTRAINTS,NODAL_FORCES,DISTRIBUTED_PRESSURES,SOLVER,MATERIALS);
    tpoststart = tic;
    [ LA_Data ] = AQUINAS_post_process(LA_Delta,LA_Vector,LA_KStorage,LA_FStorage,LA_OStorage,SEGMENTS,MATERIALS,ANALYSIS,SOLVER,[]);
    tpost = toc(tpoststart);
    if SOLVER.console_output; disp(['    Post Processing time: ',num2str(tpost),' [s]']); disp(' '); end
    if strcmp(ANALYSIS.type,'LA')
        for D=DOF_TRACKERS.keys; DOF_TRACKERS(D{1}) = DOF_TRACKERS(D{1}).extract_DOF_LPF_pair(LA_Delta,1.0); end
        for I=1:length(OUTPUT_MODES)
            if (strcmp(OUTPUT_MODES(I).mode,'WriteFile'))
                OUTPUT_MODES(I).OpenWriteClose(LA_Data,ANALYSIS(1)); % If Output Mode object is WriteFile then write the results in a file
            elseif (strcmp(OUTPUT_MODES(I).mode,'VarArgOut'))
                varargout = {AQUINAS_Output_Mode_Object.LA_structured_output(ANALYSIS, DOF_TRACKERS, SEGMENTS, LA_Data, LA_OStorage)};  % If Output Mode object is VarArgOut then fill a struct with the results
            end
        end
    end
end


% LBA step
if strcmp(ANALYSIS.type,'LBA')
    if (SOLVER.surrogate_optimisation)
        [ LBA_Eigvectors, LBA_Eigvalues, ANALYSIS ] = AQUINAS_SO(SEGMENTS,CONSTRAINTS,SOLVER,MATERIALS,ANALYSIS,[],LA_Data{1},LA_Data{2},LA_OStorage);
    else
        [ LBA_Eigvectors, LBA_Eigvalues ] = AQUINAS_LBA(SEGMENTS,CONSTRAINTS,SOLVER,MATERIALS,ANALYSIS,[],LA_Data{1},LA_Data{2},LA_OStorage);
    end
    LBA_Data = {real(LBA_Eigvalues),real(LBA_Eigvectors)};
    for I=1:length(OUTPUT_MODES)
        if (strcmp(OUTPUT_MODES(I).mode,'WriteFile'))
            OUTPUT_MODES(I).OpenWriteClose(LBA_Data,ANALYSIS(1)); % If Output Mode object is WriteFile then write the results in a file
        elseif (strcmp(OUTPUT_MODES(I).mode,'VarArgOut'))
            varargout = {AQUINAS_Output_Mode_Object.LBA_structured_output(ANALYSIS,SEGMENTS,LBA_Data,LA_OStorage)};
        end
    end
end


% GMNA step
if strcmp(ANALYSIS.type,'GNA') || strcmp(ANALYSIS.type,'MNA') || strcmp(ANALYSIS.type,'GMNA')
    [ GMNA_Data, OStorage ] = AQUINAS_GMNA(MATERIALS,SEGMENTS,CONSTRAINTS,NODAL_FORCES,DISTRIBUTED_PRESSURES,DOF_TRACKERS,SOLVER,ANALYSIS,OUTPUT_MODES);
    [ GMNA_Data ] = AQUINAS_post_process([],[],[],[],OStorage,SEGMENTS,MATERIALS,ANALYSIS,SOLVER,GMNA_Data);
    for I=1:length(OUTPUT_MODES)
        if (strcmp(OUTPUT_MODES(I).mode,'WriteFile'))
            OUTPUT_MODES(I).OpenWriteClose(GMNA_Data,ANALYSIS(1)); % If Output Mode object is WriteFile then write the results in a file
        elseif (strcmp(OUTPUT_MODES(I).mode,'VarArgOut'))
            varargout = {AQUINAS_Output_Mode_Object.GMNA_structured_output(OUTPUT_MODES(I).level,ANALYSIS,SOLVER,SEGMENTS,DOF_TRACKERS,GMNA_Data,OStorage)};
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