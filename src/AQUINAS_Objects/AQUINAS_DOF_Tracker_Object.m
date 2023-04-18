classdef AQUINAS_DOF_Tracker_Object
    % AQUINAS_DOF_Tracker_Object Class Definition for a DOF Tracker object
    % to trace the load path a specific DOF in the global system
    % AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.

    properties
        class_type = 'AQUINAS_DOF_Tracker_Object';
        name % The unique name of the DOF Tracker object, as defined by the User
        ID % The unique ID of the DOF Tracker object, as defined during the analysis
        r % Global r coordinate of the node where the DOF to be tracked lives
        z % Global z coordinate of the node where the DOF to be tracked lives
        active_end % The bottom or top end of the associated segment on which the desired DOF lives
        form % The desired form of the DOF to be plotted in ( 'alg' for algebraic , 'abs' for absolute, 'flip' for flipped (mulitplication with -1))
        DOF_type % The type of the DOF to be tracked ( 'u' , 'w' , 'b'/'beta' )
        segID % The ID of the associated Segment object, as defined during the analysis
        global_DOF % The uncondensed DOF of the global system that corresponds to the DOF to be tracked
        DOF = [] % A vector containing the value of the desired Degree Of Freedom (per step in the case on incremental type of analyses)
        LPF = [] % A vector containing the LPF value that was used to obtain the corresponding entry in the DOF vector property (per step in the case on incremental type of analyses)
        colour % RGB triplet for plotting purposes of the DOF to be tracked
        marker % Marker type for plotting purposes of the DOF to be tracked
        markerSize % Marker size for plotting purposes of the DOF to be tracked
    end

    methods
    % Initiator method for a DOF Tracker Object
        function obj = AQUINAS_DOF_Tracker_Object(options)
            arguments
                options.name { mustBeNonzeroLengthText }
                options.segment { mustBeA(options.segment,"AQUINAS_Segment_Object") }
                options.activeEnd { mustBeNonzeroLengthText }
                options.typeOfDOF { mustBeNonzeroLengthText }
                options.colour { mustBeReal } = [rand rand rand]
                options.marker { mustBeNonzeroLengthText } = 'o'
                options.markerSize { mustBeNonnegative } = 3
                options.form { mustBeNonzeroLengthText } = 'alg'
            end
            if ~isfield(options,"name"); error('AQUINAS Error: DOF Tracker Object definition must include a name for the tracker'); end
            if ~isfield(options,"segment"); error('AQUINAS Error: DOF Tracker Object definition must include the segment that the tracker is applied on'); end
            if ~isfield(options,"activeEnd"); error('AQUINAS Error: DOF Tracker Object definition must include the end of the segment that the tracker is applied on'); end
            if ~isfield(options,"typeOfDOF"); error('AQUINAS Error: DOF Tracker Object definition must include the type of DOF to be tracked'); end
            obj.name = options.name;
            obj.active_end = options.activeEnd;
            obj.DOF_type = options.typeOfDOF;
            obj.colour = options.colour;
            obj.marker = options.marker;
            obj.markerSize = options.markerSize;
            obj.form = options.form;
            if strcmp(obj.active_end,'top')
                obj.r = options.segment.rtop; obj.z = options.segment.ztop;
            elseif strcmp(obj.active_end,'bottom')||strcmp(obj.active_end,'bot')
                obj.r = options.segment.rbot; obj.z = options.segment.zbot;
            end
        end

        % Method to assign a unique global ID to the DOF Tracker Object
        function obj = assign_global_ID(obj,ID)
            obj.ID = ID;
        end
        % Method to associate the DOF Tracker with the Segment object on which the desired DOF lives
        function obj = associate_Segment_Object(obj,SID)
            obj.segID=SID;
        end

        % Method to associate the DOF Tracker with the serial number of a
        % DOF in the global system (after condensation of dvds and 'gamma')
        function obj = global_Tracker(obj,dof)
            global_Node_ID=1+(dof(1)-1)/4;
            switch obj.DOF_type
            case 'u'
                obj.global_DOF=1+4*(global_Node_ID-1);
            case 'w'
                obj.global_DOF=3+4*(global_Node_ID-1);
            case { 'b' , 'beta' }
                obj.global_DOF=4+4*(global_Node_ID-1);
            otherwise
                error(['AQUINAS Error: : Unsupported DOF ', obj.DOF_type ,' for Tracker ', obj.name]);
            end
        end

        % Method to extract the Tracked DOF value from a given vector and store its value into the DOF cell array
        function obj = extract_DOF_LPF_pair(obj,vector,LPF)
            obj.DOF(end+1) = vector(obj.global_DOF);
            obj.LPF(end+1) = LPF;
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