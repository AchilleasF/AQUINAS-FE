classdef AQUINAS_Analysis_Object
    % AQUINAS_Analysis_Object An object used for defining the type of analysis to be executed
    % AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.

    properties
        class_type = 'AQUINAS_Analysis_Object';
        type % Analysis type, corresponding to EN 1993-1-6 analysis terminology ( 'LA' , 'LBA' )
        circumferential_modes = [] % A vector that contains the harmonic modes that need to be considered
        no_eigenvalues = [] % A vector containing the number of eigenvalues to be evaluated for each requested harmonic
        normalize_eigs % An option for normalizing the eigenforms. Can be either '-', 'SecondNorm' or 'MaxValue'.
    end

    methods
        % Constructor method for an Analysis Object
        function obj = AQUINAS_Analysis_Object(options)
            arguments
                options.type { mustBeNonzeroLengthText }
                options.circumferentialModes { mustBeInteger,mustBeNonnegative }
                options.noEigenvalues { mustBeInteger, mustBePositive } = 1
                options.normalizeEigs { mustBeNonzeroLengthText } = '-'
            end

            if ~isfield(options,'type'); error('AQUINAS Error: Analysis Object definition must include the type of analysis to be executed'); end
            obj.type = options.type;
            if strcmp(options.type,'LBA') || strcmp(options.type,'GNA') || strcmp(options.type,'GMNA')
                if isfield(options,'circumferentialModes')
                    options.circumferentialModes(options.circumferentialModes < 0) = [];
                    obj.circumferential_modes = unique(options.circumferentialModes);
                end
                obj.no_eigenvalues = options.noEigenvalues;
                obj.normalize_eigs = options.normalizeEigs;
            end
            if strcmp(options.type,'GNA') || strcmp(options.type,'GMNA'); obj.no_eigenvalues = 1; end
        end

        % Method to assign an integer regarding whether or not to include nonlinear strains in the material stiffness matrix, useful for passing data into C++
        function int = includeNLstrains_toInt(obj)
            if strcmp(obj.type,'LBA')
               int = 0;
            elseif strcmp(obj.type,'GNA')
                int = 1;
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