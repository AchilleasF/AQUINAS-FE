function [MTV,VCV,offsets] = AQUINAS_initialise_offsets(seg_els,dof_sum,Kc_siz,Kc_hsiz)
% Function to initialise the sparse matrix storage array (MTV), vector
% storage array (VCV) and offset dictionary (offsets) for a particular segment object.
% AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.

MTV = zeros(seg_els * Kc_siz * Kc_siz - (seg_els - 1) * Kc_hsiz * Kc_hsiz, 3);
VCV = zeros((seg_els + 1) * Kc_hsiz, 2);

Ldiag_0 = (seg_els + 1) * Kc_hsiz; % Length of main diagonal (= number of element nodes * dofs per node)
Ldiag_1 = Ldiag_0 - 1; % Length of first upper / lower off-diagonal
Ldiag_2 = Ldiag_1 - 1; % Length of second upper / lower off-diagonal
Ldiag_3 = Ldiag_2 - 1; % Length of third upper / lower off-diagonal
Ldiag_4 = Ldiag_3 - 1; % Length of fourth upper / lower off-diagonal
Ldiag_5 = 3 * seg_els; % Length of fifth upper / lower off-diagonal
Ldiag_6 = 2 * seg_els; % Length of sixth upper / lower off-diagonal
Ldiag_7 = seg_els; % Length of seventh upper / lower off-diagonal

offsets = containers.Map();
offsets('D0') = 0;
offsets('U1') = Ldiag_0;
offsets('L1') = offsets('U1') + Ldiag_1;
offsets('U2') = offsets('L1') + Ldiag_1;
offsets('L2') = offsets('U2') + Ldiag_2;
offsets('U3') = offsets('L2') + Ldiag_2;
offsets('L3') = offsets('U3') + Ldiag_3;
offsets('U4') = offsets('L3') + Ldiag_3;
offsets('L4') = offsets('U4') + Ldiag_4;
offsets('U5') = offsets('L4') + Ldiag_4;
offsets('L5') = offsets('U5') + Ldiag_5;
offsets('U6') = offsets('L5') + Ldiag_5;
offsets('L6') = offsets('U6') + Ldiag_6;
offsets('U7') = offsets('L6') + Ldiag_6;
offsets('L7') = offsets('U7') + Ldiag_7;
offsets('TT') = offsets('L7') + Ldiag_7;

% explicit row coordinates
MTV(offsets('D0')+1:offsets('U1'), 1) = dof_sum + (1:Ldiag_0); % Main diagonal
MTV(offsets('U1')+1:offsets('L1'), 1) = dof_sum + (1:Ldiag_1); % First upper off-diagonal
MTV(offsets('L1')+1:offsets('U2'), 1) = dof_sum + (2:Ldiag_0); % First lower off-diagonal
MTV(offsets('U2')+1:offsets('L2'), 1) = dof_sum + (1:Ldiag_2); % Second upper off-diagonal
MTV(offsets('L2')+1:offsets('U3'), 1) = dof_sum + (3:Ldiag_0); % Second lower off-diagonal
MTV(offsets('U3')+1:offsets('L3'), 1) = dof_sum + (1:Ldiag_3); % Third upper off-diagonal
MTV(offsets('L3')+1:offsets('U4'), 1) = dof_sum + (4:Ldiag_0); % Third lower off-diagonal
MTV(offsets('U4')+1:offsets('L4'), 1) = dof_sum + (1:Ldiag_4); % Fourth upper off-diagonal
MTV(offsets('L4')+1:offsets('U5'), 1) = dof_sum + (5:Ldiag_0); % Fourth lower off-diagonal
tmp = 1:(Ldiag_4-1); tmp(4:4:end) = [];
MTV(offsets('U5')+1:offsets('L5'), 1) = dof_sum + (tmp); % Fifth upper off-diagonal
tmp = 6:Ldiag_0; tmp(4:4:end) = [];
MTV(offsets('L5')+1:offsets('U6'), 1) = dof_sum + (tmp); % Fifth lower off-diagonal
tmp = 1:(Ldiag_4-2); tmp([3:4:end 4:4:end]) = [];
MTV(offsets('U6')+1:offsets('L6'), 1) = dof_sum + (tmp); % Sixth upper off-diagonal
tmp = 7:Ldiag_0; tmp([3:4:end 4:4:end]) = [];
MTV(offsets('L6')+1:offsets('U7'), 1) = dof_sum + (tmp); % Sixth lower off-diagonal
MTV(offsets('U7')+1:offsets('L7'), 1) = dof_sum + (1:4:(Ldiag_4-3)); % Seventh upper off-diagonal
MTV(offsets('L7')+1:offsets('TT'), 1) = dof_sum + (8:4:Ldiag_0); % Seventh lower off-diagonal

% explicit column coordinates
MTV(offsets('D0')+1:offsets('U1'), 2) = dof_sum + (1:Ldiag_0); % Main diagonal
MTV(offsets('U1')+1:offsets('L1'), 2) = dof_sum + (2:Ldiag_0); % First upper off-diagonal
MTV(offsets('L1')+1:offsets('U2'), 2) = dof_sum + (1:Ldiag_1); % First lower off-diagonal
MTV(offsets('U2')+1:offsets('L2'), 2) = dof_sum + (3:Ldiag_0); % Second upper off-diagonal
MTV(offsets('L2')+1:offsets('U3'), 2) = dof_sum + (1:Ldiag_2); % Second lower off-diagonal
MTV(offsets('U3')+1:offsets('L3'), 2) = dof_sum + (4:Ldiag_0); % Third upper off-diagonal
MTV(offsets('L3')+1:offsets('U4'), 2) = dof_sum + (1:Ldiag_3); % Third lower off-diagonal
MTV(offsets('U4')+1:offsets('L4'), 2) = dof_sum + (5:Ldiag_0); % Fourth upper off-diagonal
MTV(offsets('L4')+1:offsets('U5'), 2) = dof_sum + (1:Ldiag_4); % Fourth lower off-diagonal
tmp = 6:Ldiag_0; tmp(4:4:end) = [];
MTV(offsets('U5')+1:offsets('L5'), 2) = dof_sum + (tmp); % Fifth upper off-diagonal
tmp = 1:(Ldiag_4-1); tmp(4:4:end) = [];
MTV(offsets('L5')+1:offsets('U6'), 2) = dof_sum + (tmp); % Fifth lower off-diagonal
tmp = 7:Ldiag_0; tmp([3:4:end 4:4:end]) = [];
MTV(offsets('U6')+1:offsets('L6'), 2) = dof_sum + (tmp); % Sixth upper off-diagonal
tmp = 1:(Ldiag_4-2); tmp([3:4:end 4:4:end]) = [];
MTV(offsets('L6')+1:offsets('U7'), 2) = dof_sum + (tmp); % Sixth lower off-diagonal
MTV(offsets('U7')+1:offsets('L7'), 2) = dof_sum + (8:4:Ldiag_0); % Seventh upper off-diagonal
MTV(offsets('L7')+1:offsets('TT'), 2) = dof_sum + (1:4:(Ldiag_4-3)); % Seventh lower off-diagonal

VCV(1:Ldiag_0) = dof_sum + (1:Ldiag_0);

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
