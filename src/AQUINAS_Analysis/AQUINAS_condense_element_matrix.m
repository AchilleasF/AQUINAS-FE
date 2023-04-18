function [MTV,VCV] = AQUINAS_condense_element_matrix(stiff_type,KT,KG,F,phi1,phi2,el,n,offsets,MTV,VCV)
% Function to condense a 12-size element system into an 8-size system
% AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.

% Taking a 12x12 matrix and condensing it into a 8x8 matrix
% The full 12x1 vector of nodal dofs is, in that order:
% [u1, (du/ds)1, v1, (dv/ds)1, w1, (dw/ds)1, u2, (du/ds)2, v2, (dv/ds)2, w2, (dw/ds)2]
% After transformation (Eq. 63; Rotter & Teng, Vol. 31), this becomes the 12x1 vector:
% [u1, beta1, v1, (dv/ds)1, w1, gamma1, u2, beta2, v2, (dv/ds)2, w2, gamma2]
% The condensed 8x1 vector of retained nodal dofs will be (Eq. 64; Rotter & Teng, Vol. 31), in that order:
% [u1, v1, w1, beta1, u2, v2, w2, beta2] or entries [1,3,5,2,7,9,11,8] of the above
% The 4x1 vector of eliminated nodal dofs will be, in that order:
% [(dv/ds)1, gamma1, (dv/ds)2, gamma2] or entries [4,6,10,12] of the above

% KTbar = KRR - KRC*KCC\KRC';  % Eq. 67 8x8 condensed element tangent stiffness matrix - Rotter & Teng (Vol. 31)
% PHIbar = PHIR - KRC*KCC\PHIC; % Eq. 68 8x1 condensed element vector of residual forced - Rotter & Teng (Vol. 31)
% where
% KRR = KT([1,3,5,2,7,9,11,8],[1,3,5,2,7,9,11,8])
% KRC = KT([1,3,5,2,7,9,11,8],[4,6,10,12])
% KCR = KRC' = KT([4,6,10,12],[1,3,5,2,7,9,11,8])
% KCC = KT([4,6,10,12],[4,6,10,12])
% PHIR = PHI([1,3,5,2,7,9,11,8])
% PHIC = PHI([4,6,10,12])

c1 = cos(phi1); s1 = sin(phi1); c2 = cos(phi2); s2 = sin(phi2);

% KT coming in () is an already-integraded 12x12 matrix corresponding to the full 12x1 vector as above
% The transformed version of this is (where c = cos(phi) and s = sin(phi)):
KT = [KT(:,1)   -s1*KT(:,2)-c1*KT(:,6)     KT(:,3)   KT(:,4)   KT(:,5)   c1*KT(:,2)-s1*KT(:,6)     KT(:,7)   -s2*KT(:,8)-c2*KT(:,12)     KT(:,9)   KT(:,10)   KT(:,11)   c2*KT(:,8)-s2*KT(:,12)];
KT = [KT(1,:)'  (-s1*KT(2,:)-c1*KT(6,:))'  KT(3,:)'  KT(4,:)'  KT(5,:)'  (c1*KT(2,:)-s1*KT(6,:))'  KT(7,:)'  (-s2*KT(8,:)-c2*KT(12,:))'  KT(9,:)'  KT(10,:)'  KT(11,:)'  (c2*KT(8,:)-s2*KT(12,:))'];
% F = [F(1)  -s1*F(2)-c1*F(6)  F(3)  F(4)  F(5)  c1*F(2)-s1*F(6)  F(7)  -c2*F(12)-s2*F(8)  F(9)  F(10)  F(11)  -s2*F(12)+c2*F(8)]';
if ~isempty(KG)
    KG = [KG(:,1)   -s1*KG(:,2)-c1*KG(:,6)     KG(:,3)   KG(:,4)   KG(:,5)   c1*KG(:,2)-s1*KG(:,6)     KG(:,7)   -s2*KG(:,8)-c2*KG(:,12)     KG(:,9)   KG(:,10)   KG(:,11)   c2*KG(:,8)-s2*KG(:,12)];
    KG = [KG(1,:)'  (-s1*KG(2,:)-c1*KG(6,:))'  KG(3,:)'  KG(4,:)'  KG(5,:)'  (c1*KG(2,:)-s1*KG(6,:))'  KG(7,:)'  (-s2*KG(8,:)-c2*KG(12,:))'  KG(9,:)'  KG(10,:)'  KG(11,:)'  (c2*KG(8,:)-s2*KG(12,:))'];
end
if strcmp(stiff_type,'material') || strcmp(stiff_type,'tangent')
    if (n==0)
        K = KT([1,3,5,2,7,9,11,8],[1,3,5,2,7,9,11,8]) - KT([1,3,5,2,7,9,11,8],[6,12])*(KT([6,12],[6,12])\KT([6,12],[1,3,5,2,7,9,11,8]));
    else
        K = KT([1,3,5,2,7,9,11,8],[1,3,5,2,7,9,11,8]) - KT([1,3,5,2,7,9,11,8],[4,6,10,12])*(KT([4,6,10,12],[4,6,10,12])\KT([4,6,10,12],[1,3,5,2,7,9,11,8]));
    end
    % P = F([1,3,5,2,7,9,11,8]) - KT([1,3,5,2,7,9,11,8],[4,6,10,12])*(KT([4,6,10,12],[4,6,10,12])\F([4,6,10,12]));
elseif strcmp(stiff_type,'geometric')
    if (n==0)
        K = KG([1,3,5,2,7,9,11,8],[1,3,5,2,7,9,11,8]);
        K = K - KG([1,3,5,2,7,9,11,8],[6,12])*(KT([6,12],[6,12])\KT([6,12],[1,3,5,2,7,9,11,8]));
        K = K - (KG([1,3,5,2,7,9,11,8],[6,12])*(KT([6,12],[6,12])\KT([6,12],[1,3,5,2,7,9,11,8])))';
        K = K + KT([1,3,5,2,7,9,11,8],[6,12])*(KT([6,12],[6,12])\(KG([6,12],[6,12])*(KT([6,12],[6,12])\KT([6,12],[1,3,5,2,7,9,11,8]))));
    else
        K = KG([1,3,5,2,7,9,11,8],[1,3,5,2,7,9,11,8]);
        K = K - KG([1,3,5,2,7,9,11,8],[4,6,10,12])*(KT([4,6,10,12],[4,6,10,12])\KT([4,6,10,12],[1,3,5,2,7,9,11,8]));
        K = K - (KG([1,3,5,2,7,9,11,8],[4,6,10,12])*(KT([4,6,10,12],[4,6,10,12])\KT([4,6,10,12],[1,3,5,2,7,9,11,8])))';
        K = K + KT([1,3,5,2,7,9,11,8],[4,6,10,12])*(KT([4,6,10,12],[4,6,10,12])\(KG([4,6,10,12],[4,6,10,12])*(KT([4,6,10,12],[4,6,10,12])\KT([4,6,10,12],[1,3,5,2,7,9,11,8]))));
    end

end

% Map matrix & vector row & column positions to column entries
el_off = (el - 1) * 4; el_off1 = (el - 1) * 3; el_off2 = (el - 1) * 2; el_off3 = (el - 1);
if ~isempty(MTV)
    MTV((1:8) + el_off + offsets('D0'), 3) = MTV((1:8) + el_off + offsets('D0'), 3) + diag(K, 0); % Main diagonal
    MTV((1:7) + el_off + offsets('U1'), 3) = MTV((1:7) + el_off + offsets('U1'), 3) + diag(K, 1); % First upper off-diagonal
    MTV((1:7) + el_off + offsets('L1'), 3) = MTV((1:7) + el_off + offsets('L1'), 3) + diag(K,-1); % First lower off-diagonal
    MTV((1:6) + el_off + offsets('U2'), 3) = MTV((1:6) + el_off + offsets('U2'), 3) + diag(K, 2); % Second upper off-diagonal
    MTV((1:6) + el_off + offsets('L2'), 3) = MTV((1:6) + el_off + offsets('L2'), 3) + diag(K,-2); % Second lower off-diagonal
    MTV((1:5) + el_off + offsets('U3'), 3) = MTV((1:5) + el_off + offsets('U3'), 3) + diag(K, 3); % Third upper off-diagonal
    MTV((1:5) + el_off + offsets('L3'), 3) = MTV((1:5) + el_off + offsets('L3'), 3) + diag(K,-3); % Third lower off-diagonal
    MTV((1:4) + el_off + offsets('U4'), 3) = MTV((1:4) + el_off + offsets('U4'), 3) + diag(K, 4); % Fourth upper off-diagonal
    MTV((1:4) + el_off + offsets('L4'), 3) = MTV((1:4) + el_off + offsets('L4'), 3) + diag(K,-4); % Fourth lower off-diagonal
    MTV((1:3) + el_off1 + offsets('U5'), 3) = MTV((1:3) + el_off1 + offsets('U5'), 3) + diag(K, 5); % Fifth upper off-diagonal
    MTV((1:3) + el_off1 + offsets('L5'), 3) = MTV((1:3) + el_off1 + offsets('L5'), 3) + diag(K,-5); % Fifth lower off-diagonal
    MTV((1:2) + el_off2 + offsets('U6'), 3) = MTV((1:2) + el_off2 + offsets('U6'), 3) + diag(K, 6); % Sixth upper off-diagonal
    MTV((1:2) + el_off2 + offsets('L6'), 3) = MTV((1:2) + el_off2 + offsets('L6'), 3) + diag(K,-6); % Sixth lower off-diagonal
    MTV(1 + el_off3 + offsets('U7'), 3) = MTV(1 + el_off3 + offsets('U7'), 3) + K(1,8); % Seventh upper off-diagonal
    MTV(1 + el_off3 + offsets('L7'), 3) = MTV(1 + el_off3 + offsets('L7'), 3) + K(8,1); % Seventh lower off-diagonal
end

if ~isempty(F)
    F = [F(1)  -s1*F(2)-c1*F(6)  F(3)  F(4)  F(5)  c1*F(2)-s1*F(6)  F(7)  -c2*F(12)-s2*F(8)  F(9)  F(10)  F(11)  -s2*F(12)+c2*F(8)]';
    if (n==0)
        P = F([1,3,5,2,7,9,11,8]) - KT([1,3,5,2,7,9,11,8],[6,12])*(KT([6,12],[6,12])\F([6,12]));
    else
        P = F([1,3,5,2,7,9,11,8]) - KT([1,3,5,2,7,9,11,8],[4,6,10,12])*(KT([4,6,10,12],[4,6,10,12])\F([4,6,10,12]));
    end
    VCV((1:8) + el_off + offsets('D0'), 2) = VCV((1:8) + el_off + offsets('D0'), 2) + P; % Nodal load vector
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