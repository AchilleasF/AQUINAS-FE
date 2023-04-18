function DELTA = AQUINAS_vaporise_element_matrix(KT,F,DELTAR,phi1,phi2,n)
% Function to vaporise (undoing of the static condensation) an 8-size vector of condensed
% element global nodal dofs into the original 12-size vector based on the retained stiffness matrix
% AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.

% The incoming 8x1 vector DELTAR has the following structure
% [u1, v1, w1, beta1, u2, v2, w2, beta2]
% The incoming KT is an already-integraded 12x12 matrix corresponding to the full 12x1 vector as above
% The transformed version of this is (where c = cos(phi) and s = sin(phi)):

c1 = cos(phi1); s1 = sin(phi1); c2 = cos(phi2); s2 = sin(phi2);

KT = [KT(:,1)   -s1*KT(:,2)-c1*KT(:,6)     KT(:,3)   KT(:,4)   KT(:,5)   c1*KT(:,2)-s1*KT(:,6)     KT(:,7)   -s2*KT(:,8)-c2*KT(:,12)     KT(:,9)   KT(:,10)   KT(:,11)   c2*KT(:,8)-s2*KT(:,12)];
KT = [KT(1,:)'  (-s1*KT(2,:)-c1*KT(6,:))'  KT(3,:)'  KT(4,:)'  KT(5,:)'  (c1*KT(2,:)-s1*KT(6,:))'  KT(7,:)'  (-s2*KT(8,:)-c2*KT(12,:))'  KT(9,:)'  KT(10,:)'  KT(11,:)'  (c2*KT(8,:)-s2*KT(12,:))'];
F = [F(1)  -s1*F(2)-c1*F(6)  F(3)  F(4)  F(5)  c1*F(2)-s1*F(6)  F(7)  -s2*F(8)-c2*F(12)  F(9)  F(10)  F(11)  c2*F(8)-s2*F(12)]';

% Recovery of original order of transformed vector, which has the structure:
% [u1, beta1, v1, (dv/ds)1, w1, gamma1, u2, beta2, v2, (dv/ds)2, w2, gamma2]
DELTA = zeros(12,1); DELTA([1,3,5,2,7,9,11,8]) = DELTAR;

% Recovery of the 4x1 vector of condensed dofs, which has the structure:
% [(dv/ds)1, gamma1, (dv/ds)2, gamma2]
if (n==0)
    % For the axisymmetric case (n=0) the values of (dv/ds)1 and (dv/ds)2 are zero by default
    DELTA([6,12]) = KT([6,12],[6,12])\(F([6,12]) - KT([6,12],[1,3,5,2,7,9,11,8])*DELTAR);
else
    DELTA([4,6,10,12]) = KT([4,6,10,12],[4,6,10,12])\(F([4,6,10,12]) - KT([4,6,10,12],[1,3,5,2,7,9,11,8])*DELTAR);
end


% Recovery of untransformed vector, which as the structure:
% [u1, (du/ds)1, v1, (dv/ds)1, w1, (dw/ds)1, u2, (du/ds)2, v2, (dv/ds)2, w2, (dw/ds)2]
DELTA = [DELTA(1)  -s1*DELTA(2)+c1*DELTA(6)  DELTA(3)  DELTA(4)  DELTA(5)  -c1*DELTA(2)-s1*DELTA(6)  DELTA(7)  -s2*DELTA(8)+c2*DELTA(12)  DELTA(9)  DELTA(10)  DELTA(11)  -c2*DELTA(8)-s2*DELTA(12)]';

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