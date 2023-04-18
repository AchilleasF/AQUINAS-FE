%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script implements a simple scaling analysis of the parallel
% performance of AQUINAS on a simple LA problem with an increasing no. of
% processor threads (strong scaling - see Amdahl's Law) and increasing mesh
% size (weak scaling - see Gustafson's Law). The script reports the
% approximate serial fraction for each problem size and the figures are
% intended to aid in the choice of an 'optimal' no. of processor threads
% for parallel stiffness matrix assembly on the user's machine.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [1] - Amdahl, Gene M. "Validity of the single processor approach to achieving
% large scale computing capabilities." Proceedings of the April 18-20, 1967,
% spring joint computer conference. 1967.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.

%%%%%%%%%%%%%%%%%
% Preliminaries %
%%%%%%%%%%%%%%%%%
clearvars, clc, close all
addpath('..','../AQUINAS_Analysis','../AQUINAS_Objects','../Maple_Source');

unaNoEls = [100 500 1000 2000];
unMaxThreads = feature('numcores');
unNoRepeats = 20;


%%%%%%%%%%
% Inputs %
%%%%%%%%%%
% Geometry
r = 500.0; % [mm] - Cylinder radius
h = 1000.0; % [mm] - Cylinder height
t = 1.0; % [mm] - Cylinder thickness

% Material
MatName = 'Steel'; % Material name
E = 2.0e5; % [N/mm2] - Young's modulus of steel
nu = 0.3; % [-] - Poisson ratio of steel

% Plot controls
FS = 20; % [-] - Font Size (pts)
FSL = 18; % [-] - Font Size for Legends (pts)
MS = 5; % [-] - Markser size (pts)
LW = 3; % [-] - Line Width (pts)


%%%%%%%%%%%%%%%%%%%%
% AQUINAS solution %
%%%%%%%%%%%%%%%%%%%%
% The user is free to input any AQUINAS commands into the code below,
% including investigating the parallel scaling of LBA and G(M)N(I)A assembly.
M1 = AQUINAS_Material_Object('name',MatName,'E',E,'nu',nu,'G',E/(2*(1+nu)));
C1 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'u'});
C2 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',0.0,'dofs',{'w'});
C3 = AQUINAS_Constraint_Object('rcoord',r,'zcoord',h,'dofs',{'u'});
A1 = AQUINAS_Analysis_Object('type','LA');

daTimes_AV = zeros(length(unaNoEls),unMaxThreads);
daTimes_SE = zeros(length(unaNoEls),unMaxThreads);

for I = 1:length(unaNoEls)
    S1 = AQUINAS_Segment_Object('type','Cone','material',MatName,'thickness',t,'rzbot',[r 0.0],'rztop',[r h],'els',unaNoEls(I));

    for T = 1:unMaxThreads
        SOL1 = AQUINAS_Solver_Object('compiler','C++','noThreads',T);

        times = zeros(unNoRepeats,1);
        for J = 1:unNoRepeats
            disp([num2str(unaNoEls(I)),' elements with ',num2str(T),' threads run ',num2str(J),' of ',num2str(unNoRepeats)]);

            tstart = tic;
            AQUINAS_Protocol(M1,S1,C1,C2,C3,SOL1,A1);
            times(J) = toc(tstart);
        end
        daTimes_AV(I,T) = mean(times); % Average time
        daTimes_SE(I,T) = std(times)/sqrt(unNoRepeats); % Standard error
    end
end

%%%%%%%%%%%%%%%%%%%%
% Scaling analysis %
%%%%%%%%%%%%%%%%%%%%
daScaling_AV = daTimes_AV(:,1)./daTimes_AV; % Average scaling
daScaling_SE = daScaling_AV.*sqrt((daTimes_SE./daTimes_AV).^2 + (daTimes_SE(:,1)./daTimes_AV(:,1)).^2); % Standard error

% Simple Gauss-Newton fit of Amdahl's Law [1] to each of the scaling relationships
% See https://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm for procedure
daSerialFractions = 0.01*ones(length(unaNoEls),1); % Initial values of serial fraction - 1%
dTol = 1e-9; % Accepted tolerance
for I = 1:length(unaNoEls)
    while true
        dSerFracOld = daSerialFractions(I);
        J = (((1:unMaxThreads).*((1:unMaxThreads)-1))./(1+((1:unMaxThreads)-1)*dSerFracOld).^2)';
        dSerFracNew = dSerFracOld - (J'*J)\J'*(daScaling_AV(I,:) - ((1:unMaxThreads)./(dSerFracOld*(1:unMaxThreads)+(1-dSerFracOld))))';
        daSerialFractions(I) = dSerFracNew;
        if abs(dSerFracNew - dSerFracOld) < dTol
            break
        end
    end
end


%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%
%% Parallel scaling analysis
figure('Name','Parallel scaling analysis plot','NumberTitle','off','WindowState','Maximized','color','w');
cols = {'k', 'r', 'b', 'g'}; markers = {'o', 's', 'd', '^'};
lgnd = cell(length(unaNoEls),1);

subplot(1,2,1);
for I = 1:length(unaNoEls)
    errorbar(1:unMaxThreads, daTimes_AV(I,:), daTimes_SE(I,:),'LineStyle','-','LineWidth',LW,'Color',cols{I},'Marker',markers{I},'MarkerSize',MS); grid on; hold on;
    lgnd{I} = [num2str(unaNoEls(I)),' elements'];
end
legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('No. of threads','interpreter','latex','fontsize',FS);
ylabel('Runtime [s]','interpreter','latex','fontsize',FS);
set(gca,'YScale','log','ticklabelinterpreter','latex','fontsize',FS,'tickdir','out');


subplot(1,2,2);
lgnd = cell(length(unaNoEls),1);
for I = 1:length(unaNoEls)
    plot(1:unMaxThreads, (1:unMaxThreads)./(daSerialFractions(I)*(1:unMaxThreads)+(1-daSerialFractions(I))),'LineStyle','--','LineWidth',LW,'Color',cols{I}); grid on; hold on;
    lgnd{I} = ['Serial fraction = ',num2str(round(daSerialFractions(I)*100,2)),' $\%$'];
end
for I = 1:length(unaNoEls)
    errorbar(1:unMaxThreads, daScaling_AV(I,:), daScaling_SE(I,:),'LineStyle','-','LineWidth',LW,'Color',cols{I},'Marker',markers{I},'MarkerSize',MS);
end
legend(lgnd,'interpreter','latex','fontsize',FSL,'location','best');
xlabel('No. of threads','interpreter','latex','fontsize',FS);
ylabel('Parallel scaling and Amdahl''s Law fits','interpreter','latex','fontsize',FS);
set(gca,'ticklabelinterpreter','latex','fontsize',FS,'tickdir','out');

% Only the very expensive matrix assembly is currently parallelised in
% AQUINAS. This is returned in sparse format and subsequently used with the
% Matlab sparse linear algebra capabilities. For axisymmetric FE, the
% solution cost of the sparse linear system is typically negligible
% compared to the assembly cost, so this portion is currently not
% parallelised.

% A typical plot will show an improvement of overall performance with
% increasing no. of process threads assigned to the matrix assembly (strong
% scaling), but with an upper bound put on this scaling that is the inverse
% of the serial fraction as per Amdahl's Law. The serial fraction is
% typically 10 to 20% of the runtime and is not currently parallelised,
% meaning that parallelisation will give you at best a 5 to 10-fold speedup
% even with an infinite no. of processors. Increasing the mesh size will
% typically reduce the relative serial fraction which will improve your
% parallel scaling as per Gustafson's Law.


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