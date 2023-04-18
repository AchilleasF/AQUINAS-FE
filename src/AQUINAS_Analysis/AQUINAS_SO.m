function [ EIGVECTORS, EIGVALUES, ANALYSIS ] = AQUINAS_SO(SEGMENTS,CONSTRAINTS,SOLVER,MATERIALS,ANALYSIS,DTs,elDOFs,Sigmas,OStorage)
% AQUINAS_SO - Surrogate Optimisation for the computation of the crtical bifurcation load and associated harmonic.
% AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.

    % Initialization of necessary fields for the Gaussian Process Regression algorithm
    noise = SOLVER.so_noise;
    Xpoints = (SOLVER.so_circWave_lb:SOLVER.so_circWave_ub)';
    Xtrain = unique(round(linspace(min(Xpoints),max(Xpoints),SOLVER.so_No_Initial_Points)))'; Ytrain = nan*Xtrain;
    stepHyperparameters = {};
    trainEigenValues = containers.Map('KeyType','uint32','ValueType','any');
    trainEigenVectors = containers.Map('KeyType','uint32','ValueType','any');
    dof_count = 0;
    for S = 1:length(SEGMENTS); dof_count = dof_count + (SEGMENTS(S).els + 1) * 4; end
    for C = 1:length(CONSTRAINTS); dof_count = dof_count + 1; end

    for I=1:length(Xtrain)
        ANALYSIS.circumferential_modes = Xtrain(I);
        [EIGVECTORS,EIGVALUES] = AQUINAS_LBA(SEGMENTS,CONSTRAINTS,SOLVER,MATERIALS,ANALYSIS,DTs,elDOFs,Sigmas,OStorage);
        trainEigenVectors(Xtrain(I)) = EIGVECTORS;
        trainEigenValues(Xtrain(I)) = EIGVALUES;
        Ytrain(I) = EIGVALUES(1);
    end

    XtrainNormalized = (Xtrain - mean(Xtrain))./std(Xtrain);
    XpointsNormalized = (Xpoints - mean(Xtrain))./std(Xtrain); XpointsNormalizedRange = max(XpointsNormalized) - min(XpointsNormalized);
    YtrainNormalized = (Ytrain - mean(Ytrain))./std(Ytrain);

    if (SOLVER.so_Hyperparameter_Optimisation)
        hyperparameters_lb = [SOLVER.so_Hyperparameters_Optimisation_lbf(1)*XpointsNormalizedRange, SOLVER.so_Hyperparameters_Optimisation_lbf(2)*std(YtrainNormalized)];
        hyperparameters_ub = [SOLVER.so_Hyperparameters_Optimisation_ubf(1)*XpointsNormalizedRange, SOLVER.so_Hyperparameters_Optimisation_ubf(2)*std(YtrainNormalized)];
        [hyperparameters,~,~,~] = evalThetas(XtrainNormalized,YtrainNormalized,hyperparameters_lb,hyperparameters_ub,noise,SOLVER.so_Hyperparameters_Optimisation_Grid_Stations);
    else
        hyperparameters = SOLVER.so_Hyperparameters;
    end
    [ypred,ysd] = gpr(XtrainNormalized,YtrainNormalized,XpointsNormalized,hyperparameters,noise);
    stepHyperparameters{end+1} = hyperparameters;

    for I = 1:SOLVER.so_Max_Steps

        % Learning Function
        [new_sample, EI, maxEI] = EI_learning(Xpoints, Xtrain, ypred, ysd);
        if (maxEI<SOLVER.so_EI_tol)
            minInd = find(Ytrain==min(Ytrain)); if length(minInd)>1; break; end
            if ~any(ismember(Xtrain,Xtrain(minInd)-1)) && any(ismember(Xpoints,Xtrain(minInd)-1))
                new_sample = Xtrain(minInd) - 1;
            elseif ~any(ismember(Xtrain,Xtrain(minInd)+1)) && any(ismember(Xpoints,Xtrain(minInd)+1))
                new_sample = Xtrain(minInd) + 1;
            else
                break;
            end
        end

        % Update and sorting training data
        Xtrain(end+1) = new_sample;
        ANALYSIS.circumferential_modes = new_sample;
        [EIGVECTORS,EIGVALUES] = AQUINAS_LBA(SEGMENTS,CONSTRAINTS,SOLVER,MATERIALS,ANALYSIS,DTs,elDOFs,Sigmas,OStorage);
        trainEigenVectors(new_sample) = EIGVECTORS;
        trainEigenValues(new_sample) = EIGVALUES;
        Ytrain(end+1) = EIGVALUES(1);
        [~,Xsort] = sort(Xtrain);
        Xtrain = Xtrain(Xsort);  Ytrain = Ytrain(Xsort);
        % Update training space
        if SOLVER.so_Auto_Bounds
            Xpoints = Xpoints(1) : max(Xpoints(end),Xtrain(end) + SOLVER.so_circWave_bound_step);
        end
        % Normalisation
        XtrainNormalized = (Xtrain - mean(Xtrain))./std(Xtrain);
        XpointsNormalized = (Xpoints - mean(Xtrain))./std(Xtrain); XpointsNormalizedRange = max(XpointsNormalized) - min(XpointsNormalized);
        YtrainNormalized = (Ytrain - mean(Ytrain))./std(Ytrain);
        % Hyperparameter Optimisation
        if (SOLVER.so_Hyperparameter_Optimisation)
            hyperparameters_lb = [SOLVER.so_Hyperparameters_Optimisation_lbf(1)*XpointsNormalizedRange, SOLVER.so_Hyperparameters_Optimisation_lbf(2)*std(YtrainNormalized)];
            hyperparameters_ub = [SOLVER.so_Hyperparameters_Optimisation_ubf(1)*XpointsNormalizedRange, SOLVER.so_Hyperparameters_Optimisation_ubf(2)*std(YtrainNormalized)];
            [hyperparameters,mu,sigma,lval] = evalThetas(XtrainNormalized,YtrainNormalized,hyperparameters_lb,hyperparameters_ub,noise,SOLVER.so_Hyperparameters_Optimisation_Grid_Stations);
        end
        % Gaussian Process Regression
        [ypred,ysd] = gpr(XtrainNormalized,YtrainNormalized,XpointsNormalized,hyperparameters,noise);
        stepHyperparameters{end+1} = hyperparameters;

    end
    % Separate consideration of the axisymmetric buckling mode (0th harmonic)
    ANALYSIS.circumferential_modes = 0;
    [EIGVECTORS,EIGVALUES] = AQUINAS_LBA(SEGMENTS,CONSTRAINTS,SOLVER,MATERIALS,ANALYSIS,DTs,elDOFs,Sigmas,OStorage);
    trainEigenVectors(0) = EIGVECTORS;
    trainEigenValues(0) = EIGVALUES;

    ANALYSIS.circumferential_modes = [ 0 Xtrain' ];
    EIGVALUES = nan*zeros(ANALYSIS.no_eigenvalues,trainEigenValues.Count);
    EIGVECTORS = nan*zeros(dof_count,ANALYSIS.no_eigenvalues,trainEigenVectors.Count);
    EIGVALUES(:,1) = trainEigenValues(0);
    EIGVECTORS(:,:,1) = trainEigenVectors(0);
    for I = 2:trainEigenValues.Count
        EIGVALUES(:,I) = trainEigenValues(Xtrain(I-1));
        EIGVECTORS(:,:,I) = trainEigenVectors(Xtrain(I-1));
    end

end

function [new_sample, EI, maxEI] = EI_learning(Xpoints, Xtrain, ypred, ysd)

    ysd = ysd + 1e-10; % to prevent zero standard deviation value

    % 1-Find the current minimum
    min_current_objective = min(ypred);

    % 2-Calculate the EI values of the candidate samples
    EI = (min_current_objective-ypred).*normcdf((min_current_objective-ypred)./ysd) + ysd.*normpdf((min_current_objective-ypred)./ysd);

    % 3-Select a new sample
    [~,candidates] = setdiff(Xpoints,Xtrain);
    maxEI = max(EI(candidates));
    new_sample = Xpoints(find(EI==maxEI,1));

end

function [mu,s] = gpr(Xtrain,Ytrain,X,hyperparameters,noise)

    l = hyperparameters(1);
    sigma = hyperparameters(2);

    K = kernel(Xtrain,Xtrain,l,sigma);
    L = chol(K+sigma^2*noise*eye(length(Xtrain)))';
    Lk = L\kernel(Xtrain,X,l,sigma);

    % predicted
    mu = Lk'*(L\Ytrain);
    K_ = kernel(X,X,l,sigma);
    s2 = diag(K_) - sum(Lk.^2)';
    s = sqrt(s2);

end

function [K] =  kernel(X1, X2, l, sigma)

    K = zeros(length(X1),length(X2));
    for I = 1:length(X1)
        for J = 1:length(X2)
            K(I,J) = (sigma^2)*exp(-((X1(I)-X2(J))^2)/2/l/l); % Radial Basis Function kernel
        end
    end

end

function [thetas,mu,sigma,lval] = evalThetas(X,Y,hyperparameters_lb,hyperparameters_ub,noise,grid_stations)

    [grid_l,grid_sigma] = meshgrid(logspace(log10(hyperparameters_lb(1)),log10(hyperparameters_ub(1)),grid_stations(1)),logspace(log10(hyperparameters_lb(2)),log10(hyperparameters_ub(2)),grid_stations(2)));
    f = @(thetas) -logMarginalLikelihood(Y,X,X,thetas,noise);
    gridLML = nan*grid_l;
    for I = 1:size(grid_l,1)
        for J = 1:size(grid_l,2)
            gridLML(I,J) = f([grid_l(I,J) grid_sigma(I,J) noise]);
        end
    end
    gridLML = real(gridLML);
    [I,J] = find(gridLML == min(min(gridLML)),1);
    hyperparameters = [grid_l(I,J) grid_sigma(I,J)];

    [lval] = logMarginalLikelihood(Y,X,X,hyperparameters,noise);

    K = kernel(X,X,hyperparameters(1),hyperparameters(2));
    [n,~] = size(Y);

    s = K\ones(n,1);

    r = K\Y;

    mu = sum(r)/sum(s);
    var = (1/(n))*(Y - mu*ones(n,1))'*r;
    sigma = sqrt(var);

    thetas = [hyperparameters(1),hyperparameters(2)];

end

function [lval] = logMarginalLikelihood(y,X1,X2,hyperparameters,noise)

    l = hyperparameters(1);
    sigma = hyperparameters(2);

    [n,~] = size(y);
    K = kernel(X1, X2, l, sigma);
    Ky = K + noise*eye(n);

    lval = -1/2*(y'*(Ky\y) + log(det(Ky)) + n*log(2*pi));

end

% The following functions are coded to avoid dependency to MATLAB's 'Statistics and Machine Learning Toolbox'
function pdf = normpdf(x)
    % Normal Gaussian probability density function
    pdf = 1/sqrt(2*pi).*exp(-x.*x/2);
end

function cdf = normcdf(x)
    % Normal Gaussian cumulative distribution function
    cdf = 0.5*(1+erf(x/sqrt(2)));
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