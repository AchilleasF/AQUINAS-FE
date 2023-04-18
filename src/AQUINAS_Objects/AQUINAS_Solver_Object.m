classdef AQUINAS_Solver_Object
    % Class definition for a Solver Object - used for defining the parameters of the
    % solution procedure.
    % Using the AQUINAS_Solver_Object the User can decide on the number
    % of integration points to be used according to the Gauss Legendre
    % integration scheme, or the compiler to be used in the numerical
    % integration of the stiffnesses.
    % In addition, the AQUINAS_Solver_Object is responsible for the
    % definition of the parameters needed for finding the critical
    % circumferential wavenumber that the axisymmetric shell will buckle into,
    % using a Gaussian Process Regression algorithm.
    % AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.

    properties

        % Miscellaneous
        class_type = 'AQUINAS_Solver_Object';
        compiler % Option for choosing between Matlab and C++ as the compiling language for more computationally demanding operations
        console_output % A boolean value that is true if the User wants to view some basic output displayed in the console window and false otherwise
        No_Threads % No. of OpenMP threads (C++ interface only)

        % Numerical integration
        No_Gauss_Stations % Number of integration stations for Gauss-Legendre integration
        No_Simpson_Stations % Number of through-thickness stations for application of the Simpson integration rule
        Gauss_Nodes % Gauss nodes
        Gauss_Weights % Associated Gauss weights
        Lobatto_Nodes % Lobatto nodes
        Lobatto_Weights % Associated Lobatto weights

        % Bifurcation analysis
        eigsSigma % Value to be used as the sigma parameter during the call to built-in function eigs for computing the eigenvalues closest to sigma
        eigsTolerance % Value to be used in defining the 'Tolerance' property of the Matlab built in eigs function
        eigsSubspaceDimensionMultiplier % Factor used alongside the ANALYSIS.no_eigenvalues property to define the 'SubspaceDimension' of the Matlab built in eigs function
        eigsNegativeTreatment % Option that dictates whether negative eigenvalues should be kept or dropped during bifurcation analysis (takes values of 'Keep' or 'Drop' accordingly)
        eigsForceReal % Boolean option that allows the LBA solver to only keep the real part of the eigenvalues and eigenvectors, in the case of complex solutions returned by eigs
        Constr_Langr_Multip_Factor % The Factor that the terms corresponding to the enforcement of Constraints, according to the Lagrange Multipliers method, will be multiplied with.
                                   % This multiplication will lead to eigenvalues corresponding to the Constraints equal to -1/Constr_Langr_Multip_Factor,
                                   % giving the analyst the control to lead those eigenvalues as far as possible from the physically important ones.

        % Nonlinear analysis
        flag % Integer flag that is zero if no problem occurs during the analysis, positive for a warning and negative if an error is encountered. Counter-measures will be undertaken from different parts of the solver depending on the value of the flag.
        nonlinear_solver % Option to decide whether to use a Newton-Raphson or Arc-Length (Riks) method for the nonlinear computations of the equilibrium path
        arclength_method % Option that determines whether a Crisfield-Riks [Crisfield (1981)] or Ramm's [Ramm (1981)] method is to be employed for determining the load increment according to the arc length method. Can be either 'Ramm' or 'Crisfield'
        no_max_steps % Maximum number of steps allowed for incrementing the load in the nonlinear analysis
        Jmax % Maximum number of iterations allowed to achieve convergence of each applied load step
        Jd % Desired number of iterations to achieve convergence of each applied load step
        LPF % The Load Proportionality Factor that the nonlinear solver has reached at any point in the analysis
        dLPF % The increment of the externally applied load that is applied per step of the nonlinear analysis. This can be updated while the analysis proceeds, based on the type (Newton-Raphson or Arc-Length) of the solver and other options
        LPFmax % The maximum Load Proportionality Factor that the nonlinear solver is allowed to reach. Execution will stop once this LPF has been excedeed.
        dLPFmax % The maximum Load Proportionality Factor increment that the nonlinear solver is allowed to attempt. Higher dLPF increment (in absolute value) will be reduced to this value. This works in parallel to the lmax variable, relevant to the arc length solver.
        l % The increment of the arc length method that the nonlinear solver will attempt to solve for during each load step. Only relevant for the arc length solver.
        lmax % The maximum arc length increment that the nonlinear solver is allowed to attempt. Only relevant for the arc length solver. This is evaluated in the first iteration of the first increment using the 'dLPFmax' option provided through the constructor.
        zeta_t % Tolerance to test for convergence of an applied load increment
        zeta_d % Tolerance to test for divergence of an applied load increment
        ksi % Tolerance parameter for comparison of eigenvalues with unity, according to Gupta's method.
        epsilon_s % Reference strain to determine number of sub increments for materially nonlinear analyses
        max_epsilon_bar % Maximum effective strain increment that may be attempted to be solved for with the sub-increment method (corresponding to a number of sub-increments Nsb = max_epsilon_bar/epsilon_s sub-increments)
        check_axisym_stability % Boolean that determines whether the eigenvalue of the tangent stiffness matrix after solution of a load increment is to be examined. A negative eigenvalue means a loss of stability at the axisymmetric prebuckling path has occurred
        axisym_path_eig % The eigenvalue of the tangent stiffness matrix at the axisymmetric prebuckling path. Only relevant if check_axisym_stability is set to true
        simultaneous_bifurcation_treatment % Option that determines how to deal with simultaneous bifurcation modes occuring for eigenvalues in the range 1-ksi < lambda < 1+ksi. Possible values are 'none', 'pickSmallestWithinRange' and 'divideKsiBy10'
        check_bifurcation_after_LPF % The LPF after which the shell will be submitted for a bifurcation check for every applied load step
        max_attempts % Number of consecutive failed (not converged) steps to be reattempted, before giving up and stopping execution
        reattempt_cutback % Factor to multiply the dLPF (for Newton-Raphson) or arc-length (for arc-length) with, when a load step is deemed as failed and going to be re-attempted. Must be in the range (0,1)
        cip % struct containing all the necessary parameters and variables required to check for the termination of a nonlinear analysis based on the development of the Convergence Indicator Plot (CIP)
        visualise_nonlinear % Boolean that determines whether a figure window will be generated to supervise the progress of the nonlinear solvers
        termination_conditions % Character 1D array that holds letters corresponding to the criteria that will determine the termination of the nonlinear analysis. Potential conditions, with their corresponding code letters, are:
                               % A  -  Maximum number of load steps reached
                               % B  -  Maximum LPF exceeded
                               % C  -  Bifurcation check, terminate at first event. The circumferential modes to be considered should either be provided through the AQUINAS_Analysis_Object or found by turning on the minimisation algorithm (Gaussian Process Regression, see following properties)
                               % D  -  Limit point, where the global stiffness matrix is not positive definite any more
                               % E  -  Prediction of plateauing region, based on the development of the Convergence Indicator Plot (CIP). For methodology and references see the convergence_indicator_plot_termination_check method below.
                               % N.B that some of the termination condition require the additional definition of another property, for example for applying the 'A' condition the maximum number of load steps allowed needs also to be provided through the 'noMaxSteps' name-value pair.
                               % Any combination of termination conditions can be provided in a character array, for example 'AC' will test for the maximum number of load steps reached and the possibility of bifurcation per converged load step

        % Minimisation related properties, using a Surrogate optimisation method (Gaussian Process Regression)
        surrogate_optimisation % A boolean value that is true if the GPR algorithm needs to be employed for finding the critical circumferential wavenumber
        so_Auto_Bounds % A boolean value that is true if the User does not provide bounds for searching the critical circumferential wavenumber
        so_circWave_lb % Lower bound of circumferential wavenumbers to be considered (cannot be omitted if so_Auto_Bounds is false). N.B. the axisymmetric mode (zero circumferential wavenumber) is treated separately and not included in the lower bound.
        so_circWave_ub % Upper bound of circumferential wavenumbers to be considered (cannot be omitted if so_Auto_Bounds is false)
        so_circWave_bound_step % The step to be used for extending the upper bound of circumferential waves to be considered, in the case of auto bound identification
        so_Hyperparameter_Optimisation % A boolean value that is true if the User wishes for an optimisation of the hyperparameters to be used in the GPR algorithm by maximizing the
                                        % log marginal likelihood function. If this value is false then the provided hyperparameters will be treated as constant and used throughout the algorithm
        so_Hyperparameters % The Hyperparameters to be used in the GPR algorithm. If the so_Hyperparameter_Optimisation property is true, then the Hyperparameter values will be changed as the minimisation process progresses.
                            % The Hyperparameters, in the order they need to be provided as input, are: [ l (length scale of RBF kernel), sigma (sd of RBF kernel), noise (the noise of the train data, should be really small, for example 1e-10) ]
        so_Hyperparameters_Optimisation_Grid_Stations % Stations to be used for the grid search during the hyperparameter optimisation
        so_Hyperparameters_Optimisation_lbf % Lower bound factors for hyperparameters. These factors are multiplied with the range of normalized Xpoints and with the standard deviation
                                             % of the normalized Y trained data to set the lower bounds for the length scale and sigma hyperparameters accordingly. Only relevant if so_Hyperparameter_Optimisation is true
        so_Hyperparameters_Optimisation_ubf % Upper bound factors for hyperparameters. These factors are multiplied with the range of normalized Xpoints and with the standard deviation
                                             % of the normalized Y trained data to set the upper bounds for the length scale and sigma hyperparameters accordingly. Only relevant if so_Hyperparameter_Optimisation is true
        so_noise % The assumed noise of the Y data. This should be a rather small value since the eigenvalues are obtained through a thoroughly deterministic process. The noise only exists in this case to make sure the covariance matrix is positive definite.
                  % It is important to note that the 'noise' parameter is usually treated as a hyperparameter in other implementations of the GPR algorithm.
        so_No_Initial_Points % Number of points to be used for the initial training set of data for the GPR algorithm
        so_alpha % Significance level
        so_EI_tol % Tolerance for termination of the minimisation process based on the Learning function output
        so_Max_Steps % No of max steps to be allowed before the surrogate optimisation algorithm is terminated

    end

    methods
        % Constructor method for a Solver Object
        function obj = AQUINAS_Solver_Object(options)
            arguments
                % Miscellaneous inputs
                options.compiler { mustBeNonzeroLengthText } = 'C++'
                options.consoleOutput { mustBeNumericOrLogical } = false
                options.noThreads { mustBeInteger,mustBePositive } = 1

                % Numerical integration related inputs
                options.noGaussStations { mustBeInteger,mustBePositive } = 3
                options.noSimpsonStations { mustBeInteger,mustBePositive } = 7

                % Bifurcation analysis related inputs
                options.eigsSigma { mustBeReal} = 0
                options.eigsTolerance { mustBeReal,mustBePositive } = 1e-6
                options.eigsSubspaceDimensionMultiplier { mustBeInteger,mustBePositive } = 2
                options.eigsNegativeTreatment { mustBeNonzeroLengthText } = 'Keep'
                options.eigsForceReal { mustBeNumericOrLogical } = true
                options.constrLagrMultF { mustBeReal } = 1

                % Nonlinear analysis related inputs
                options.NonlinearSolver { mustBeNonzeroLengthText } = 'ArcLength'
                options.ArclengthMethod  { mustBeNonzeroLengthText } = 'Ramm'
                options.noMaxSteps { mustBeInteger,mustBePositive } = []
                options.Jd { mustBeInteger,mustBePositive } = 3
                options.Jmax { mustBeInteger,mustBePositive } = 20
                options.dLPF { mustBeReal,mustBePositive } = []
                options.LPFmax { mustBeReal,mustBePositive } = Inf
                options.dLPFmax { mustBeReal,mustBePositive } = Inf
                options.zetaT { mustBeReal,mustBePositive } = 1e-3
                options.zetaD { mustBeReal,mustBePositive } = 2
                options.epsilon_s { mustBeReal,mustBePositive } = 2e-4
                options.max_epsilon_bar { mustBeReal,mustBePositive } = 2e-2
                options.ksi { mustBeReal,mustBePositive } = 1e-3
                options.checkAxisymStability = false;
                options.simultaneousBifurcationTreatment { mustBeNonzeroLengthText } = 'none'
                options.bifurcationAfterLPF { mustBeReal,mustBeNonnegative } = 0
                options.maxAttempts { mustBeInteger,mustBeNonnegative } = 0
                options.reattemptCutback { mustBeReal,mustBePositive } = 0.25
                options.DOFTrackerForCIP {mustBeNonzeroLengthText} = '-'
                options.visualiseNonlinear { mustBeNumericOrLogical } = false
                options.terminationConditions { mustBeNonzeroLengthText } = '-'

                % Surrogate optimisation related inputs
                options.surrogate_optimisation { mustBeNumericOrLogical } = false
                options.so_Auto_Bounds { mustBeNumericOrLogical } = true
                options.so_circWave_lb { mustBeInteger,mustBePositive } = 1
                options.so_circWave_ub { mustBeInteger,mustBePositive } = 100
                options.so_circWave_bound_step { mustBeInteger,mustBePositive } = 10
                options.so_Hyperparameter_Optimisation { mustBeNumericOrLogical } = true
                options.so_Hyperparameters = []
                options.so_Hyperparameters_Optimisation_Grid_Stations = [50 50]
                options.so_Hyperparameters_Optimisation_lbf = [0.05 0.5]
                options.so_Hyperparameters_Optimisation_ubf = [0.25 2.5]
                options.so_noise { mustBeReal,mustBePositive } = 1e-10
                options.so_No_Initial_Points { mustBeInteger,mustBePositive } = 10
                options.so_alpha { mustBeReal,mustBePositive } = 0.05
                options.so_EI_tol { mustBeReal,mustBePositive } = 1e-5
                options.so_Max_Steps { mustBeInteger,mustBePositive } = 1e5
            end
            % Miscellaneous
            obj.compiler = options.compiler;
            obj.console_output = options.consoleOutput;
            obj.No_Threads = options.noThreads;
            % Numerical integration
            obj.No_Simpson_Stations = options.noSimpsonStations;
            obj.No_Gauss_Stations = options.noGaussStations;
            % Bifurcation analysis
            obj.eigsSigma = options.eigsSigma;
            obj.eigsTolerance = options.eigsTolerance;
            obj.eigsSubspaceDimensionMultiplier = options.eigsSubspaceDimensionMultiplier;
            obj.eigsNegativeTreatment = options.eigsNegativeTreatment;
            obj.eigsForceReal = options.eigsForceReal;
            obj.Constr_Langr_Multip_Factor = options.constrLagrMultF;
            % Nonlinear analysis
            obj.flag = 0;
            obj.nonlinear_solver = options.NonlinearSolver;
            obj.arclength_method = options.ArclengthMethod;
            obj.no_max_steps = options.noMaxSteps;
            obj.Jmax = options.Jmax;
            obj.Jd = options.Jd;
            obj.dLPF = options.dLPF;
            obj.LPFmax = options.LPFmax;
            obj.dLPFmax = options.dLPFmax;
            obj.zeta_t = options.zetaT;
            obj.zeta_d = options.zetaD;
            obj.ksi = options.ksi;
            obj.epsilon_s = options.epsilon_s;
            obj.max_epsilon_bar = options.max_epsilon_bar;
            obj.check_axisym_stability = options.checkAxisymStability;
            obj.simultaneous_bifurcation_treatment = options.simultaneousBifurcationTreatment;
            obj.check_bifurcation_after_LPF = options.bifurcationAfterLPF;
            obj.max_attempts = options.maxAttempts;
            obj.reattempt_cutback = options.reattemptCutback;
            obj.cip.dof_tracker = options.DOFTrackerForCIP;
            obj.visualise_nonlinear = options.visualiseNonlinear;
            obj.termination_conditions = options.terminationConditions;
            % Surrogate optimisation with Gaussian process regression
            obj.surrogate_optimisation = options.surrogate_optimisation;
            obj.so_Auto_Bounds = options.so_Auto_Bounds;
            obj.so_circWave_lb = options.so_circWave_lb;
            obj.so_circWave_ub = options.so_circWave_ub;
            obj.so_circWave_bound_step = options.so_circWave_bound_step;
            obj.so_Hyperparameter_Optimisation = options.so_Hyperparameter_Optimisation;
            obj.so_Hyperparameters = options.so_Hyperparameters;
            obj.so_Hyperparameters_Optimisation_Grid_Stations = options.so_Hyperparameters_Optimisation_Grid_Stations;
            obj.so_Hyperparameters_Optimisation_lbf = options.so_Hyperparameters_Optimisation_lbf;
            obj.so_Hyperparameters_Optimisation_ubf = options.so_Hyperparameters_Optimisation_ubf;
            obj.so_noise = options.so_noise;
            obj.so_No_Initial_Points = options.so_No_Initial_Points;
            obj.so_alpha = options.so_alpha;
            obj.so_EI_tol = options.so_EI_tol;
            obj.so_Max_Steps = options.so_Max_Steps;
            if obj.so_Auto_Bounds; obj.so_circWave_lb = 1; end
            if obj.No_Gauss_Stations > 6; error('AQUINAS Error : Up to 6 Gauss-Legendre stations are supported for numerical integration along the meridian of an element.'); end
            if options.noSimpsonStations < 3; error('AQUINAS Error: No less than 3 through thickness stations may be used for Simpson''s 1/3 integration rule. The number of stations provided must also be odd.'); end
            if mod(options.noSimpsonStations,2)==0
                warning(strcat("AQUINAS Warning: The number of stations for through thickness integration, according to Simpson''s 1/3 rule, must be odd. ",num2str(options.noSimpsonStations+1)," stations will be used instead of ",num2str(options.noSimpsonStations)," through thickness stations."));
                obj.No_Simpson_Stations = options.noSimpsonStations + 1;
            else
                obj.No_Simpson_Stations = options.noSimpsonStations;
            end

            Gauss = obj.Gauss_Legendre_nodes_weights(obj.No_Gauss_Stations);
            obj.Gauss_Nodes = Gauss.nodes; obj.Gauss_Weights = Gauss.weights;

        end

        % Method to compute dLPF increment, according to either the Newton-Raphson or the Arc Length (Riks) method
        function [obj,DELTA,dDELTA_a,dDELTA_1,Dij,zeta] = update_nonlinear_solution(obj,i,j,nnodes,uw,Dijm1,Jim1,DELTA,DELTA_I,dDELTA_R,dDELTA_a,dDELTA_1,dDelta_a_prev)

            if any(isnan(DELTA)) || any(isnan(DELTA_I)) || any(isnan(dDELTA_R)) || any(isnan(dDELTA_a)); obj.flag = -2; return; end
            % Update displacement vector and load factor (for arc-length)
            if i == 1 && j == 1
                if abs(obj.dLPF) > obj.dLPFmax; obj.dLPF = sign(obj.dLPF)*obj.dLPFmax; end
                DELTA = DELTA + obj.dLPF*DELTA_I(1:4*nnodes);
                obj.LPF = obj.LPF + obj.dLPF;
                if strcmp(obj.nonlinear_solver,'ArcLength')
                    dDELTA_a = dDELTA_a + obj.dLPF*DELTA_I(1:4*nnodes);
                    obj.l = sqrt(DELTA(uw)'*DELTA(uw)); % Initial arc-length computation, according to eq. 91 of Teng and Rotter (1989a)
                    obj.lmax = obj.dLPFmax*sqrt(DELTA_I(uw)'*DELTA_I(uw)); % Maximum arc-length to be attempted
                    if obj.l > obj.lmax; obj.l = obj.lmax; end
                end
            else
                if strcmp(obj.nonlinear_solver,'NewtonRaphson')
                    if j == 1
                        DELTA = DELTA + obj.dLPF*DELTA_I(1:4*nnodes);
                        obj.LPF = obj.LPF + obj.dLPF;
                    else
                        DELTA = DELTA + dDELTA_R(1:4*nnodes);
                    end
                elseif strcmp(obj.nonlinear_solver,'ArcLength')
                    if j == 1 && i > 1
                        if ~isnan(Jim1); obj.l = obj.l*sqrt(obj.Jd/Jim1); end % update the arc-length according to eq. 97 of Teng and Rotter (1989a)
                        obj.l = min(obj.l,obj.lmax);
                        obj.dLPF = sign(dDelta_a_prev'*DELTA_I(1:4*nnodes))*obj.l/sqrt(DELTA_I(uw)'*DELTA_I(uw)); % Update load factor increment, according to eq. 92 of Teng and Rotter (1989a)
                    elseif j > 1
                        if strcmp(obj.arclength_method,'Ramm') % LPF increment according to Ramm's method
                            obj.dLPF = - (dDELTA_1(uw)'*dDELTA_R(uw))/(dDELTA_1(uw)'*DELTA_I(uw));
                        elseif strcmp(obj.arclength_method,'Crisfield')
                            A = DELTA_I(uw)'*DELTA_I(uw);
                            B = 2*(dDELTA_a(uw) + dDELTA_R(uw))'*DELTA_I(uw);
                            C = (dDELTA_a(uw) + dDELTA_R(uw))'*(dDELTA_a(uw) + dDELTA_R(uw)) - obj.l^2;
                            Discr = B^2 - 4*A*C;
                            if Discr < 0; obj.flag = -5; Dij = nan; zeta = nan; return; end
                            dLambda1 = (-B + sqrt(Discr))/(2*A);
                            dLambda2 = (-B - sqrt(Discr))/(2*A);
                            % Select the root that maintains a positive angle between the original and accumulated displacements, or is closer to the linear solution if both said angles are positive.
                            if (dDELTA_a(uw)'*(dDELTA_a(uw) + dLambda1*DELTA_I(uw) + dDELTA_R(uw)) > 0) && (dDELTA_a(uw)'*(dDELTA_a(uw) + dLambda2*DELTA_I(uw) + dDELTA_R(uw)) < 0)
                                obj.dLPF = dLambda1;
                            elseif (dDELTA_a(uw)'*(dDELTA_a(uw) + dLambda1*DELTA_I(uw) + dDELTA_R(uw)) < 0) && (dDELTA_a(uw)'*(dDELTA_a(uw) + dLambda2*DELTA_I(uw) + dDELTA_R(uw)) > 0)
                                obj.dLPF = dLambda2;
                            elseif (dDELTA_a(uw)'*(dDELTA_a(uw) + dLambda1*DELTA_I(uw) + dDELTA_R(uw)) > 0) && (dDELTA_a(uw)'*(dDELTA_a(uw) + dLambda2*DELTA_I(uw) + dDELTA_R(uw)) > 0)
                                obj.dLPF = (-C/B);
                            else
                                obj.flag = -5; Dij = nan; zeta = nan;
                                return;
                            end
                        end
                    end
                    dDELTA_a = dDELTA_a + obj.dLPF*DELTA_I(1:4*nnodes) + dDELTA_R(1:4*nnodes);
                    DELTA = DELTA + obj.dLPF*DELTA_I(1:4*nnodes) + dDELTA_R(1:4*nnodes);
                    obj.LPF = obj.LPF + obj.dLPF;
                end
            end
            if strcmp(obj.nonlinear_solver,'ArcLength') && (j==1); dDELTA_1 = dDELTA_a; end
            % Compute nodal displacement magnitudes corrections vector, to be used in the convergence criterion
            Dij = DELTA(uw(1:2:end)).*DELTA(uw(1:2:end)) + DELTA(uw(2:2:end)).*DELTA(uw(2:2:end)); Dij((Dij.*Dij) < eps) = nan; % Ensure division with zero is avoided
            zeta = 1.0 - Dijm1./Dij; zeta(isnan(zeta)) = 0.0; % nan values correspond to division with zero from the Dij vector

        end

        % Check for termination of a nonlinear analysis based on the development of the convergence indicator plot, as presented in
        % Sadowski et al. 'A computational strategy to establish algebraic parameters for the Reference Resistance Design of metal shell structures'
        % and based on the work of Doerich and Rotter 'Accurate determination of plastic collapse loads from finite element analyses'
        % The following function is a matlab implementation of the fortran subroutine presented in the Appendix of the article by Sadowski et al (2017)
        function [obj,cip_check] = convergence_indicator_plot_termination_check(obj,DOF_TRACKER)

            % Preliminaries
            cip_check = false;
            all_ms = abs(DOF_TRACKER.LPF(2:end)./DOF_TRACKER.DOF(2:end)); % x axis of modified southwell plot
            cv = std(all_ms)/mean(all_ms); % coefficient of variation of modified southwell x axis points
            if cv > obj.cip.tol_cv && length(DOF_TRACKER.LPF) > 3
                 xbar = mean(all_ms(end-2:end)); ybar = mean(DOF_TRACKER.LPF(end-2:end));
                 mn1 = (all_ms(end) - xbar)*(DOF_TRACKER.LPF(end) - ybar); md1 = (all_ms(end) - xbar)^2;
                 mn2 = (all_ms(end-1) - xbar)*(DOF_TRACKER.LPF(end-1) - ybar); md2 = (all_ms(end-1) - xbar)^2;
                 mn3 = (all_ms(end-2) - xbar)*(DOF_TRACKER.LPF(end-2) - ybar); md3 = (all_ms(end-2) - xbar)^2;
                 m = (mn1 + mn2 + mn3)/(md1 + md2 + md3); obj.cip.counter = obj.cip.counter + 1; obj.cip.all_pms(obj.cip.counter) = ybar - m*xbar;
                 obj.cip.all_omega_bar(obj.cip.counter) = (obj.cip.all_pms(obj.cip.counter) - DOF_TRACKER.LPF(end))/DOF_TRACKER.LPF(end);
                 if obj.cip.counter > 1
                    sumx = sum(obj.cip.all_omega_bar); sumy = sum(obj.cip.all_pms); sumxy = sum(obj.cip.all_omega_bar.*obj.cip.all_pms); sumxx = sum(obj.cip.all_omega_bar.^2);
                    curr_cip_rpl = (sumy*sumxx - sumx*sumxy)/(obj.cip.counter*sumxx - sumx*sumx);
                    change_cip_rpl = abs((curr_cip_rpl - obj.cip.prev_cip_rpl)/obj.cip.prev_cip_rpl);
                 end
            end
            % Check for the change of the Rpl projection of the CIP plot to determine whether it is below the tolerance threshold
            if length(DOF_TRACKER.LPF) > 3 && obj.cip.counter > 1
                if change_cip_rpl >= obj.cip.tol_cip_rpl
                    obj.cip.prev_cip_rpl = curr_cip_rpl;
                else
                    obj.cip.Rpl = curr_cip_rpl;
                    cip_check = true;
                end
            end

        end

        % A method that determines whether a nonlinear analysis is to be terminated, depending on the termination conditions provided and the progress of the solution
        function [obj,termination_cause] = nonlinear_analysis_termination(obj,i,j,step_converged,increment_attempts,zeta,RESIDUAL,EIGVALUES,DOF_TRACKERS,GMNAstep)

            % Criteria that correspond to an unexpected termination of an analysis are given a negative flag, while a termination condition being met leads to a positive flag indicator
            if increment_attempts > obj.max_attempts
                termination_cause = "Analysis aborted. Number of allowed attempts at a load step has been exceeded.";
                if obj.check_axisym_stability && obj.axisym_path_eig < 0
                    termination_cause = strcat(termination_cause," Negative eigenvalue found, loss of stability of axisymmetric prebuckling path occurred during last uncorverged attempt at an increment.");
                end
                if j==obj.Jmax
                    obj.flag = -1;
                    termination_cause =  strcat(termination_cause," Maximum number of iterations Jmax for achieving convergence of a load step exhausted.");
                elseif (j > 1 && norm(zeta) > obj.zeta_d) || any(isnan(RESIDUAL))
                    obj.flag = -2;
                    termination_cause = strcat(termination_cause," Divergence has occurred.");
                elseif any(ismember(obj.termination_conditions,'C')) && (obj.LPF >= obj.check_bifurcation_after_LPF && min(EIGVALUES) < (1 - obj.ksi))
                    obj.flag = -3;
                    termination_cause = strcat(termination_cause," Critical eigenvalue is below the 1-ksi threshold.");
                elseif obj.LPF >= obj.check_bifurcation_after_LPF && ((strcmp(obj.simultaneous_bifurcation_treatment,'divideKsiBy10') && min(EIGVALUES) > (1 - obj.ksi)))
                    obj.flag = -4;
                    termination_cause = strcat(termination_cause," More than one eigenvalue exist within the 1-ksi < lambda < 1+ksi range.");
                elseif obj.flag == -5
                    termination_cause = strcat(termination_cause," Both roots obtained during the arc-length constraint application, for updating the load factor increment, produce a negative angle between the original and updated accumulated displacements.");
                elseif obj.flag == -6
                    termination_cause = strcat(termination_cause," Effective strain increment at least one material point is over the maximum allowed value of max_epsilon_bar.");
                end
            elseif step_converged && any(ismember(obj.termination_conditions,'A')) && ~isempty(obj.no_max_steps) && i >= obj.no_max_steps
                obj.flag = 1;
                termination_cause = "Termination condition met. Reached maximum number of load steps.";
            elseif step_converged && any(ismember(obj.termination_conditions,'B')) && obj.LPF > obj.LPFmax
                obj.flag = 2;
                termination_cause = "Termination condition met. Maximum LPF to be considered has been exceeded.";
            elseif step_converged && any(ismember(obj.termination_conditions,'C')) &&  ~isempty(GMNAstep{i}.BuckledIntoMode)
                obj.flag = 3;
                termination_cause = "Termination condition met. Bifurcation has occurred.";
            elseif step_converged && any(ismember(obj.termination_conditions,'D')) &&  (i > 1 && GMNAstep{end}.LPF < GMNAstep{end-1}.LPF)
                obj.flag = 4;
                termination_cause = "Termination condition met. Limit point detected.";
            elseif step_converged && any(ismember(obj.termination_conditions,'E'))
                [obj,cip_check] = obj.convergence_indicator_plot_termination_check(DOF_TRACKERS(obj.cip.dof_tracker));
                if cip_check
                    obj.flag = 5;
                    termination_cause = "Termination condition met. Convergence indicator plot prediction for Rpl obtained.";
                else
                    termination_cause = "Unterminated analysis.";
                    obj.flag = 0;
                end
            else
                termination_cause = "Unterminated analysis.";
                obj.flag = 0;
            end

        end

    end

    methods(Static)

        function Gauss = Gauss_Legendre_nodes_weights(noGaussStations)

            switch noGaussStations
                case 1
                    Gauss.nodes   = 0.0;
                    Gauss.weights = 2.0;
                case 2
                    Gauss.nodes   = [-0.577350269189626 0.577350269189626];
                    Gauss.weights = [ 1.0               1.0              ];
                case 3
                    Gauss.nodes   = [-0.774596669241483 0.0               0.774596669241483];
                    Gauss.weights = [ 0.555555555555556 0.888888888888889 0.555555555555556];
                case 4
                    Gauss.nodes   = [-0.861136311594053 -0.339981043584856 0.339981043584856 0.861136311594053];
                    Gauss.weights = [ 0.347854845137454  0.652145154862546 0.652145154862546 0.347854845137454];
                case 5
                    Gauss.nodes   = [-0.906179845938664 -0.538469310105683 0.0               0.538469310105683 0.906179845938664];
                    Gauss.weights = [ 0.236926885056189  0.478628670499366 0.568888888888889 0.478628670499366 0.236926885056189];
                case 6
                    Gauss.nodes   = [-0.9324695142031521 -0.6612093864662645 -0.2386191860831969 0.2386191860831969 0.6612093864662645 0.9324695142031521];
                    Gauss.weights = [ 0.1713244923791704  0.3607615730481386  0.4679139345726910 0.4679139345726910 0.3607615730481386 0.1713244923791704];
                case 7
                    Gauss.nodes   = [-0.9491079123427585 -0.7415311855993945 -0.4058451513773972 0.0                0.4058451513773972 0.7415311855993945 0.9491079123427585];
                    Gauss.weights = [ 0.1294849661688697  0.2797053914892766  0.3818300505051189 0.4179591836734694 0.3818300505051189 0.2797053914892766 0.1294849661688697];
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