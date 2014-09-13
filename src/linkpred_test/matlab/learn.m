function [J, w] = learn ()
    % set the path to the java classes
    javaaddpath('/home/skyra/matlab/javaclasses');
    clear java;
    
    % set parameters
    g = 2;
    n = 100;
    fs = [2,3];
    s = 0;
    alpha = 0;
    b = 1e-6;
    lambda = 1;
    param = [1, -1, 0.5, -2, 1];
    
    % create an object for the problem
    o = init_linkpred(g, n, fs, s, alpha, b, lambda, param, learningrate, );
    
    % generates random numbers between -3 and 3
    w_temp = -3+rand(f, 1).*(3+3);
   
    % set optimization options 
    options = optimset('GradObj', 'on', 'MaxIter', 100, ...
        'FinDiffRelStep', 0.01);
    
    % create the optimization problem
    problem = createOptimProblem('fmincon', ...
        'objective', @(x)linkpred_cost(o, x), ...
        'x0', w_temp, ...        
        'lb', [-3;-3], ...
        'ub', [3;3], ...
        'options', options);
    
    % create multistart
    ms = MultiStart('UseParallel', 'always', ...
                    'Display', 'iter', ...
                    'StartPointsToRun', 'bounds');
                    %'PlotFcns', @gsplotbestf);
    
    % run parallel, 8 cores used
    % matlabpool('open', 2);
    
    %run the algorithm from 15 different points
    tic
    [w, J] = run(ms, problem, 15);  
    toc
    
    % close the pool
    % matlabpool('close');
    
end
