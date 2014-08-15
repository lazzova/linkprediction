function [J, gradient] = linkpred_cost (obj, param)
    obj.runProblem(param);
    J = obj.getCost();
    gradient = obj.getGradient();
end