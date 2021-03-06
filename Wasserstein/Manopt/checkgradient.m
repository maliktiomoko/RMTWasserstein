function checkgradient(problem, x, d)
% Checks the consistency of the cost function and the gradient.
%
% function checkgradient(problem)
% function checkgradient(problem, x)
% function checkgradient(problem, x, d)
%
% checkgradient performs a numerical test to check that the gradient
% defined in the problem structure agrees up to first order with the cost
% function at some point x, along some direction d. The test is based on a
% truncated Taylor series (see online Manopt documentation).
%
% It is also tested that the gradient is indeed a tangent vector.
% 
% Both x and d are optional and will be sampled at random if omitted.
%
% See also: checkdiff checkhessian

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 
%
%   April 3, 2015 (NB):
%       Works with the new StoreDB class system.
%
%   Nov. 1, 2016 (NB):
%       Now calls checkdiff with force_gradient = true, instead of doing an
%       rmfield of problem.diff. This became necessary after getGradient
%       was updated to know how to compute the gradient from directional
%       derivatives.

    
%pick them at random.

x = problem.M.rand();
d = problem.M.randvec(x);


    %% Check that the gradient yields a first order model of the cost.
    
    % Call checkdiff with force_gradient set to true, to force that
    % function to make a gradient call.
    checkdiff(problem, x, d, true);
    title(sprintf(['Gradient check.\nThe slope of the continuous line ' ...
                   'should match that of the dashed\n(reference) line ' ...
                   'over at least a few orders of magnitude for h.']));
    xlabel('h');
    ylabel('Approximation error');
    
    %% Try to check that the gradient is a tangent vector.
    if isfield(problem.M, 'tangent')
        storedb = StoreDB();
        key = storedb.getNewKey();
        grad = getGradient(problem, x, storedb, key);
        pgrad = problem.M.tangent(x, grad);
        residual = problem.M.lincomb(x, 1, grad, -1, pgrad);
        err = problem.M.norm(x, residual);
        fprintf('The residual should be 0, or very close. Residual: %g.\n', err);
        fprintf('If it is far from 0, then the gradient is not in the tangent space.\n');
    else
        fprintf(['Unfortunately, Manopt was unable to verify that the '...
                 'gradient is indeed a tangent vector.\nPlease verify ' ...
                 'this manually or implement the ''tangent'' function ' ...
                 'in your manifold structure.']);
    end

end
