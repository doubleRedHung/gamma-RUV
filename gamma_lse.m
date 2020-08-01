% gamma-LSE
%
% ----- INPUT -----
% y: n by 1 response vector
% x: n by p data matrix
% gamma: gamma value
% lam: penalty for stablization 
%
% ----- OUTPUT -----
% b1: robust estimate of coefficient
% s1: robust estimate of sigma^2
% cov_b: estimate of acov(b1) 
% p_val: p-value of b1
% wt: weight for each subject


function [b1, s1, cov_b, p_val, wt] = gamma_lse(y, x, gamma, lam)

    if nargin < 4
        lam = 10^-4;
    end

    [n,p] = size(x);
    d = 1;
    counter = 1;
    DI = lam*eye(p);
    
    b1 = (x'*x + n*DI)\(x'*y);    
    s1 = sum((y-x*b1).^2)/n;
    while d > 10^-3 && counter <= 100
        b0=b1;
        s0=s1;
        wt = exp(-gamma/s1/2*(y-x*b1).^2);
        b1 = (x'*diag(wt)*x + DI)\(x'*diag(wt)*y);
        s1 = (gamma+1)*sum((y-x*b1).^2.*wt)/sum(wt);
        
        d = norm([b1;s1]-[b0;s0])/norm([b0;s0]);
        counter = counter + 1;
    end
    if counter > 100
        display(['gamma-LSE not converge, [r,d] = ',num2str([gamma,d])])
    end
    cov_temp = ((x'*x+DI)\(x'*x))/(x'*x+DI);
    cov_b = (1+gamma^2/(2*gamma+1))^1.5*s1 * cov_temp;
       
    T_temp = b1./sqrt(diag(cov_b));
    p_val = 1-chi2cdf(T_temp.^2,1);
end



%     wt = (2*pi*s1)^-0.5*exp(-1/s1/2*(y-x*b1).^2);
%     wt = wt.^(2*gamma) .* (y-x*b1).^2;
%     cov_b = (1+gamma)^3*(2*pi*s1)^gamma * (x'*x+DI)\(x'*diag(wt)*x)/(x'*x+DI);
    










