% gamma-RUV
%
% ----- INPUT -----
% X: n by p data matrix
% r: gamma value
% diag_sig: 1 = diagnoalize SIgma when fitting 
%
% ----- OUTPUT -----
% Z: robust standardized data
% mu1: robust estimate of mean vector
% Sig1: robust estimate of covariance matrix
% wt: weight for each subject 

function [Z, mu1, Sig1, wt] = gamma_prewhiten(X, r, diag_sig)

    if nargin < 3
        diag_sig = 0;
    end
    [n,p] = size(X);
    d=1;
    counter=1;
    mu1=mean(X)';    
    Sig1=cov(X);
    D = 10^-8*eye(p);
    
    while d > 10^-4 && counter <= 50
        mu0 = mu1;
        Sig0 = Sig1;
        if diag_sig == 0
            wt=(X-ones(n,1)*mu0')*(Sig0 + D)^-0.5;  
%             wt=(X-ones(n,1)*mu0')*(Sig0)^-0.5;  
        elseif diag_sig == 1
            wt=(X-ones(n,1)*mu0')*diag(diag(Sig0).^-0.5);  
        end
        wt=exp(-r*diag(wt*wt')/2);
        wt=wt/sum(wt);
        mu1=X'*wt;
        Sig1=(r+1)*(X-ones(n,1)*mu0')'*diag(wt)*(X-ones(n,1)*mu0');
       
        d = norm([mu1;Sig1(:)]-[mu0;Sig0(:)])/norm([mu0;Sig0(:)]);
        counter=counter+1;
    end
    Z=(X-ones(n,1)*mu1')*(Sig1 + D)^-0.5;
    
    if counter == 101
        display(['Not converge at ',num2str(r)])
    end
end






