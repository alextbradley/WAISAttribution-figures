function pslr = get_pslr(x, slr, gamma, D, sigma_m, sigma_g, mu)

% Inputs:
% x :       (nx x 1) array
%           target sea level rise values
%
% slr :     (ng x 1) array 
%           sea level rise as a function of gamma (at a given time)
%
% gamma :   (ng x 1) array
%           gamma values which return slr values above
%
% D :       (ng x 1) array
%           errors in melting
%
% sigma_m : scalar
%           error covariance in melting
%
% sigma_g : scalar
%           error covariance in gamma
%
% mu :      scalar
%           hyperparameter of the prior


% Outputs:
%
% pslr :    (nx x 1) array
%           probability of slr as a function of x

pslr = zeros(size(x));


% for each x, compute gamma_x (the values of gamma which return slr of x)
for ix = 1:length(x)
    gamma_x = get_gamma_x(x(ix), slr, gamma);

    for ig = 1:length(gamma_x)

        % for each gamma_x, compute D(gamma_x) and therefore P(mitgcm| gamma_x)
        D_gammax = interp1(gamma,D,gamma_x(ig)); %use linear interpolation to find the value of D(gamma_x);
        P_mit_given_gammax = 1/sqrt(2*pi*sigma_m^2) * exp (-D_gammax .^2 /2/sigma_m^2);
        
        % for each gamma_x, compute P(gamma_x | mu)  
        P_gammax_given_mu =  1/sqrt(2*pi*sigma_g^2)*exp(-(gamma_x - mu).^2 /2/sigma_g^2);

        % take the product and add it
        l = P_mit_given_gammax * P_gammax_given_mu;
       
        pslr(ix) = pslr(ix) + l;
        
    end
end
