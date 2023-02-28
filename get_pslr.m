function pslr = get_pslr(x, slr, gamma, D, sigma_m, sigma_g, mu)

% Inputs:
% x :       (nx x 1) array
%           target sea level rise values (assumes regular grid)
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

% create the fit to P(mit| gamma_x)
yy = 1/sqrt(2*pi*sigma_m^2) * exp (-D .^2 /2/sigma_m^2);
dg = diff(gamma);
dg = dg(1);
yy = yy/(sum(yy)*dg);
yy = reshape(yy, [length(yy),1]); %make column vector
gamma = reshape(gamma, [length(gamma),1]);
f = fit(gamma,yy,'Gauss1');

% for each x, compute gamma_x (the values of gamma which return slr of x)
for ix = 1:length(x)
    gamma_x = get_gamma_x(x(ix), slr, gamma);


    for ig = 1:length(gamma_x)

        % for each gamma_x, compute D(gamma_x) and therefore P(mitgcm| gamma_x)
        %D_gammax = interp1(gamma,D,gamma_x(ig)); %use linear interpolation to find the value of D(gamma_x);
        %P_mit_given_gammax = 1/sqrt(2*pi*sigma_m^2) * exp (-D_gammax .^2 /2/sigma_m^2);
        P_mit_given_gammax = f(gamma_x(ig));


        % for each gamma_x, compute P(gamma_x | mu)  
        P_gammax_given_mu =  1/sqrt(2*pi*sigma_g^2)*exp(-(gamma_x(ig) - mu).^2 /2/sigma_g^2);

        % take the product and add it
        l = P_mit_given_gammax * P_gammax_given_mu;
        
       
         try
        pslr(ix) = pslr(ix) + l;
        catch
            l
        end
        
    end
end

% normalize
dx = diff(x); dx = dx(1); 
pslr = pslr / (sum(pslr) * dx);

