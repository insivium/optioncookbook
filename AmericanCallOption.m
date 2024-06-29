function C_american = AmericanCallOption(S, X, T, r, delta, sigma)
    % Parameters
    % S - Spot price of the commodity
    % Critical asset value ?
    % X - Strike price
    % T - Time to maturity
    % r - Risk-free interest rate
    % delta - Dividend yield or convenience yield
    % sigma - Volatility

    % cost of carry
    b = r - delta;
    9
    % Calculate the European Call Price using Black-Scholes formula
    c = BlackScholesCall(S, X, T, r, b, sigma);
    fprintf('European call price: %f\n', c);
    
    % Define the function to find S* ( obtaining complex roots using fzero, using bisection method..)
    func = @(S_star) real(S_star - X - (BlackScholesCall(S_star, X, T, r, b, sigma) + ...
        (1 - exp((b-r)*T) * normcdf(real(d1(S_star, X, T, b, sigma)))) * real(S_star^q2(S_star, r, b, sigma, T))));

    % Initial interval for S* ( We start from exercise)
    a = 0.5 * X;
    b = 1.5 * X;

    % Bisection method to find the root
    tol = 1e-6;
    max_iter = 100;
    iter = 0;
    while iter < max_iter
        S_star = (a + b) / 2;
        f_val = func(S_star);
        f_val_a = func(a);
        
        if abs(f_val) < tol || (b - a) / 2 < tol
            break;
        end
        
        if sign(f_val) == sign(f_val_a)
            a = S_star;
        else
            b = S_star;
        end
        
        iter = iter + 1;
    end
    
    if iter == max_iter
        error('Bisection method did not converge');
    end

    % Debugging: Print S*
    fprintf('Calculated S*: %f\n', S_star);

    % Calculate the call price using the Quadratic Approximation Method
    q2_value = real(q2(S_star, r, b, sigma, T));
    d1_value_S_star = real(d1(S_star, X, T, b, sigma)); % Calculate d1 using S_star
    
    % Debugging: Print intermediate values
    fprintf('q2_value: %f\n', q2_value);
    fprintf('d1_value_S_star: %f\n', d1_value_S_star);
    
    % Calculate the American call price
    C_american = c + real(S_star^q2_value) * (1 - exp((r-delta)*T) * normcdf(d1_value_S_star)) * real((S / S_star)^q2_value);
    fprintf('American call price: %f\n', C_american);
end

function c = BlackScholesCall(S, X, T, r, b, sigma)
    d1_value = real(d1(S, X, T, b, sigma));
    d2_value = real(d1_value - sigma * sqrt(T));
    c = S * exp((b-r) * T) * normcdf(d1_value) - X * exp(-r * T) * normcdf(d2_value);
end

function d1_value = d1(S, X, T, b, sigma)
    d1_value = (log(S / X) + (b + sigma^2 / 2) * T) / (sigma * sqrt(T));
end

function q2_value = q2(S_star, r, b, sigma, T)
    k1 = 2 * r / sigma^2;
    k2 = 2 * b / sigma^2;
    discriminant = sqrt((k2 - 1)^2 + 4 * k1);
    q2_value = 0.5 * (-k2 + 1 + discriminant);
end
