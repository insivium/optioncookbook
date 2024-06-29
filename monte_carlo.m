function [call_price, payoff_variance, confidence_interval,call_bs] = monte_carlo(S0, K, r, T, sigma, N)
    % S0: initial stock price
    % K: strike price
    % r: risk-free interest rate
    % T: time to maturity (in years)
    % sigma: volatility
    % N: number of simulations

    % Simulate end-of-period stock prices
    Z = randn(N, 1); % Generate N standard normal random variables
    ST = S0 * exp((r - 0.5 * sigma^2) * T + sigma * sqrt(T) * Z); % Stock price at maturity

    % Calculate the payoff for each simulation
    payoffs = max(ST - K, 0);

    % Discount the payoffs to present value
    discounted_payoffs = exp(-r * T) * payoffs;

    % Estimate the call price
    call_price = mean(discounted_payoffs);

    % Calculate the sample variance of the discounted payoffs
    a_M = call_price;
    b_M = var(discounted_payoffs, 1); %  biased estimator with 1/N

    payoff_variance = b_M;

    % Calculate the standard error
    SE = sqrt(b_M / N);

    % Determine the critical value for a 95% confidence interval
    z = 1.96;

    % Construct the confidence interval
    confidence_interval = [call_price - z * SE, call_price + z * SE];

     % Black-Scholes option pricer
    d1 = (log(S0 / K) + (r + 0.5 * sigma^2) * T) / (sigma * sqrt(T));
    d2 = d1 - sigma * sqrt(T);
    call_bs = S0 * normcdf(d1) - K * exp(-r * T) * normcdf(d2);

    fprintf('The estimated European call option price is: %.4f\n', call_price);
fprintf('The variance of the discounted payoffs is: %.4f\n', payoff_variance);
fprintf('The 95%% confidence interval is: [%.4f, %.4f]\n', confidence_interval(1), confidence_interval(2))
fprintf('The Black-Scholes call option is: %.4f\n', call_bs);
end


