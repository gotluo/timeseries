function [logS logFq, Hq, tq, Dq, aq, fq] = mfdfa(X, m, q)
% MFDFA  Multifractal detrended fluctuation analysis.
% 
%   Input:
%       X (vector): time series.
%       m (single value): polynomial order for fitting local trend, e.g., 1  for linear trend.
%       q (either single value or vector): scaling moment.
%   Output:
%       logS (vector): log(s), where s is a vector of scales.
%       logFq (vector): log(F(q)), where F is a vector of fluctuation  functions.
%
% This script is implimented by by Dr. Yu Zhou and Dr. Ming Luo (luom38mail.sysu.edu.cn). 
%
% For those using this script, your are recommended to cite the following 
%   papers in your work: 
%     Luo, M., Leung, Y., Zhou, Y., Zhang, W., 2015. Scaling behaviors of 
%         global sea surface temperature. Journal of Climate, 28: 3122-3132,
%         https://doi.org/10.1175/JCLI-D-13-00743.1
%
%

warning off;

N = length(X); % length of time series
num_s = 30; % number of scales
s1 = zeros(num_s, 1); % scales    

X_ave = nanmean(X);
Y = zeros(1, N);
for i = 1 : N
    Y(i) = nansum(X(1 : i) - X_ave);
end
    
F_q = zeros(num_s, length(q));

% s_min = m + 3;
% s_max = round(N / 3);
s_min = 12;
s_max = round(N / 4); 

scales = linspace(log(s_min), log(s_max), num_s);
for i = 1 : num_s
    s = round(exp(scales(i)));
    s1(i) = log10(s);
    N_s = fix(N / s);
      
    F = zeros(1, 2 * N_s);
    for v = 1 : N_s
        p = polyfit([1 : s], Y((v - 1) * s + 1 : v * s), m);
        y = polyval(p, [1 : s]);
        F(v) = nansum((Y((v - 1) * s + 1 : v * s) - y) .^ 2) / s;
    end
    for v = N_s + 1 : 2 * N_s
        p = polyfit([1 : s], Y(N - (v - N_s) * s + 1 : N - (v - N_s - 1) * s) ,m);
        y = polyval(p, [1 : s]);
        F(v) = nansum((Y(N - (v - N_s) * s + 1 : N - (v - N_s - 1) * s) - y) .^ 2) / s;
    end
    
    for nq = 1:length(q);
        if q(nq) == 0
            F_q(i, nq) = exp(nanmean(log(F)) / 2);
        else
            F_q(i, nq) = (nansum(F .^(q(nq) / 2)) / (2 * N_s)) ^ (1 / q(nq));
            
            % when q==1
            F_q1(i, 1) = (nansum(F .^ (1 / 2)) / (2 * N_s)) ^ (1 / 1);
        end
    end
end

logS = s1; 
logFq = log10(F_q);
logFq1 = log10(F_q1);


%% multifratality
if(1)
    logS = unique(logS, 'rows', 'stable');
    logFq = unique(logFq, 'rows', 'stable');
    logFq1 = unique(logFq1, 'rows', 'stable');
    
    %% H(q)
    crossover = [1 25];
    for nq = 1 : length(q)
        C = polyfit(logS(crossover(1): crossover(2)), logFq(crossover(1): crossover(2), nq), 1);
        Hq(nq) = C(1);
        qRegLine{nq} = polyval(C, logS(crossover(1): crossover(2)));
    end
    % when q==1
    C = polyfit(logS(crossover(1): crossover(2)), logFq1(crossover(1): crossover(2)), 1);
    H1 = C(1);
    
%     plot(q, Hq, 'o');
%     xlabel('{\itq}','fontweight','normal','fontsize',18);
%     ylabel('{\ith}({\itq})','fontweight','normal','fontsize',18);
%     set(gca, 'fontsize', 18);
%     grid on;
    
    %% t(q)
    tq=Hq.*q-1;
    % figure;
    % plot(q, tq, 'o');
    % xlabel('{\itq}','fontweight','normal','fontsize', 18);
    % ylabel('{\it\tau}({\itq})','fontweight','normal','fontsize',18);
    % set(gca, 'fontsize', 18);
    % grid on;
    
    %% t(q) by Zhou et al (2011)
    tq = Hq .* q - q * (H1 - 1) - 1;
    % figure;
    % plot(q, tq, 'o');
    % xlabel('{\itq}','fontweight','normal','fontsize', 18);
    % ylabel('{\it\tau}({\itq})','fontweight','normal','fontsize',18);
    % set(gca, 'fontsize', 18);
    % grid on;
    
    %% f(q)
    aq = diff(tq) ./ diff(q);
    fq = q(1 : end - 1) .* aq - tq(1 : end - 1);
    % figure;
%     plot(aq, fq, '-', 'linewidth', 2);
    % xlabel('{\it\alpha}','fontweight','normal','fontsize', 18);
    % ylabel('{\itf}({\it\alpha})','fontweight','normal','fontsize', 18);
    % set(gca, 'fontsize', 18);
    % grid on;
    
    %% D(q)
    Dq = tq ./ (q - 1);
    % figure;
    % plot(q, Dq, 'o');
    % xlabel('{\itq}','fontweight','normal','fontsize', 18);
    % ylabel('{\itD}({\itq})','fontweight','normal','fontsize', 18);
    % set(gca, 'fontsize', 18);
    % grid on;
end 

