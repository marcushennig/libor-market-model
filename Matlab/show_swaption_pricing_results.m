%--------------------------------------------------------------------------
% Clean up workspace
%--------------------------------------------------------------------------

clear
clc
close all

%**************************************************************************
%
% Caplet
%
%**************************************************************************
%--------------------------------------------------------------------------
% Caplet result
%--------------------------------------------------------------------------

pricingData = './PricingData/swaption_pricing.txt';
data = load(pricingData);

strikes = data(:, 1);
mcPrices = data(:, 2 : end);
prices = mean(mcPrices, 2);
confidence = 0.96;
error = norminv((1+confidence) / 2) * std(mcPrices, 0, 2);

%--------------------------------------------------------------------------
% Figure caplet
%--------------------------------------------------------------------------

fig =figure;
hold on
set(gcf, 'Color', [ 1 1 1 ]);

% Monte-Carlo pricing result 
plot(strikes, prices, '-', 'Color', [ 0.6 0 0.5 ]);
h(1) = plot(strikes(1 : 4 : end), prices(1 : 4 : end), ...
            '.', ...
            'MarkerSize',  15, ...
            'Color', [ 0.6 0 0.5 ]);

% confidence interval 
h(2) = plot(strikes, prices - error, '--', 'Color', [ 0.8 0.4 0.8 ]);
plot(strikes, prices + error, '--',  'Color', [ 0.8 0.4 0.8 ])


title('Price of a swaption')

legend(h, {'Monte-Carlo', ...
           '96% confidence interval'})
legend('boxoff')

xlabel('strike [%]')
ylabel('Price [$]')
box on

%--------------------------------------------------------------------------
% Export figure
%--------------------------------------------------------------------------

%folder = 'C:\Users\d90789\Documents\Oxford MSc in Mathematical Finance\Modules\Module 5 Advanced Modelling Topics 1\Assignments\LIBOR Market Model\Paper\figures\';
folder = 'C:\figures\';
fileName = 'swaption_pricing.eps';
filePath = [folder fileName];

saveas(fig, filePath, 'eps2c');
system(['epstopdf ' filePath]);