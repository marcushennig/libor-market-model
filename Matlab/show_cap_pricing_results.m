%**************************************************************************
%
% Cap
%
%**************************************************************************
%--------------------------------------------------------------------------
% Cap result
%--------------------------------------------------------------------------

pricingData = './PricingData/cap_pricing.txt';
data = load(pricingData);

strikes = data(:, 1);
prices = data(:, 2);
blackPrices = data(:, 3);

%--------------------------------------------------------------------------
% Figure cap
%--------------------------------------------------------------------------

fig = figure;
hold on
set(gcf, 'Color', [ 1 1 1 ]);

plot(100 * strikes, 100 * prices)
plot(100 * strikes, 100 * blackPrices, 'r--')

title('Price of a cap')

legend({'MC', 'Black-caplet-formula'})
legend('boxoff')

xlabel('cap rate [%]')
ylabel('Price of the cap [$]')
box on

%--------------------------------------------------------------------------
% Export figure
%--------------------------------------------------------------------------

%folder = 'C:\Users\d90789\Documents\Oxford MSc in Mathematical Finance\Modules\Module 5 Advanced Modelling Topics 1\Assignments\LIBOR Market Model\Paper\figures\';
folder = 'C:\figures\';
fileName = 'cap_pricing.eps';
filePath = [folder fileName];

saveas(fig, filePath, 'eps2c');
system(['epstopdf ' filePath]);