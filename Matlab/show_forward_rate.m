%--------------------------------------------------------------------------
% Clean up workspace
%--------------------------------------------------------------------------

clear
clc
close all

%**************************************************************************
%
% Foward rate time series 
%
%**************************************************************************
%--------------------------------------------------------------------------
% Forward rate data
%--------------------------------------------------------------------------

forwardRateFile = './PricingData/lmm_forward_rates.txt';
data = load(forwardRateFile);
maturity = data(1,2 : end);
t = data(2:end, 1);
j = 1 : size(data, 2) - 1;
[Maturity, Time ] = meshgrid(maturity, t);
forwardRates = data(2:end, 2:end);

%--------------------------------------------------------------------------
% Forward rate data
%--------------------------------------------------------------------------

fig = figure;

surf(Time, Maturity, 100 * forwardRates)
set(gca, 'YDir', 'reverse')
title('6M forward LIBOR time series')

xlabel('time [years]')
ylabel('maturity [years]')
zlabel('6M forward libor rate [%/year]')
box on 

%--------------------------------------------------------------------------
% Export figure
%--------------------------------------------------------------------------

%folder = 'C:\Users\d90789\Documents\Oxford MSc in Mathematical Finance\Modules\Module 5 Advanced Modelling Topics 1\Assignments\LIBOR Market Model\Paper\figures\';
folder = 'C:\figures\';
fileName = 'forward_rate_time_series.eps';
filePath = [folder fileName];

saveas(fig, filePath, 'eps2c');
system(['epstopdf ' filePath]);

%**************************************************************************
%
% Yield curve 
%
%**************************************************************************

%--------------------------------------------------------------------------
% yield curve data
%--------------------------------------------------------------------------

yielCurveFile = './PricingData/yield_curve.txt';
data = load(yielCurveFile);
t = [0; data(:, 1)];
yield = [0; data(:, 2)];

%--------------------------------------------------------------------------
% yield curve figure
%--------------------------------------------------------------------------

fig = figure;
hold on
title('Yield curve')
plot(t, 100 * yield, '-', 'Color', [ 0.6 0 0.5 ]);
plot(t, 100 * yield, '.', 'MarkerSize',  15, 'Color', [ 0.6 0 0.5 ]);
xlabel('time [years]')
ylabel('Price [%/year]')
box on

%--------------------------------------------------------------------------
% Export figure
%--------------------------------------------------------------------------

%folder = 'C:\Users\d90789\Documents\Oxford MSc in Mathematical Finance\Modules\Module 5 Advanced Modelling Topics 1\Assignments\LIBOR Market Model\Paper\figures\';
folder = 'C:\figures\';
fileName = 'yield_curve.eps';
filePath = [folder fileName];

saveas(fig, filePath, 'eps2c');
system(['epstopdf ' filePath]);