%--------------------------------------------------------------------------
% Clean up workspace
%--------------------------------------------------------------------------

clear
clc
close all

%--------------------------------------------------------------------------
% Discount curve 
%--------------------------------------------------------------------------
% 20-year IR curve 

curveTimes = [ 0 1 : 5 7 10 20 ]';
zeroRates = [0 .01 .018 .024 .029 .033 .034 .035 .036]';

fig = figure;
hold on 
plot(curveTimes, 100 * zeroRates);
plot(curveTimes, 100 * zeroRates, '.', 'MarkerSize', 20);

%curveDates = daysadd(settleDate, 360 * curveTimes, 1 );

%curve = IRDataCurve('Zero', settleDate, curveDates, zeroRates);
%--------------------------------------------------------------------------
% Export figure
%--------------------------------------------------------------------------

%folder = 'C:\Users\d90789\Documents\Oxford MSc in Mathematical Finance\Modules\Module 5 Advanced Modelling Topics 1\Assignments\LIBOR Market Model\Paper\figures\';
folder = 'C:\figures\';
fileName = 'yield_curve.eps';
filePath = [folder fileName];

saveas(fig, filePath, 'eps2c');
system(['epstopdf ' filePath])