%--------------------------------------------------------------------------
% Clean up workspace
%--------------------------------------------------------------------------

clear
clc
close all

%--------------------------------------------------------------------------
% Load external C#-Pricing library
%--------------------------------------------------------------------------

folder = 'C:\Users\d90789\Documents\Oxford MSc in Mathematical Finance\Modules\C# Projects\PricingEngine\PricingLibrary\bin\Release\';
% If .NET assembly is modified MATLAB has to be restarted as there is no 
% way yet to unload the .NET assembly
NET.addAssembly( [ folder 'PricingLibrary.dll' ]);

%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------

% Consider a ten-year semi-annual European swaption (10y x 5y), exercisable on each 
% reset date, starting on 1-October-2020
% settle or valuation date
settleDate = datetime( 2007, 12, 15 );

% 1 EUR 
nominal = 1;
  
% Exercise date of the option = start date of swap
exerciseDate = settleDate + years(1);

% Maturity of the underlying swap
maturityDate = exerciseDate + years(10);
            
% Strike (fixed interest rate payed by fixed leg) 
strike = 0.045;
            
% Tenor of the underlying swap in units of [month]
% = floating leg pays 6M-Libor 
tenor = 6;

%--------------------------------------------------------------------------
% Zero curve observable in the market
%--------------------------------------------------------------------------

curveTimes = [ 1 : 5 7 10 20 ]';
zeroRates = [.01 .018 .024 .029 .033 .034 .035 .034]';
curveDates = daysadd(settleDate, 360 * curveTimes, 1 );

%--------------------------------------------------------------------------
% Volatility
%--------------------------------------------------------------------------

% Use the volatility term structure populaized by Riccardo Rebonato
% sigma_i(t) =  (a * tau + b) * Math.Exp(-c * tau) + d; with tau = T[i] - t
% t is assumed to be in units of [years]
volatiltyParameters = [ 0.3, -0.02, 0.7, 0.14 ];

% Instantaneous correlations between the Brownian motions 
% and hence the forward LIBOR rates
% rho(i,j) =   exp(-beta * |T[i] - T[j]|)
beta = 0.08;

%--------------------------------------------------------------------------
% Initialize LIBOR market model using matlab library
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Convert matlab variables to .NET variables
%--------------------------------------------------------------------------

exerciseDateNet =  System.DateTime( year(exerciseDate), month(exerciseDate), day(exerciseDate));
maturityDateNet = System.DateTime( year(maturityDate), month(maturityDate), day(maturityDate));

settleDateNet = System.DateTime( year(settleDate), month(settleDate), day(settleDate));

volatiltyParametersNet = NET.convertArray(volatiltyParameters, 'System.Double', length(volatiltyParameters));
yieldsNet = NET.convertArray(zeroRates, 'System.Double', length(zeroRates));
curveDatesNet = NET.createArray('System.DateTime', length(curveDates));

for i = 1 : length(curveDates)
    curveDatesNet(i) = System.DateTime(year(curveDates(i)), ...
                                          month(curveDates(i)), ...
                                          day(curveDates(i)));
end

%--------------------------------------------------------------------------
% Build swaption (.NET object)
%--------------------------------------------------------------------------   

swaption =  PricingLibrary.FinancialInstruments.Swaption(settleDateNet, ...
                                                         exerciseDateNet, ...
                                                         maturityDateNet, ...
                                                         tenor, strike, nominal);

%--------------------------------------------------------------------------
% Build zero curve (.NET object)
%--------------------------------------------------------------------------


% YieldCurve(DateTime settleDate,  DateTime[] maturityDates, double[] yield)
zeroCurve = PricingLibrary.MarketData.YieldCurve(settleDateNet, curveDatesNet, yieldsNet);

%--------------------------------------------------------------------------
% Build LIBOR market model (.NET object)
%--------------------------------------------------------------------------

seed = 41;
liborMarketModel = PricingLibrary.MarketModels.LiborMarketModel(zeroCurve, ...
                                                                beta, ...
                                                                volatiltyParametersNet, ... 
                                                                swaption.TenorStructure, ...
                                                                seed);
                                                               
%--------------------------------------------------------------------------
% Pricing of payer-swaption for different strikes using the LIBOR market 
% model. The buyer pays fixed interest rate and receives LIBOR
%--------------------------------------------------------------------------

numberOfMonteCarloPaths = 100;
strike = 0 : 0.0025 : 0.1;
price = zeros(size(strike));

for i = 1 : numberOfMonteCarloPaths
    % Simulate forward LIBOR rate paths using Monto carlo simulation up to the
    % exercise date of the swaption 
    paths = liborMarketModel.SimulateForwardRatePathes(swaption.ExerciseDate);
    
    for j = 1 : length(strike)
        swaption.Strike = strike(j);
        price(j) = price(j) + swaption.DiscountedPayoff(paths);
    end
end
price = price / numberOfMonteCarloPaths;

%--------------------------------------------------------------------------
% Figure
%--------------------------------------------------------------------------

figure
set(gcf, 'Color', [ 1 1 1 ]);
plot(100 * strike, 100 * price)
title('Price of a european payer swaption')
xlabel('Fixed Interest rate (strike) [%]')
ylabel('Price of swaption [%]')                                                               