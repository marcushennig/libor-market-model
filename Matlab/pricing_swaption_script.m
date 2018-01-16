%--------------------------------------------------------------------------
% Clean up workspace
%--------------------------------------------------------------------------

clear
clc
close all

%--------------------------------------------------------------------------
% Construct LIBOR Market Model Object 
%--------------------------------------------------------------------------

plotFigures = 0;

%--------------------------------------------------------------------------
% Tenor structure  
%--------------------------------------------------------------------------

% Consider a ten-year semi-annual European swaption (10y x 5y), exercisable on each 
% reset date, starting on 1-October-2020
% settle or valuation date
settleDate = datetime( 2007, 12, 15 );

% Exercise date of the option = start date of swap
exerciseDate = settleDate + years(1);

% Maturity of the underlying swap
maturityDate = exerciseDate + years(10);
            
% Tenor of the underlying swap in units of [years]
% = floating leg pays 6M-Libor 
tenor = 6;

% tenorStructure
T = datenum(tenorStructure(exerciseDate, maturityDate, tenor));

% number of Libor rates L(T[i], T[i+1])
nRates = length(T) - 1;

%--------------------------------------------------------------------------
% Discount curve 
%--------------------------------------------------------------------------

curveTimes = [ 1 : 5 7 10 20 ]';
zeroRates = [.01 .018 .024 .029 .033 .034 .035 .034]';
curveDates = daysadd(settleDate, 360 * curveTimes, 1 );

curve = IRDataCurve('Zero', settleDate, curveDates, zeroRates);

%--------------------------------------------------------------------------
% plot discount curve 
%--------------------------------------------------------------------------

if(plotFigures)

    figure
    set(gcf, 'Color', [ 1 1 1 ] )
    plot( curveDates, curve.getDiscountFactors(curveDates))
    datetick
    title( ['Discount factors for ' datestr(settleDate)]);
    
end

%--------------------------------------------------------------------------
% volatility function
%--------------------------------------------------------------------------

volatilityFunc = @( a, t ) ( a(1) * t + a(2)) .* exp(-a(3) * t) + a(4);
volatiltyParameters = [ 0.3, -0.02, 0.7, 0.14 ];
  
volatilityFunction( 1 : nRates ) = { @(t) volatilityFunc(volatiltyParameters, t) };

% Plot the volatility function
if(plotFigures)
    figure
    set(gcf, 'Color', [ 1 1 1 ] )
    fplot( volatilityFunction{1}, [0 20] )
    title('Volatility Function')
end

%--------------------------------------------------------------------------
% Correlation between the Brownian motions driving the LIBOR rates
%--------------------------------------------------------------------------
beta = .08;

% use act/360 
corrFuncunctoin = @( i, j, beta) exp(-beta * abs( yearfrac(T(i), T(j), 2 )));

correlationMatrix = corrFuncunctoin( meshgrid( 1 : nRates )', ...
                                     meshgrid( 1 : nRates ), beta);
% disp('Correlation Matrix')
% fprintf([repmat('%1.3f ', 1, length(correlationMatrix)) ' \n'], correlationMatrix)

%--------------------------------------------------------------------------
% Swaption parameter
%--------------------------------------------------------------------------

% how may time steps for simulation
nPeriods = 5;

% time step between simulation
deltaTime = 1;
% dates where states of LMM are calcuted, t[i] in [0, T]
simDates = daysadd(settleDate, 360 * deltaTime * ( 0 :  nPeriods),  1 );
% compute time steps t[i+1]-t[i]
simTimes = diff(yearfrac(simDates(1), simDates) );

% Tenor structure of the contract in units of years []
% Ti - T0, contradicts with exercise and maturity => 5 years
% number of forward rates is 10
Tenor = ( 1 : 10 )';
% year fractopm T[i]-T[0]
alpha = yearfrac(T(1), T)';

% For 1 year periods and an evenly spaced tenor, the exercise row will be
% the sixth row and the swaption maturity will be the 5th column
expirationRow = nPeriods + 1;
maturityCol = 17;

%--------------------------------------------------------------------------
% Pricing simulation
%-------------------------------------------------------------------------
%  The default is 2, meaning forward rates are spaced at 0, .5, 1, 1.5, and so on. 
% Possible values for Period are: 1, 2, 4, and 12.
period = 1;
liborMarketModel = LiborMarketModel(curve, volatilityFunction, correlationMatrix, 'Period', period);

% simulates 10000 future zero curve paths
nTrials = 1000;

% result is [nPeriods+1]-by-[nTenors]-by-[nTrials] matrix of simulated zero-rate term structures.
[zeroRate, forwardRate] = liborMarketModel.simTermStructs(nPeriods, 'nTrials', nTrials);

% forward rates at maturity
F = squeeze(forwardRate(end,:,:));

% annual yield, remove singelton dimensions using squeeze
Y = squeeze(zeroRate(end,:,:));

% need discount factor at excercise date t0 of the swaption;
% calculates the discount factor B(T0,Tk) = exp(-(Tk-tT0)*Y(t0,Tk))
B = exp(-bsxfun(@times,  yearfrac(T(1), T)', Y));

strike = 0 : 0.0025 : 0.1;
price = zeros(size(strike));
for i = 1 : length(strike)
    price(i) = mean(realizedDF .* payoffValue );
end

% % Swap rate at t = T0: S(T0) = (B(T0,T0)-B(T0,Tn))/sum(ai*B(T0,Ti+1)), what
% % the tenor ai=const???
% swapRate = (1 - discountFactors(expirationRow, maturityCol, : ))./sum(bsxfun(@times,1, discountFactors(expirationRow,1 : maturityCol,:)));
% 
% nominal = 1;
% strike = 0 : 0.0025 : 0.1;
% price = zeros(size(strike));
% for i = 1 : length(strike)
%     payoffValue = nominal * max( swapRate - strike(i), 0).*sum(bsxfun(@times,1,discountFactors(expirationRow,1:maturityCol,:)));
%     realizedDF = prod(exp(bsxfun(@times,-zeroRates( 2 : expirationRow, 1, : ), simTimes(1:expirationRow-1))),1);
%     
%     price(i) = mean(realizedDF .* payoffValue );
% end
% 
% %--------------------------------------------------------------------------
% % Figure
% %--------------------------------------------------------------------------
% 
% figure
% set(gcf, 'Color', [ 1 1 1 ]);
% 
% plot(100 * strike, 100 * price)
% title('Price of a european payer swaption')
% 
% xlabel('Fixed Interest rate (strike) [%]')
% ylabel('Price of swaption [%]') 
% 
% %--------------------------------------------------------------------------
% % Evolution of zero curve
% %--------------------------------------------------------------------------
%  
% % trailIndex = 1;
% % 
% % if(plotFigures)
% %     
% %     figure
% %     tmpPlotData = zeroRates( :, :, trailIndex);
% %     tmpPlotData(tmpPlotData == 0) = NaN;
% % 
% %     surf(Tenor, simDates, tmpPlotData)
% % 
% % 
% %     title(['Evolution of the Zero Curve for Trial:' num2str(trailIndex) ' of LIBOR Market Model'])
% %     xlabel('Tenor (Years)')
% % end