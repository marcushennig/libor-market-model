using System;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Random;

namespace PricingLibrary.FinancialInstruments
{
    public class AmericanOption
    {
        #region Properties
        #endregion Properties

        #region Methods
        /// <summary>
        /// Based on code originally downloaded from 
        /// http://www.mathworks.com/matlabcentral/fileexchange/
        ///        16476-pricing-american-options/content/AmericanOptLSM.m
        ///
        /// AmericanPutLSM - Price an american put option via Longstaff-Schwartz Method
        ///
        /// Inputs:
        ///
        ///   S0      Initial asset price
        ///   K       Strike Price
        ///   r       Interest rate
        ///   T       Time to maturity of option
        ///   sigma   Volatility of underlying asset
        ///   N       Number of timesteps
        ///   M       Number of paths 
        /// </summary>
        public static void TestFunctions()
        {
            var S0  = 1.0;   // initial asset price
            var K   = 1.0;   // strike
            var r   = 0.05;  // risk-free interest rate
            var T   = 1.0;   // time to maturity
            var sig = 0.2;   // volatility
            var nTimes   = 64;    // number of timesteps
            var nPaths   = 100000;   // number of paths
            var dt = T/nTimes;

            const int randomSeed = 42; 
            var randomGenerator = new MersenneTwister(randomSeed);
            // Normally distributed independent random numbers 
      
            // dW = sqrt(dt)*randn(N,M);                                    
            var dW = Matrix<double>.Build
                .Random(nTimes, nPaths, Normal.WithMeanStdDev(0, 1))
                .Multiply(Math.Sqrt(dt));

            // dX = exp((r-sig^2/2)*dt+sig*dW);
            var dX = dW.Multiply(sig)
                       .Add((r - 0.5 * sig * sig) * dt)
                       .PointwiseExp();
            
            // Y = repmat(S0,1,M); % 1 x M initial vector 
            // Z = [Y; X];
            // S  = cumprod(Z);  % paths
            // Initialize forward rate with its at time t = 0 value; 
            var S = Matrix<double>.Build.Dense(nTimes, nPaths);
            
            S.SetRow(0, Vector<double>.Build.Dense(nPaths, S0));
            for (var i = 0; i < nTimes; i++)
            {
                var prod = Vector<double>.Build.Dense(nPaths, 1);
                for (var j = 0; j < i; j++)
                {
                    prod = prod.PointwiseMultiply(dX.Row(j));
                }
                S.SetRow(i, prod);        
            }
            Console.WriteLine(S);

            // put option payoff
            // P = max(K-S(end,:),0)
            var P = S.Row(S.RowCount - 1)
                     .Multiply(-1)
                     .Add(K)
                     .Map(p => Math.Max(p, 0));
            
            // loop backwards in time
            for (var n = S.RowCount - 1; n > 0; n--)
            {
                P = P.Multiply(Math.Exp(-r*dt));        
                // use all paths
                // var ind  = 1: M; 
                // X = S.Row(n).
                // X    = S[n,ind]';
                //  Y    = P(ind)';                                 % important payoff
            }

            //for n = N : -1 : 2
            //  P    = exp(-r*dt)*P;                            % discount 
            //  ind  = 1:M;                                     % use all paths
            //  X    = S(n,ind)';
            //  Y    = P(ind)';                                 % important payoffs 
            //  A    = [ ones(size(X)) (1-X) 1/2*(2-4*X-X.^2)]; % 3 basis functions (nice numerical properties of the corresponding matrix)
            //  beta = A\Y;                                     % linear regression (check backslash operator in MATLAB)
            //  betas(n,:) = beta';
            //  C    = A*beta;                                  % continuation value

            //  ind  = ind(find(K-X > C));                       % identifies the indices of the optimal strategy boundary
            //  P(ind) = K - S(n,ind);
            //end

            //P   = exp(-r*dt)*P;
            //val = sum(P)/M;
            //sd  = sqrt((sum(P.^2)/M - val^2)/M);

            //fprintf(' Longstaff-Schwartz Monte Carlo method \n');
            //fprintf(' phase 1: option value = %f, std dev = %f \n',val, sd);
        }










        #endregion Methods

    }
}
