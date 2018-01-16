// ReSharper disable once CheckNamespace

using System;
using System.Linq;
using MathNet.Numerics;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Random;
using PricingLibrary.General;
using PricingLibrary.MarketData;

namespace PricingLibrary.MarketModels
{
    /// <summary>
    /// LIBOR market model
    /// F(i) = LIBOR(t, T[i], T[i+1]) is the forward LIBOR rate accruing from T[i] to T[i+1],
    /// Check http://numerics.mathdotnet.com/ for numerical functions used in the model 
    /// </summary>
    public class LiborMarketModel
    {
        #region Properties

        /// <summary>
        /// Times relevant of the chosen derivative
        /// </summary>
        public Vector<double> T { private set; get; }

        /// <summary>
        /// // day count fraction of the accrual period the LIBOR
        /// </summary>
        public Vector<double> Alpha { private set; get; }
        
        /// <summary>
        /// 
        /// </summary>
        public Vector<double> Delta { private set; get; }
        
        /// <summary>
        /// Number of Libor rates 
        /// </summary>
        public int NumRates { private set; get; }

        /// <summary>
        /// Instantaneous volatility funcitons 
        /// </summary>
        public Func<double, double>[] Volatility { private set; get; }

        /// <summary>
        /// Instantaneous volatilities between Brownian motions
        /// </summary>
        private Matrix<double> Correlation { set; get; }
        
        /// <summary>
        /// Settle date  
        /// </summary>
        public  DateTime SettleDate { private set; get; }
   
        /// <summary>
        /// Initial forward LIBOR rate L(0, T[i-1],T[i]) ;
        /// </summary>
        public Vector<double> InitialForwardRate { private set; get; }

        /// <summary>
        /// Intial numeraire adjusted discount factor B(0,T[i])
        /// where B(t,T) is the price of a Bond with maturity T and time t
        /// </summary>
        public Vector<double> InitialDiscountFactor { private set; get; }

        /// <summary>
        /// Random generator 
        /// </summary>
        public MersenneTwister RandomGenerator { private set; get; }

        /// <summary>
        /// Index of the numeraire 
        /// </summary>
        public int NumeraireIndex { private set; get; }

        /// <summary>
        /// Underlying yield curve 
        /// </summary>
        public YieldCurve YieldCurve { private set; get; }

        #endregion Properties

        #region Constructor 

        /// <summary>
        /// Standard constructor for Libor market model 
        /// </summary>
        /// <param name="zeroCurve"></param>
        /// <param name="correlationBeta"></param>
        /// <param name="volatilityParamaters"></param>
        /// <param name="tenorStructure"></param>
        /// <param name="randomSeed"></param>
        public LiborMarketModel(YieldCurve zeroCurve, double correlationBeta, double[] volatilityParamaters, DateTime[] tenorStructure, int randomSeed = 42) 
            : this(zeroCurve, null, null, tenorStructure, randomSeed)
        {
            #region Correlation

            Func<int, int, double> correlationFunction = (i, j) => Math.Exp(-correlationBeta * (Math.Abs(T[i] - T[j])));

            // Build corellation matrix dW(t) * dW(t)^T = correlation * dt
            Correlation = Matrix<double>.Build
                                        .Dense(NumRates, NumRates)
                                        .MapIndexed((i, j, p) => correlationFunction(i, j));
            #endregion Correlation

            #region Volatility

            // Check volatility parameters 
            if (volatilityParamaters == null || volatilityParamaters.Count() != 4)
            {
                throw new ArgumentException("There are four volatility parameters expected");
            }

            // Use the volatility term structure populaized by Riccardo Rebonato
            // sigma_i(t) =  (a * tau + b) * Math.Exp(-c * tau) + d; with tau = T[i] - t
            // t is assumed to be in units of [years]
            Func<double, double> volatilityFunc = t =>
            {
                var a = volatilityParamaters[0];
                var b = volatilityParamaters[1];
                var c = volatilityParamaters[2];
                var d = volatilityParamaters[3];

                return (a * t + b) * Math.Exp(-c * t) + d;
            };

         
            //  Vector of volatility functions, each corresponds to forward 
            // LIBOR rate of the same index 
            Volatility = new Func<double, double>[NumRates].Select(p => volatilityFunc)
                                                           .ToArray();
            #endregion Volatility 
        }

        /// <summary>
        /// More general Constructor
        /// </summary>
        /// <param name="yieldCurve"></param>
        /// <param name="correlation"></param>
        /// <param name="volatility"></param>
        /// <param name="tenorStructure"></param>
        /// <param name="randomSeed"></param>
        public LiborMarketModel(YieldCurve yieldCurve, 
                                Matrix<double> correlation , 
                                Func<double, double>[] volatility, 
                                DateTime[] tenorStructure, 
                                int randomSeed = 42 )
        {

            YieldCurve = yieldCurve;

            RandomGenerator = new MersenneTwister(randomSeed);
            
            Correlation = correlation;
            Volatility = volatility;
            SettleDate = yieldCurve.SettleDate;
            // Use continuous time in units of years for tenor structure
            T = Vector<double>.Build
                              .DenseOfEnumerable(Utilities.ConvertToDayCountFraction(SettleDate, tenorStructure));

            // Initial discount factors: T[i] :-> B(0, T[i]) i =0, ..., n = length(T)-1
            InitialDiscountFactor = yieldCurve.GetDiscountFactors(tenorStructure);

            #region Tenors

            // tenors (day-count fraction) in units of [years]
            // Alpha[i] = T[i]-T[i+1], i = 0,...,n - 1
            Alpha = Vector<double>.Build.DenseOfEnumerable(T.Skip(1)) - Vector<double>.Build.DenseOfEnumerable(T.Take(T.Count - 1));
            
            // Choose terminal measure, can be an value between 0,..., n = length(T)-1
            NumeraireIndex = T.Count - 1;
           
            #endregion Tenors

            #region Initial Forward rate
            
            NumRates = T.Count - 1;

            // B[i]
            var zeroBond0 = Vector<double>.Build
                                          .DenseOfEnumerable(InitialDiscountFactor.Take(NumRates));
            
            // B[i+1]            
            var zeroBond1 = Vector<double>.Build
                                          .DenseOfEnumerable(InitialDiscountFactor.Skip(1));

            // Initial forward LIBOR rate 
            // L[i] = L(0, T[i],T[i+1]) = 1/alpha[i] * (B[i]/B[i+1] - 1); i = 0,..., n - 1
            // L[i](t) : 0 <= t <= T[i]
            InitialForwardRate = (zeroBond0.PointwiseDivide(zeroBond1).Subtract(1)).PointwiseDivide(Alpha);
            
            #endregion Initial Forward rate

            // Delta for smile 
            Delta = Vector<double>.Build.Dense(NumRates);

        }

        #endregion Constructor 

        #region Methods
        
        /// <summary>
        /// Simulate forward rates and corresponding zero rates on the
        /// interval [settleDate, tmax]
        /// step in [months]
        /// </summary>
        public LiborMarketModelPath SimulateForwardRates(DateTime maxDate, int step=3)
        {
            #region Numeraire 

            // Numeraire at time 0: N(0) = B(0, T[k])
            var intialNumeraire = InitialDiscountFactor[NumeraireIndex];
            
            #endregion Numeraire

            #region Time discertization

            // Initialize rates 
            //t = linspace(0, length(model.T(end-1)), 100);
            
            // Use 3 Month time step
            var dates = Utilities.TenorStructure(SettleDate, maxDate, step);
 
            // Convert to continous time in units of [years] 
            var t = Utilities.ConvertToDayCountFraction(SettleDate, dates);

            #endregion Time discertization

            #region Forward rate

            // Initialize forward rate with its at time t = 0 value; 
            var forwaredRates = Matrix<double>.Build.Dense(NumRates, t.Length);

            forwaredRates.SetColumn(0, InitialForwardRate);

            // Propagate forward rates 
            for (var j = 0; j < t.Length - 1; j++)
            {
                var dt = t[j + 1] - t[j];

                var forwardrate = PropagteForwardRate(t[j], forwaredRates.Column(j), NumeraireIndex, dt);
  
                forwaredRates.SetColumn(j + 1, forwardrate);
            }

            #endregion Forward rate

            #region Numeraire adjusted discount factor

            // Initialize numeraire-adjusted-discount factor with its at time t = 0 value; 
            // D[i] = B[i]/B[k], where B[k] is the numeraire , i = 0, ..., n = length(T)-1
            // D[i]: 0 <= t <= min(T[i], T[k])
            // numeraireAdjustedDiscountFactors
            var discountFactors = Matrix<double>.Build.Dense(T.Count, t.Length);
             
            discountFactors.SetColumn(0, InitialDiscountFactor.Divide(intialNumeraire));

            // Numeraire-adjusted discount factors 
            for (var j = 1; j < t.Length; j++)
            {
                // forward rate for time t[j]
                var forwardRate = forwaredRates.Column(j);
                var discountFactor = Vector<double>.Build.Dense(T.Count);

                for (var i = 0; i < discountFactor.Count; i++)
                {
                    #region Check time 

                    // Discount factor is not defined 
                    if (t[j] > Math.Min(T[i], T[NumeraireIndex]))
                    {
                        discountFactor[i] = double.NaN;
                        continue;
                    }
                    
                    #endregion Check time

                    #region Compute numeraire adjusted discount factors from forward rate
                    if (i == NumeraireIndex)
                    {
                        discountFactor[i] = 1;
                    }
                    else if (i < NumeraireIndex)
                    {
                        var prod = 1.0;
                        for (var k = i; k <= NumeraireIndex - 1; k++)
                        {
                            prod *= 1 + Alpha[k] * forwardRate[k];
                        }
                        discountFactor[i] = prod;
                    }
                    else if (i > NumeraireIndex) 
                    {
                        var prod = 1.0;
                        for (var k = NumeraireIndex; k <= i - 1; k++)
                        {
                            prod *= 1 + Alpha[k] * forwardRate[k];
                        }
                        discountFactor[i] = 1.0 / prod;
                    }
                    #endregion Compute numeraire adjusted discount factors from forward rate
                }
                discountFactors.SetColumn(j, discountFactor);
            }

            #endregion Numeraire adjusted discount factor

            #region Result

            return new LiborMarketModelPath
            {
                SettleDate = SettleDate,
                Dates = dates,
                Tenors = Alpha,
                TenorStructure = Utilities.ConvertFromDayCountFraction(SettleDate, T.ToArray()),
                ForwardRates = forwaredRates,
                NumeraireAdjustedDiscountFactors = discountFactors,
                InitialNumeraire = intialNumeraire
            };

            #endregion Result
        }

        /// <summary>
        /// Compute correlation matrix  
        /// </summary>
        /// <param name="t"></param>
        /// <param name="dt"></param>
        /// <returns></returns>
        private Matrix<double> CorrelationMatrix(double t, double dt)
        {
            var C = Matrix<double>.Build.Dense(NumRates, NumRates);

            for (var i = 0; i < NumRates; i++)
            {
                for (var j = 0; j <= i; j++)
                {
                    var i1 = i;
                    var j1 = j;

                    Func<double, double> func = x => Correlation[i1,j1] * Volatility[i1](T[i1] - x) * Volatility[j1](T[j1] - x);
                    
                    var corr = Integrate.OnClosedInterval(func, t, t + dt);
                    
                    C[i, j] = corr;
                    C[j, i] = corr;
                }
            }
            return C;
        }

        /// <summary>
        /// Propagates the i-the forward LIBOR rate  
        /// </summary>
        /// <param name="i"></param>
        /// <param name="Li"></param>
        /// <param name="totalDrift"></param>
        /// <param name="Cii"></param>
        /// <param name="Mi"></param>
        /// <returns></returns>
        private double PropagteForwardRate(int i, double Li, double totalDrift, double Cii, double Mi)
        {
            return (Li + Delta[i]) * Math.Exp(totalDrift - 0.5 * Cii + Mi) - Delta[i];
        }

        /// <summary>
        /// Take care that F_i(t) is only defined for t le T[i] 
        /// Step in the T_k-forward measure for the SDE
        /// dF[i](t) / (F[i] + delta[i]) = mu[i](t,F(t)) * dt + sigma[i](t) * dW(t)
        /// </summary>
        /// <param name="t"></param>
        /// <param name="forwardRate"></param>
        /// <param name="k"></param>
        /// <param name="dt"></param>
        /// <returns></returns>
        private Vector<double> PropagteForwardRate(double t, Vector<double> forwardRate, int k, double dt)
        {
            #region Correlated random numbers 

            // Integrated instantaneous correlation from t 
            // to t+h. Approximate the Ito-Integral 
            // lz := int_t^t+h { rho[i,j] * sigma[i](s) * sigma[j](s) * dW(s) }
            var C = CorrelationMatrix(t, dt);

            // Normally distributed independent random numbers 
            var r = new double[NumRates];
            Normal.Samples(RandomGenerator, r, 0, 1);
            
            // Correlated Gaussiam random variable with 
            // M ~ N(0, C)
            var normal01 = Vector<double>.Build.DenseOfArray(r);
            var lower = C.Cholesky().Factor;
            
            var M = lower.Multiply(normal01);

            #endregion Correlated random numbers

            #region Predictor step

            var predicatorL = Vector<double>.Build.Dense(NumRates);
            var drift = Vector<double>.Build.Dense(NumRates);
            for (var i = 0; i < NumRates; i++)
            {
                if (t + dt > T[i]) continue;

                drift[i] = Drift(t, forwardRate, i, k);
                
                // Esitimate the total drift (integrated over [t, t + dt]) 
                var totalDrift = drift[i] * dt;
                predicatorL[i] = PropagteForwardRate(i, forwardRate[i], totalDrift, C[i,i], M[i]);
            }

            #endregion Predictor step

            #region Corrector step

            var correctorL = Vector<double>.Build.Dense(NumRates);
            for (var i = 0; i < NumRates; i++)
            {
                if (t + dt > T[i]) continue;
            
                // Esitimate the total drift (integrated over [t,t+dt]) 
                // using trapezoidal rule 
                var totalDrift = (Drift(t + dt, predicatorL, i, k) + drift[i]) * dt / 2;
                correctorL[i] = PropagteForwardRate(i, forwardRate[i], totalDrift, C[i, i], M[i]);
            }

            #endregion Corrector step

            return correctorL;
        }

        /// <summary>
        /// Instantaneous drift of the F_i in  the T_k-forward-measure
        ///-F is a n-dim vector of the LIBOR forward rates
        /// -alpha is a n-dim vector containing the tenors of the LIBOR rates
        /// -k indicates that the T_k forward measure is used
        /// - rho is n x n - Matrix containg the correlation betweeen the BM driving
        /// the LIBOR rates
        /// </summary>
        /// <param name="t"></param>
        /// <param name="forwardRate"></param>
        /// <param name="i"></param>
        /// <param name="k"></param>
        /// <returns></returns>
        private double Drift(double  t, Vector<double> forwardRate, int i, int k)
        {
            if (i + 1 == k)
            {
                return 0;
            }
            
            if(i + 1 > k)
            {
                double sum = 0;
                var si = Volatility[i](T[i] - t);
                for (var j = i + 1; j <= k - 1; j++)
                {
                    // TODO: Check volatility of L[j]
                    var sj = Volatility[j](T[j] - t);

                    sum += Alpha[j] * (forwardRate[j] + Delta[j]) / (1 + Alpha[j] * forwardRate[j]) * Correlation[i, j] * si * sj;
                }
                return -sum;
            }
            else
            {
                double sum = 0;
                var si = Volatility[i](T[i] - t);

                for (var j = k; j <= i; j++)
                {
                    // L[t,T[i], T[i+1]] => volatility sigma(T[i]-t)
                    var sj = Volatility[j](T[j] - t);

                    sum += Alpha[j] * (forwardRate[j] + Delta[j]) / (1 + Alpha[j] * forwardRate[j]) * Correlation[i, j] * si * sj;
                }
                return sum;
            }
        }  

        #endregion Methods
    }
    
}