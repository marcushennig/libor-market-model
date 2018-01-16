using System;
using System.Security.Policy;
using MathNet.Numerics;
using PricingLibrary.MarketModels;
using PricingLibrary.General;
using MathNet.Numerics.Distributions;
using PricingLibrary.MarketData;

namespace PricingLibrary.FinancialInstruments
{
    /// <summary>
    /// A caplet is an assurance against high interest rate 
    /// </summary>
    public class Caplet
    {
        #region Properties

        /// <summary>
        /// cap rate in units of [1 / year] 
        /// </summary>
        public double CapRate
        {
            set; 
            get;
        }
        
        /// <summary>
        /// Notional amount (nominal)
        /// </summary>
        public double Principal { set; get; }
        
        /// <summary>
        /// Currency of the notioal amount 
        /// </summary>
        public Currency NominalCurreny { set; get; }
        
        /// <summary>
        /// Expiry data of caplet
        /// </summary>
        public DateTime FixingDate { set; get; }

        /// <summary>
        /// Expiry date of caplet
        /// </summary>
        public DateTime  ExpiryDate { set; get; }
        
        /// <summary>
        ///  Day count fraction corresponding to the 
        /// period to which L applies
        /// </summary>
        public double Tenor { set; get; }

        /// <summary>
        /// Tenor structure of the caplet 
        /// </summary>
        public DateTime[] TenorStructure { private set; get; }

        #endregion Properties

        #region Constructor

        /// <summary>
        /// Standard constructror 
        /// </summary>
        /// <param name="fixingDate"></param>
        /// <param name="tenor"></param>
        /// <param name="capRate"></param>
        /// <param name="principal"></param>
        public Caplet(DateTime fixingDate, int tenor, double capRate, double principal)
        {
            // T
            FixingDate = fixingDate;
            
            // T + alpha
            ExpiryDate = FixingDate.AddMonths(tenor);

            CapRate = capRate;
            Principal = principal;
            NominalCurreny = Currency.EUR;
            Tenor = tenor;

            TenorStructure = new [] {FixingDate, ExpiryDate};
        }

        #endregion Constructor 

        #region Methods

        /// <summary>
        /// Determine the value of the caplet using the Black-Caplet-formula 
        /// </summary>
        /// <returns></returns>
        public double Value(YieldCurve yieldCurve, double delta, Func<double, double> volatilityFunc)
        {
            // valuation date 
            var valuationDate = yieldCurve.SettleDate;

            var T = Utilities.ConvertToDayCountFraction(valuationDate, FixingDate);

            Func<double, double> func = t => volatilityFunc(T - t) * volatilityFunc(T - t);

            var sigma = Math.Sqrt(1 / T * Integrate.OnClosedInterval(func, 0, T));
            
            var tenor = Utilities.ConvertToDayCountFraction(FixingDate, ExpiryDate);
            var bond = yieldCurve.GetDiscountFactors(new [] {FixingDate, ExpiryDate});
            
            var kappa = CapRate/(Tenor*Principal) + delta;

            var spotForwardRate = 1 / tenor * (bond[0] / bond[1] - 1);

            var logLK = Math.Log((spotForwardRate + delta) / (kappa));
            var gamma = T * sigma * sigma;
            var dPlus = (logLK + 0.5 * gamma) / Math.Sqrt(gamma);
            var dMinus = (logLK - 0.5 * gamma) / Math.Sqrt(gamma);

            return bond[1] * Principal * Tenor * (spotForwardRate * Normal.CDF(0, 1, dPlus) - kappa * Normal.CDF(0, 1, dMinus));
        }

        /// <summary>
        /// Discounted payoff of  caplet
        /// Payoff P(T+alpha) = max(N * alpha * L(T,T,T+alpa) - K, 0)
        /// </summary>
        /// <param name="paths"></param>
        /// <returns>Payoff at exercise date</returns>
        public double DiscountedPayoff(LiborMarketModelPath paths)
        {
            // Choose L(T,T, T + alpha)
            var j = Array.IndexOf(paths.Dates, FixingDate);
            var i = Array.IndexOf(paths.TenorStructure, FixingDate);

            var forwardRate = paths.ForwardRates[i, j];
            
            // Choose D(T+alpha) = B(T+alpha,T+alpha)/B(T+alpha,Tk)
            var jj= Array.IndexOf(paths.Dates, ExpiryDate);
            var ii = Array.IndexOf(paths.TenorStructure, ExpiryDate);

            var numeraireAdjustedDiscountFactor = paths.NumeraireAdjustedDiscountFactors[ii, jj];

            // payoff of the swaption
            return paths.InitialNumeraire * numeraireAdjustedDiscountFactor * Math.Max(Tenor * Principal * forwardRate - CapRate, 0);
        }


        #endregion Methods
    }
}
