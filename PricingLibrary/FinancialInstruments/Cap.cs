using System;
using System.Collections.Generic;
using System.Linq;
using PricingLibrary.General;
using PricingLibrary.MarketData;
using PricingLibrary.MarketModels;

namespace PricingLibrary.FinancialInstruments
{
    /// <summary>
    /// Class to decribe an interest rate cap
    /// </summary>
    public class Cap
    {
        #region Fields

        private double _capRate;
        private double _principal;

        #endregion Fields 

        #region Properties

        /// <summary>
        /// Underling caplets 
        /// </summary>
        public List<Caplet> Caplets { private set; get; }

        /// <summary>
        /// Strike or cap rate in units of [% / year] 
        /// </summary>
        public double CapRate
        {
            set
            {
               _capRate = value;
               
                if (Caplets == null) return;

                foreach (var caplet in Caplets)
                {
                    caplet.CapRate = value;
                }
            }
            get { return _capRate; }
        }

        /// <summary>
        /// Notional amount (nominal)
        /// </summary>
        public double Principal
        {
            set
            {
                _principal = value;
                if (Caplets == null) return;
                foreach (var caplet in Caplets)
                {
                    caplet.CapRate = value;
                }
            }
            get { return _principal;  }
        }

        /// <summary>
        /// Currency of the principal
        /// </summary>
        public Currency PrincipalCurreny { private set; get; }

        /// <summary>
        /// first Fixing date 
        /// </summary>
        public DateTime FirstFixingDate { private set; get; }

        /// <summary>
        /// Maturity date of the cap 
        /// </summary>
        public DateTime MaturityDate 
        {
            get
            {
                return TenorStructure[Caplets.Count];
            }
        }

        /// <summary>
        /// Number of caplets
        /// </summary>
        public int NumberOfCaplets 
        {
            get { return Caplets.Count;  }
        }
        
        /// <summary>
        /// Tenor structure of the caplet 
        /// </summary>
        public DateTime[] TenorStructure { private set; get; }
        
        /// <summary>
        /// Tenor  
        /// </summary>
        public int Tenor { private set; get; }

        #endregion Properties

        #region Constructor 

        /// <summary>
        /// Standard constructor of interest cap 
        /// </summary>
        /// <param name="fixingDate"></param>
        /// <param name="tenor"></param>
        /// <param name="numberOfCaplets"></param>
        /// <param name="capRate"></param>
        /// <param name="principal"></param>
        public Cap(DateTime fixingDate, int tenor, int numberOfCaplets, double capRate, double principal)
        {
            Caplets = new List<Caplet>();
            TenorStructure = new DateTime[numberOfCaplets+1];

            var T = fixingDate;
            for (var  i = 0; i < numberOfCaplets; i++)
            {
                TenorStructure[i] = T;

                Caplets.Add(new Caplet(T, tenor, capRate, principal));
                T = T.AddMonths(tenor);
            }
            TenorStructure[numberOfCaplets] = T;
            
        
            FirstFixingDate = fixingDate;

            CapRate = capRate;
            Principal = principal;
            PrincipalCurreny = Currency.EUR;
            Tenor = tenor;
        }

        #endregion Constructor

        #region Methods

        /// <summary>
        /// Determine the value of the caplet using the Black-Caplet-formula 
        /// </summary>
        /// <returns></returns>
        public double Value(YieldCurve yieldCurve, double delta, Func<double, double> volatilityFunction)
        {
           return Caplets.Sum(caplet => caplet.Value(yieldCurve, delta, volatilityFunction));
        }

        /// <summary>
        /// Discounted payoff of  caplet
        /// Payoff P(T+alpha) = max(N * alpha * L(T,T,T+alpa) - K, 0)
        /// </summary>
        /// <param name="paths"></param>
        /// <returns>Payoff at exercise date</returns>
        public double DiscountedPayoff(LiborMarketModelPath paths)
        {
            return Caplets.Sum(caplet => caplet.DiscountedPayoff(paths));
        }

        #endregion Methods
    }
}