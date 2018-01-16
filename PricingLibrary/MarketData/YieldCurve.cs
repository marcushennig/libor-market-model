using System;
using System.IO;
using MathNet.Numerics.Interpolation;
using MathNet.Numerics.LinearAlgebra;
using PricingLibrary.General;

namespace PricingLibrary.MarketData
{
    /// <summary>
    /// Yield curve 
    /// </summary>
    public class YieldCurve
    {
        #region Properties
        
        /// <summary>
        /// Settle date 
        /// </summary>
        public DateTime SettleDate { private set; get; }

        /// <summary>
        /// Curve dates, vectors of DateTimes are not supported by Math.Net yet 
        /// </summary>
        public DateTime[] MaturityDates { private set; get;  }

        /// <summary>
        /// Rates 
        /// </summary>
        public Vector<double> Yield { private set; get; }
        
        /// <summary>
        /// B(0,T) devided by its face value (one unit of money)
        /// </summary>
        public Vector<double> DiscountFactor { private set; get; }
       

        /// <summary>
        /// Maturity in contineous time in units of years
        /// </summary>
        private Vector<double> T { set; get; }

        #endregion Properties

        #region Constructor 

        /// <summary>
        /// Constructor  
        /// </summary>
        public YieldCurve(DateTime settleDate,  DateTime[] maturityDates, double[] yield)
        {
            // Settle date 
            SettleDate = settleDate.Date;
            
            // Maturity dates
            MaturityDates = maturityDates;

            // Convert maturities to contineous time [year]
            T = Vector<double>.Build.DenseOfEnumerable(Utilities.ConvertToDayCountFraction(SettleDate, MaturityDates));

            // yield curve T -> Y(0,T) = -1/T*log[B(0,T)] in units of [1/year]
            Yield = Vector<double>.Build.DenseOfArray(yield);

            // Discount curve T->B(0,T) = exp(-T*B(0,T)) [1]
            DiscountFactor = (-T.PointwiseMultiply(Yield)).PointwiseExp();
        }

        #endregion

        #region Methods

        /// <summary>
        /// Export to file 
        /// </summary>
        /// <param name="filePath"></param>
        public void Export(string filePath)
        {

            try
            {
                var file = new StreamWriter(filePath);
                
                // save tenor structure in first line 
                for (var i = 0; i < T.Count; i++)
                {
                    file.WriteLine("{0}\t{1}", T[i], Yield[i]);
                }
                file.Close();
            }
            catch (Exception exception)
            {
                Console.WriteLine(exception.Message);
            }
        }

        /// <summary>
        /// Discount factors
        /// </summary>
        /// <param name="dates"></param>
        /// <returns></returns>
        public Vector<double> GetDiscountFactors(DateTime[] dates)
        {
            // Convert dates to contineous time in untis of years 
            var t = Vector<double>.Build.DenseOfArray(Utilities.ConvertToDayCountFraction(SettleDate, dates));

            var interpolation = new BulirschStoerRationalInterpolation(T.ToArray(), DiscountFactor.ToArray());

            var result = Vector<double>.Build.Dense(t.Count);
                
            t.Map(p => interpolation.Interpolate(p), result);

            return result;
        }

        /// <summary>
        /// Discount factor 
        /// </summary>
        /// <param name="date"></param>
        /// <returns></returns>
        public double GetDiscountFactors(DateTime date)
        {
            // Convert dates to contineous time in untis of years 
            var t = Utilities.ConvertToDayCountFraction(SettleDate, date);

            var interpolation = new BulirschStoerRationalInterpolation(T.ToArray(), DiscountFactor.ToArray());

            return interpolation.Interpolate(t);
        }

        #endregion Methods
    }
}
