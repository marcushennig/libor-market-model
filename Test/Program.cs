using System;
using System.IO;
using System.Linq;
using PricingLibrary.FinancialInstruments;
using PricingLibrary.General;
using PricingLibrary.MarketData;
using PricingLibrary.MarketModels;

namespace Test
{
    static class Program
    {
        const string DataFolder = @"C:\Users\d90789\Documents\Oxford MSc in Mathematical Finance\Modules\C# Projects\PricingEngine\Matlab\PricingData";
  
        /// <summary>
        /// Export 
        /// </summary>
        private static void TestLiborMarketModel(string fowardrateFile, string yieldCurveFile)
        {
            #region Swaption

            // 1 EUR 
            var nominal = 1;

            // Consider a ten-year semi-annual European swaption (10y x 5y), exercisable on each 
            // reset date, starting on 1-October-2020 
            var valuationDate = new DateTime(2015, 10, 1);

            // Exercise date of the option = start date of swap
            var exerciseDate = valuationDate.AddYears(2);

            // Maturity of the underlying swap
            var maturityDate = exerciseDate.AddYears(10);

            // Strike (fixed interest rate payed by fixed leg) 
            var strike = 0.045;

            // Tenor of the underlying swap in units of [months]
            // = floating leg pays 6M-Libor 
            var tenor = 6;

            var swaption = new Swaption(valuationDate, exerciseDate, maturityDate, tenor, strike, nominal);

            #endregion Swaption

            #region Initialize LIBOR market model

            #region Initial yield curve

            // Sample points of curve in units of [years]
            var curveTimes = new[] { 1, 2, 3, 4, 5, 7, 10, 20 };

            // annual interest rate in units of [1/year]
            var yields = new[] { 0.01, 0.018, 0.024, 0.029, 0.033, 0.034, 0.035, 0.036 };
            var maturityDates = curveTimes.Select(y => valuationDate.AddYears(y)).ToArray();

            var yieldCurve = new YieldCurve(valuationDate, maturityDates, yields);

            #endregion  Initial yield curve

            #region Model parameters

            // Use the volatility term structure populaized by Riccardo Rebonato
            // sigma_i(t) =  (a * tau + b) * Math.Exp(-c * tau) + d; with tau = T[i] - t
            // t is assumed to be in units of [years]
            var volatiltyParameters = new[] { 0.3, -0.02, 0.7, 0.14 };

            // Instantaneous correlations between the Brownian motions 
            // and hence the forward LIBOR rates
            // rho(i,j) =   exp(-beta * |T[i] - T[j]|)
            const double beta = 0.08;

            #endregion Model parameters

            const int seed = 41;
            var liborMarketModel = new LiborMarketModel(yieldCurve, beta, volatiltyParameters, swaption.TenorStructure, seed);

            #endregion Initialize LIBOR market model

            var paths = liborMarketModel.SimulateForwardRates(swaption.ExerciseDate, 1);

            liborMarketModel.YieldCurve.Export(Path.Combine(DataFolder, yieldCurveFile));
            paths.Export(Path.Combine(DataFolder, fowardrateFile));   
           
        }

        /// <summary>
        /// Test cap pricing function 
        /// </summary>
        private static void TestCapPricing(string resultFile)
        {
            #region Cap

            // 1 EUR 
            var principal = 1;

            // Consider a ten-year semi-annual European swaption (10y x 5y), exercisable on each 
            // reset date, starting on 1-October-2020 
            var valuationDate = new DateTime(2015, 10, 1);

            // Exercise date of the option = start date of swap
            var startDate = valuationDate.AddYears(2);

            // Strike (fixed interest rate payed by fixed leg) 
            var capRate = 0.045;

            // Tenor of the underlying swap in units of [months]
            // = floating leg pays 6M-Libor 
            var tenor = 6;

            const int numberOfCaplets = 10;

            var cap = new Cap(startDate, tenor, numberOfCaplets, capRate, principal);
           
            #endregion Cap

            #region Initialize LIBOR market model

            #region Initial yield curve

            // Sample points of curve in units of [years]
            var curveTimes = new[] { 1, 2, 3, 4, 5, 7, 10, 20 };

            // annual interest rate in units of [1/year]
            var yields = new[] { 0.01, 0.018, 0.024, 0.029, 0.033, 0.034, 0.035, 0.034 };
            var maturityDates = curveTimes.Select(y => valuationDate.AddYears(y)).ToArray();

            var yieldCurve = new YieldCurve(valuationDate, maturityDates, yields);

            #endregion  Initial yield curve

            #region Model parameters

            // Use the volatility term structure populaized by Riccardo Rebonato
            // sigma_i(t) =  (a * tau + b) * Math.Exp(-c * tau) + d; with tau = T[i] - t
            // t is assumed to be in units of [years]
            var volatilityParamaters = new[] { 0.3, -0.02, 0.7, 0.14 };

            // Instantaneous correlations between the Brownian motions 
            // and hence the forward LIBOR rates
            // rho(i,j) =   exp(-beta * |T[i] - T[j]|)
            const double beta = 0.08;

            #endregion Model parameters

            var seed = DateTime.Now.Millisecond;

            var liborMarketModel = new LiborMarketModel(yieldCurve, beta, volatilityParamaters, cap.TenorStructure, seed);

            #endregion Initialize LIBOR market model

            #region Pricing

            #region Cap rates

            const int nCapRates = 100;
            const double capRateMin = 0;
            const double capRateMax = 0.4;
            const double step = (capRateMax - capRateMin) / (nCapRates - 1);

            var capRates = new double[nCapRates];
            for (var j = 0; j < capRates.Length; j++)
            {
                capRates[j] = capRateMin + step * j;
            }

            #endregion Cap rates

            #region Monte Carlo

            const int numberOfMonteCarloPaths = 100;
            var prices = new double[nCapRates];

            for (var i = 0; i < numberOfMonteCarloPaths; i++)
            {
                // Simulate forward LIBOR rate paths using Monto carlo simulation up to the
                // exercise date of the swaption 
                var paths = liborMarketModel.SimulateForwardRates(cap.MaturityDate);

                for (var j = 0; j < capRates.Length; j++)
                {
                    cap.CapRate = capRates[j];
                    prices[j] += cap.DiscountedPayoff(paths) / numberOfMonteCarloPaths;
                }
                // Utilities.DrawTextProgressBar(i + 1, numberOfMonteCarloPaths, "Monte Carlo");
            }

            #endregion Monte Carlo

            #region Black-Caplet formula

            // Use the volatility term structure populaized by Riccardo Rebonato
            // sigma_i(t) =  (a * tau + b) * Math.Exp(-c * tau) + d; with tau = T[i] - t
            // t is assumed to be in units of [years]
            Func<double, double> volatilityFunction = t =>
            {
                var a = volatilityParamaters[0];
                var b = volatilityParamaters[1];
                var c = volatilityParamaters[2];
                var d = volatilityParamaters[3];

                return (a * t + b) * Math.Exp(-c * t) + d;
            };
            var blackPrices = new double[nCapRates];

            for (var i = 0; i < capRates.Length; i++)
            {
                cap.CapRate = capRates[i];
                blackPrices[i] = cap.Value(yieldCurve, 0, volatilityFunction);

                Utilities.DrawTextProgressBar(i + 1, capRates.Length, "Black-Caplet formula");
            }

            #endregion Black-Caplet formula

            #endregion Pricing

            #region Export
            
            try
            {
                // create a writer and open the file
                var file = new StreamWriter(resultFile);
                for (var i = 0; i < capRates.Length; i++)
                {
                    file.WriteLine("{0}\t{1}\t{2}", capRates[i], prices[i], blackPrices[i]);

                    Utilities.DrawTextProgressBar(i + 1, capRates.Length, "Export results");
                }
                file.Close();
            }
            catch (Exception exception)
            {
                Console.WriteLine(exception.Message);
            }

            #endregion Export
        }
        
        /// <summary>
        /// Test caplet pricing function 
        /// </summary>
        private static void TestCapletPricing(string resultFile)
        {
            
            #region Caplet 
            
            // 1 EUR 
            var principal = 1;

            // Consider a ten-year semi-annual European swaption (10y x 5y), exercisable on each 
            // reset date, starting on 1-October-2020 
            var valuationDate = new DateTime(2015, 10, 1);

            // Exercise date of the option = start date of swap
            var fixingDate = valuationDate.AddYears(2);

            // Strike (fixed interest rate payed by fixed leg) 
            var capRate = 0.045;

            // Tenor of the underlying swap in units of [months]
            // = floating leg pays 6M-Libor 
            var tenor = 6;

            var caplet = new Caplet(fixingDate, tenor, capRate, principal);

            #endregion Caplet

            #region Initialize LIBOR market model

            #region Initial yield curve

            // Sample points of curve in units of [years]
            var curveTimes = new[] { 1, 2, 3, 4, 5, 7, 10, 20 };

            // annual interest rate in units of [1/year]
            var yields = new[] { 0.01, 0.018, 0.024, 0.029, 0.033, 0.034, 0.035, 0.034 };
            var maturityDates = curveTimes.Select(y => valuationDate.AddYears(y)).ToArray();

            var yieldCurve = new YieldCurve(valuationDate, maturityDates, yields);

            #endregion  Initial yield curve

            #region Model parameters

            // Use the volatility term structure populaized by Riccardo Rebonato
            // sigma_i(t) =  (a * tau + b) * Math.Exp(-c * tau) + d; with tau = T[i] - t
            // t is assumed to be in units of [years]
            var volatilityParamaters = new[] { 0.3, -0.02, 0.7, 0.14 };

            // Instantaneous correlations between the Brownian motions 
            // and hence the forward LIBOR rates
            // rho(i,j) =   exp(-beta * |T[i] - T[j]|)
            const double beta = 0.08;

            #endregion Model parameters

            var seed = DateTime.Now.Millisecond;

            var lmm = new LiborMarketModel(yieldCurve, beta, volatilityParamaters, caplet.TenorStructure, seed);

            #endregion Initialize LIBOR market model

            #region Pricing

            #region Cap rates 

            const int nCapRates = 100;
            const double capRateMin = 0;
            const double capRateMax = 0.4;
            const double step = (capRateMax - capRateMin)/(nCapRates - 1);
            
            var capRates = new double[nCapRates];
            for (var j = 0; j < capRates.Length; j++)
            {
                capRates[j] = capRateMin + step * j;
            }

            #endregion Cap rates

            #region Monte Carlo

            const int numberOfMonteCarloPaths = 100;
            const int numberOfMonteCarloRuns = 100;


            var prices = new double[numberOfMonteCarloRuns, nCapRates];
         
            for (var k = 0; k < numberOfMonteCarloRuns; k++)
            {
                Utilities.DrawTextProgressBar(k + 1, numberOfMonteCarloRuns, "Monte Carlo run");
                for (var i = 0; i < numberOfMonteCarloPaths; i++)
                {
                    // Simulate forward LIBOR rate paths using Monto carlo simulation up to the
                    // exercise date of the swaption 
                    var paths = lmm.SimulateForwardRates(caplet.ExpiryDate);
                    
                    for (var j = 0; j < capRates.Length; j++)
                    {
                        caplet.CapRate = capRates[j];

                        var price = caplet.DiscountedPayoff(paths);

                        prices[k, j] += price/numberOfMonteCarloPaths;
                    }
                }
            }

            #endregion Monte Carlo

            #region Black-Caplet formula

            // Use the volatility term structure populaized by Riccardo Rebonato
            // sigma_i(t) =  (a * tau + b) * Math.Exp(-c * tau) + d; with tau = T[i] - t
            // t is assumed to be in units of [years]
            Func<double, double> volatilityFunction = t =>
            {
                var a = volatilityParamaters[0];
                var b = volatilityParamaters[1];
                var c = volatilityParamaters[2];
                var d = volatilityParamaters[3];

                return (a * t + b) * Math.Exp(-c * t) + d;
            };
            var blackPrices = new double[nCapRates];

            for (var i = 0; i < capRates.Length; i++)
            {
                caplet.CapRate = capRates[i];
                blackPrices[i] = caplet.Value(yieldCurve, 0, volatilityFunction);
                
                Utilities.DrawTextProgressBar(i + 1, capRates.Length, "Black-Caplet formula");
            }

            #endregion Black-Caplet formula

            #endregion Pricing

            #region Export

            try
            {
                // create a writer and open the file
                var file = new StreamWriter(resultFile);
                for (var j = 0; j < capRates.Length; j++)
                {
                    file.Write("{0}\t", capRates[j]);

                    for (var k = 0; k < numberOfMonteCarloRuns; k++)
                    {
                        file.Write("{0}\t", prices[k, j]);
                    }
                    file.Write("{0}\n", blackPrices[j]);
             

                    Utilities.DrawTextProgressBar(j + 1, capRates.Length, "Export results");
                }
                file.Close();
            }
            catch (Exception exception)
            {
                Console.WriteLine(exception.Message);
            }
            
            #endregion Export
        }

        /// <summary>
        /// TODO: Get this code into MATLAB using the .NET interface
        /// Test swaption pricing 
        /// </summary>
        private static void TestSwaptionPricing(string resultFile)
        {
            #region Swaption

            // 1 EUR 
            var nominal = 1;

            // Consider a ten-year semi-annual European swaption (10y x 5y), exercisable on each 
            // reset date, starting on 1-October-2020 
            var valuationDate = new DateTime(2015, 10, 1);

            // Exercise date of the option = start date of swap
            var exerciseDate = valuationDate.AddYears(2);

            // Maturity of the underlying swap
            var maturityDate = exerciseDate.AddYears(10);

            // Strike (fixed interest rate payed by fixed leg) 
            var strike = 0.045;

            // Tenor of the underlying swap in units of [months]
            // = floating leg pays 6M-Libor 
            var tenor = 6;

            var swaption = new Swaption(valuationDate, exerciseDate, maturityDate, tenor, strike, nominal);

            #endregion Swaption

            #region Initialize LIBOR market model

            #region Initial yield curve

            // Sample points of curve in units of [years]
            var curveTimes = new[] { 1, 2, 3, 4, 5, 7, 10, 20 }; 

            // annual interest rate in units of [1/year]
            var yields = new[] { 0.01, 0.018, 0.024, 0.029, 0.033, 0.034, 0.035, 0.034 };
            var maturityDates = curveTimes.Select(y => valuationDate.AddYears(y)).ToArray();

            var yieldCurve = new YieldCurve(valuationDate, maturityDates, yields);
            
            #endregion  Initial yield curve

            #region Model parameters

            // Use the volatility term structure populaized by Riccardo Rebonato
            // sigma_i(t) =  (a * tau + b) * Math.Exp(-c * tau) + d; with tau = T[i] - t
            // t is assumed to be in units of [years]
            var volatiltyParameters = new[] { 0.3, -0.02, 0.7, 0.14 };

            // Instantaneous correlations between the Brownian motions 
            // and hence the forward LIBOR rates
            // rho(i,j) =   exp(-beta * |T[i] - T[j]|)
            const double beta = 0.08;

            #endregion Model parameters

            const int seed = 41;
            var liborMarketModel = new LiborMarketModel(yieldCurve, beta, volatiltyParameters, swaption.TenorStructure, seed);
           
            #endregion Initialize LIBOR market model
            
            #region Strikes

            const int nStrikes = 100;
            const double strikeMin = 0;
            const double strikeMax = 0.1;
            const double step = (strikeMax - strikeMin) / (nStrikes - 1);

            var strikes = new double[nStrikes];
            for (var j = 0; j < strikes.Length; j++)
            {
                strikes[j] = strikeMin + step * j;
            }

            #endregion Strikes

            #region Pricing

            const int numberOfMonteCarloPaths = 100;
            const int numberOfMonteCarloRuns = 100;


            var prices = new double[numberOfMonteCarloRuns, nStrikes];

            for (var k = 0; k < numberOfMonteCarloRuns; k++)
            {
                Utilities.DrawTextProgressBar(k + 1, numberOfMonteCarloRuns, "Monte Carlo run");

                for (var i = 0; i < numberOfMonteCarloPaths; i++)
                {
                    // Simulate forward LIBOR rate paths using Monto carlo simulation up to the
                    // exercise date of the swaption 
                    var paths = liborMarketModel.SimulateForwardRates(swaption.ExerciseDate);

                 
                    for (var j = 0; j < strikes.Length; j++)
                    {
                        swaption.Strike = strikes[j];

                        prices[k, j] += swaption.DiscountedPayoff(paths) / numberOfMonteCarloPaths;
                    }
                }
            }
            #endregion Pricing
            
            #region Export

            try
            {
                // create a writer and open the file
                var file = new StreamWriter(resultFile);
                for (var j = 0; j < strikes.Length; j++)
                {
                    file.Write("{0}\t", strikes[j]);

                    for (var k = 0; k < numberOfMonteCarloRuns; k++)
                    {
                        file.Write("{0}", prices[k, j]);
                        file.Write(k < numberOfMonteCarloRuns - 1 ? "\t" : "\n");
                    }


                    Utilities.DrawTextProgressBar(j + 1, strikes.Length, "Export results");
                }
                file.Close();
            }
            catch (Exception exception)
            {
                Console.WriteLine(exception.Message);
            }

            #endregion Export
        }

        /// <summary>
        /// Main entry point
        /// </summary>
        /// <param name="args"></param>
        static void Main(string[] args)
        {
            AmericanOption.TestFunctions();

            //TestSwaptionPricing(Path.Combine(DataFolder, "swaption_pricing.txt"));
            //TestLiborMarketModel(Path.Combine(DataFolder, "lmm_forward_rates.txt"), Path.Combine(DataFolder, "yield_curve.txt"));
            //TestCapletPricing(Path.Combine(DataFolder, "caplet_pricing.txt"));
            //TestCapPricing(Path.Combine(DataFolder, "cap_pricing.txt"));

            //Console.WriteLine("\n\nPress any key");
            //Console.ReadKey();
        }
    }
}
