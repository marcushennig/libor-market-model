using System;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;
using PricingLibrary.MarketModels;
using PricingLibrary.General;

namespace PricingLibrary.FinancialInstruments
{
/// <summary>
/// European option on an interest rate swap
/// </summary>
public class Swaption
{
    #region Properties

    /// <summary>
    /// Strike [rate of the fixed leg of the underlying swap]
    /// of the swaption in units of [1/year]
    /// </summary>
    public double Strike { set; get; }

    /// <summary>
    /// Notional amount (nominal)
    /// </summary>
    public double Nominal { set; get; }

    /// <summary>
    /// Currency of the notioal amount 
    /// </summary>
    public Currency NominalCurreny { private set; get; }

    /// <summary>
    /// Tenor of the underlying swap in units of [months] 
    /// </summary>
    public int Tenor { private set; get; }

    /// <summary>
    /// Tenor structure of swaption
    /// </summary>
    public DateTime[] TenorStructure { private set; get; }

    /// <summary>
    /// Settle date
    /// </summary>
    public DateTime SettleDate { private set; get; }

    /// <summary>
    /// Execise date of the option
    /// </summary>
    public DateTime ExerciseDate { private set; get; }
        
    /// <summary>
    /// Maturity date of the swap 
    /// </summary>
    public DateTime MaturityDate { private set; get; }

    #endregion Properties

    #region Constructor 

    /// <summary>
    /// Constructs a swaption
    /// </summary>
    public Swaption(DateTime settleDate, DateTime exerciseDate, DateTime maturityDate, int tenor, double strike, double nominal)
    {
        SettleDate = settleDate;
        ExerciseDate = exerciseDate;
        MaturityDate = maturityDate;

        Strike = strike;
        Nominal = nominal;
        NominalCurreny = Currency.EUR;
        Tenor = tenor;

        // compute tenor structure of the swaption
        TenorStructure = Utilities.TenorStructure(ExerciseDate, MaturityDate, Tenor);
    }

    #endregion Constructor

    #region Methods

    /// <summary>
    /// Present value 
    /// </summary>
    /// <param name="numeraire0"> B(0,T[N]) </param>
    /// <param name="tenor">alpha[i]</param>
    /// <param name="forwardRate">F[i](T[0])</param>
    /// <param name="numeraireAdjustedDiscountFactor">B(T[0],T[i+1])/B(T[0],T[N)</param>
    /// <param name="strike"></param>
    /// <param name="nominal"></param>
    /// <returns></returns>
    private double PresentValueUnderlyingSwap(double numeraire0, 
                                                Vector<double> tenor,
                                                Vector<double> forwardRate,
                                                Vector<double> numeraireAdjustedDiscountFactor, 
                                                double strike, double nominal)
    {
        return numeraire0 * forwardRate.Select((L, i) => numeraireAdjustedDiscountFactor[i] * nominal * tenor[i] * (L - strike))
                                        .Sum(); ;
    }

    /// <summary>
    /// Determine present value of swaption from fiven paths of 
    /// forward LIBOR  rates and  numeraire adjusted discount factors
    /// </summary>
    /// <param name="paths"></param>
    /// <returns>Payoff at exercise date</returns>
    public double DiscountedPayoff(LiborMarketModelPath paths)
    {
        // European swaption, hence there is no path dependence
        // therefore only we only need state of model at the 
        // exercise date
        var j = Array.IndexOf(paths.Dates, ExerciseDate);
        var forwardRates = paths.ForwardRates.Column(j);            
        var numeraireAdjustedDiscountFactors = paths.NumeraireAdjustedDiscountFactors.Column(j);
            
        // present value of the underling swap at exercise date
        var presentValue = PresentValueUnderlyingSwap(paths.InitialNumeraire, paths.Tenors, 
                                                        forwardRates, numeraireAdjustedDiscountFactors,
                                                        Strike, Nominal);

        // payoff of the swaption
        return Math.Max(presentValue, 0);
    }

    #endregion
}
}