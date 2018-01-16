using System;
using System.IO;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;
using PricingLibrary.General;

namespace PricingLibrary.MarketModels
{
/// <summary>
/// Result object 
/// </summary>
public class LiborMarketModelPath
{
    #region Properties

    public DateTime SettleDate { set; get; }
    public DateTime[] Dates { set; get; }
    public DateTime[] TenorStructure { set; get; }
    public Vector<double> Tenors { set; get; }
    public Matrix<double> ForwardRates { set; get; }
    public Matrix<double> NumeraireAdjustedDiscountFactors { set; get; }
    public double InitialNumeraire { set; get; }

    #endregion Properties

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

            var forwardRates = ForwardRates;
            var t = Utilities.ConvertToDayCountFraction(Dates.First(), Dates);
            var T = Utilities.ConvertToDayCountFraction(SettleDate, TenorStructure);

            // save tenor structure in first line 
            for (var i = 0; i < T.Length; i++)
            {
                file.Write("{0}", T[i]);
                file.Write(i < T.Length - 1 ? "\t" : "\n");
            }

            // loop over all time slices 
            for (var j = 0; j < t.Length; j++)
            {
                file.Write("{0}\t", t[j]);
                for (var i = 0; i < forwardRates.RowCount; i++)
                {
                    file.Write("{0}", forwardRates[i, j]);
                    file.Write(i < forwardRates.RowCount - 1 ? "\t" : "\n");
                }
            }

            file.Close();
        }
        catch (Exception exception)
        {
            Console.WriteLine(exception.Message);
        }
    }

    #endregion Methods
}
}
