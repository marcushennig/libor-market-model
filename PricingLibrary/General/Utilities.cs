using System;
using System.Collections.Generic;
using System.Linq;

namespace PricingLibrary.General
{
    public static class Utilities
    {
        #region Miscellaneous

        /// <summary>
        /// Draw progress bar in console 
        /// </summary>
        /// <param name="progress"></param>
        /// <param name="total"></param>
        public static void DrawTextProgressBar(int progress, int total, string info)
        {
            //draw empty progress bar
            Console.CursorLeft = 0;
            Console.Write("["); //start
            Console.CursorLeft = 32;
            Console.Write("]"); //end
            Console.CursorLeft = 1;

            var onechunk = 30.0f / total;

            //draw filled part
            var position = 1;
            for (var i = 0; i < onechunk * progress; i++)
            {
                Console.BackgroundColor = ConsoleColor.DarkGreen;
                Console.CursorLeft = position++;
                Console.Write(" ");
            }

            //draw unfilled part
            for (var i = position; i <= 31; i++)
            {
                Console.BackgroundColor = ConsoleColor.Green;
                Console.CursorLeft = position++;
                Console.Write(" ");
            }

            //draw totals
            Console.CursorLeft = 35;
            Console.BackgroundColor = ConsoleColor.Black;
            Console.Write("{0}: {1} of {2}        ", info, progress, total); //blanks at the end remove any excess
        }

        #endregion Miscellaneous

        #region Time structure

        /// <summary>
        /// Computes the tenor structure of
        /// </summary>
        /// <param name="startDate"></param>
        /// <param name="endDate"></param>
        /// <param name="tenor">Tenor in units of [months]</param>
        /// <returns></returns>
        public static DateTime[] TenorStructure(DateTime startDate, DateTime endDate, int tenor)
        {
            var tenorStructure = new List<DateTime>();

            var t = startDate;
            while (t <= endDate )
            {
                tenorStructure.Add(t);
                t = t.AddMonths(tenor);
            }

            if (tenorStructure.Last() != endDate)
            {
                tenorStructure.Add(endDate);
            }

            return tenorStructure.ToArray();
        }


        /// <summary>
        /// Inverse operation of ConvertToDayCountFraction
        /// </summary>
        /// <param name="startDate"></param>
        /// <param name="dayCountFraction"></param>
        /// <returns></returns>
        public static DateTime ConvertFromDayCountFraction(DateTime startDate, double dayCountFraction)
        {
            return startDate.AddDays( 360 * dayCountFraction);
        }
        /// <summary>
        /// Inverse operation of ConvertToDayCountFraction
        /// </summary>
        /// <param name="startDate"></param>
        /// <param name="dayCountFractions"></param>
        /// <returns></returns>
        public static DateTime[] ConvertFromDayCountFraction(DateTime startDate, double[] dayCountFractions)
        {
            return dayCountFractions.Select(dayCountFraction => ConvertFromDayCountFraction(startDate, dayCountFraction))
                                    .ToArray();
        }
        /// <summary>
        /// act / 360
        /// Actual number of days between start and end date divided by 360
        /// </summary>
        /// <param name="startDate"></param>
        /// <param name="endDate"></param>
        /// <returns></returns>
        public static double ConvertToDayCountFraction(DateTime startDate, DateTime endDate)
        {
            return (endDate - startDate).TotalDays / 360;
        }

        /// <summary>
        /// act / 360
        /// Actual number of days between start and end date divided by 360
        /// </summary>
        /// <param name="startDate"></param>
        /// <param name="endDates"></param>
        /// <returns></returns>
        public static double[] ConvertToDayCountFraction(DateTime startDate, DateTime[] endDates)
        {
            return endDates.Select(endDate => ConvertToDayCountFraction(startDate, endDate))
                           .ToArray();
        }

        #endregion Time structure
    }
}
