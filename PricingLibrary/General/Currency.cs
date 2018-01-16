namespace PricingLibrary.General
{
    public class Currency
    {
        #region Properties

        /// <summary>
        /// Name of the currency (3-letter code)
        /// </summary>
        private string Name;

        public static Currency EUR = new Currency("EUR");
        public static Currency USD = new Currency("USD");

        #endregion Properties

        #region Constructor 

        /// <summary>
        /// Constructor 
        /// </summary>
        /// <param name="name"></param>
        private Currency(string name)
        {
            Name = name;
        }

        #endregion Constructor 

        #region Methods

        /// <summary>
        /// String representation
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            return Name;
        }

        #endregion Methods
    }
}
