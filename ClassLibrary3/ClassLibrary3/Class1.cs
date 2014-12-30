using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using System.Collections.Concurrent;
using ExcelDna.Integration;

namespace TpVBA
{
    public class Donnees
    {

        public Donnees(double S, double T, double K, double V, double r)
        {

            nomOption = "vanilla";
            prix = S;
            maturite = T;
            strike = K;
            taux = r / 100;
            volatilite = V / 100;
            call = 0.0;
            put = 0.0;
            payoff = 0.0;
        }



        public static string nomOption { get; set; }
        public static double prix { get; set; }
        public static double maturite { get; set; }
        public static double strike { get; set; }
        public static double taux { get; set; }
        public static double volatilite { get; set; }
        public static double call { get; set; }
        public static double put { get; set; }
        public static double payoff { get; set; }

        public static double CalculCall(double S, double K, double T, double V, double r)
        {
            prix = S;
            maturite = T;
            strike = K;
            volatilite = V / 100;
            taux = r / 100;
            var normal = new MathNet.Numerics.Distributions.Normal(0.0, 1.0);
            var d1 = ((Math.Log(prix / strike) + maturite * (taux + ((volatilite * volatilite) / 2)))) / (volatilite * Math.Sqrt(maturite));
            var d2 = d1 - volatilite * Math.Sqrt(maturite);
            call = prix * normal.CumulativeDistribution(d1) - strike * Math.Exp(-taux * maturite) * normal.CumulativeDistribution(d2);
            return call;
        }
    }

    // Génération d'un nombre aléatoire suivant la loi Gaussienne. 
    public class BoxMuller
    {

        private double? maValeurActuelle;
        private double tempValeur;
        private readonly Random nombreAleatoire;

        public BoxMuller()
        {
            nombreAleatoire = new Random();
        }

        public double Aleatoire()
        {
            if (maValeurActuelle.HasValue)
            {
                maValeurActuelle = null;
                return tempValeur;
            }
            var valeurSauvegarde = Math.Sqrt(-2 * Math.Log(nombreAleatoire.NextDouble())) * Math.Cos(2.0 * Math.PI * nombreAleatoire.NextDouble());
            var valeurRetourne = Math.Sqrt(-2 * Math.Log(nombreAleatoire.NextDouble())) * Math.Sin(2.0 * Math.PI * nombreAleatoire.NextDouble());

            maValeurActuelle = valeurSauvegarde;
            tempValeur = valeurSauvegarde;

            return valeurRetourne;
        }
    }

    public class Pricer
    {

        public static ConcurrentBag<double>[] CreationMonteCarlos(double[] tabSO, double[] expDrift, Vector<double>[] x, double[] racineVolatilite, int nbr, int nbrOption)
        {
            var results = new ConcurrentBag<double>[nbrOption];
            for (int j = 0; j < nbrOption; ++j)
            {
                results[j] = new ConcurrentBag<double>();
                Parallel.For(0, nbr, i =>
                {
                    var S = tabSO[j] * Math.Exp(expDrift[j] + x[i].At(j) * racineVolatilite[j]);
                    results[j].Add(S);
                });
            }
            return results;
        }
        public static object CalculPricerCall(double spotDepart, double strike, double maturiteDepart, double volatiliteDepart, double tauxDepart, int nbr)
        {
            return ExcelAsyncUtil.Run("WebSnippetAsync", new object[] { spotDepart, strike, maturiteDepart, volatiliteDepart, tauxDepart, nbr },
            delegate
            {

                var aleaNorm = new BoxMuller();
                var nbrAleatoire = 0.0;
                var tabSpot = new double[] { spotDepart };
                var racineVolatilite = new double[1];
                var expDrift = new double[1];
                var x = new Vector<double>[nbr + 1];
                var resultat = new ConcurrentBag<double>[1];
                var payoff = new List<double>();
                var res = new double[1];


                expDrift[0] = (tauxDepart / 100.0 - (volatiliteDepart / 100.0 * volatiliteDepart / 100.0) / 2.0) * maturiteDepart;
                racineVolatilite[0] = Math.Sqrt(maturiteDepart) * (volatiliteDepart / 100.0);

                for (int i = 0; i < nbr + 1; ++i)
                {
                    x[i] = Vector<double>.Build.Dense(1);
                    nbrAleatoire = aleaNorm.Aleatoire();
                    x[i].At(0, nbrAleatoire);
                }

                resultat = CreationMonteCarlos(tabSpot, expDrift, x, racineVolatilite, nbr, 1);

                foreach (double element in resultat[0])
                {
                    payoff.Add(Math.Max(0, element - strike));
                }

                res[0] = payoff.Average() * Math.Exp(-tauxDepart * maturiteDepart / 100.0);

                return res[0];
            });
        }
    }
}
