using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Numerics;

namespace ConsoleApplication1
{

    

    public class Pair<T1, T2, T3, T4, T5>
    {
        public Pair(T1 x1, T2 x2, T3 x3, T4 x4, T5 x5)
        {
            Item1 = x1;
            Item2 = x2;
            Item3 = x3;
            Item4 = x4;
            Item5 = x5;

        }
        public T1 Item1 { get; set; }
        public T2 Item2 { get; set; }

        public T3 Item3 { get; set; }


        public T4 Item4 { get; set; }


        public T5 Item5 { get; set; }
    }

    class Program
    {




        static public void depthExperiment3(int nMin, int nMax, int dTarget, int seed)
        {
            Dictionary<string, string> plots = new Dictionary<string, string>();

            plots.Add("Argon2iBGreedy", "\\addplot coordinates {");
            plots.Add("Argon2iBNoGreedy", "\\addplot coordinates {");
            plots.Add("Argon2iBRatio", "\\addplot coordinates {");
            // plots.Add("DRSampleHybrid", "\\addplot coordinates {");
            // plots.Add("Deterministic", "\\addplot coordinates {");
            plots.Add("DRSampleGreedy", "\\addplot coordinates {");
            plots.Add("DRSampleNoGreedy", "\\addplot coordinates {");
            plots.Add("DRSampleRatio", "\\addplot coordinates {");

            plots.Add("DRSampleBRGGreedy", "\\addplot coordinates {");
            plots.Add("DRSampleBRGNoGreedy", "\\addplot coordinates {");
            plots.Add("DRSampleBRGRatio", "\\addplot coordinates {");
            //plots.Add("DRSampleBitRev", "\\addplot coordinates {");
            for (int n = nMin; n <= nMax; n++)
            {
                int radius = 8;
                double ratio = Math.Pow(2, (n - 18) / 2.0);
                int nNodes = (int)Math.Round(4000 * ratio);//increase @0419
                List<Pair<int[], string, double, int, int>> DAGs = new List<Pair<int[], string, double, int, int>>();
                DAGs.Add(new Pair<int[], string, double, int, int>(Helpers.Argon2iBDAG(n, seed), "Argon2iB", 0.0, 0, 0));
                DAGs.Add(new Pair<int[], string, double, int, int>(Helpers.DRSampleDAG(n, seed), "DRSample", 0.0, 0, 0));
                DAGs.Add(new Pair<int[], string, double, int, int>(Helpers.BitSwapDRSample(n, seed), "DRSampleBRG", 0.0, 0, 0));

                Console.WriteLine("Built The Graphs");
                foreach (Pair<int[],string,double,int,int> P in DAGs)
                {
                    int[] DAG = P.Item1;
                    Tuple<int,int,int,bool[]> T1 = Helpers.BestDepthReducingSetExcludeGreedy(DAG, dTarget);
                    int e1 = T1.Item2;
                    plots[P.Item2 + "NoGreedy"] += "(" + n + "," + e1 + ") ";
                    Tuple<int, int, int, bool[]> T2 = Greedy.myGreedyDRSetsApxNoSortingNew(DAG, dTarget,radius, nNodes, dTarget);
                    int e2 = T2.Item2;
                    plots[P.Item2 + "Greedy"] += "(" + n + "," + e2 + ") ";
                    double theRatio = e1 *1.0 /(1.0* e2);
                    plots[P.Item2+"Ratio"] += "(" + n + "," + theRatio + ") ";
                    Console.WriteLine(P.Item2 + "Greedy: (" + n + "," + e2 + ")" );
                    Console.WriteLine(P.Item2 + "NoGreedy: (" + n + "," + e1 + ")");
                    Console.WriteLine("Ratio: " + theRatio);
                }

            }
            foreach (KeyValuePair<string, string> pair in plots)
            {
                Console.WriteLine("% " + pair.Key + ", d = " + dTarget);
                Console.WriteLine(pair.Value + "};");
            }
        }



        static public void attackExperiment(int nMin, int nMax, int seed, double R = 3000.0)
        {
            Dictionary<string, string> plots = new Dictionary<string, string>();
            plots.Add("Argon2iB", "\\addplot coordinates {");
            plots.Add("DRSample", "\\addplot coordinates {");
            plots.Add("DRSampleBitRev", "\\addplot coordinates {");// @Lu BitSwapDRSample

            int[] ib_depth = new int[] {7, 7, 7, 7, 7, 60313, 29824, 30835, 26316, 24930, 31201};
            int[] dr_depth = new int[] {8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8};
            int[] drBR_depth = new int[] {5, 5, 5, 5, 64776, 64528, 16656, 30017, 117263, 46135, 352309};

            for (int n=nMin; n <= nMax; n++)
            {
                Console.WriteLine("n = " + n + ",     honestCost = "+ HonestECost(n));
                /// int[] encodes the DAG
                /// string to print
                /// String the name of the DAG (e.g., Argon2iB)
                /// double the minimum energy cost to pebble the DAG
                /// int e parameter e under optimal attack
                /// int d parameter d under optimal attack
               
                List<Pair<int[],string,double, int,int>> DAGs = new List<Pair<int[], string, double, int,int>>();
                 DAGs.Add(new Pair<int[], string, double,  int, int>(Helpers.Argon2iBDAG(n, seed),  "Argon2iB", double.PositiveInfinity, 0,0));

                int dMid = 0;
                int radius = 8;
                double ratio = Math.Pow(2, (n - 18) / 2.0);
                //int nNodes = (int)Math.Round(400 * ratio);
                int nNodes = (int)Math.Round(400 * ratio);//increase @0419
                foreach (Pair<int[], string,double,int,int> P in DAGs)
                {
                    int[] DAG = P.Item1;
                    string name = P.Item2;
                    if (string.Equals(name, "Argon2iB"))
                        dMid = (int)Math.Ceiling(Math.Log(ib_depth[n - 14], 2.0));
                    else if (string.Equals(name, "DRSample"))
                        dMid = (int)Math.Ceiling(Math.Log(dr_depth[n - 14], 2.0)); 
                    else if (string.Equals(name, "DRSampleBitRev"))
                        dMid = (int)Math.Ceiling(Math.Log(drBR_depth[n - 14], 2.0));

                    if (dMid == 0)
                       System.Environment.Exit(1);
                    
                    
                    for (int dPow=dMid-1; dPow <= dMid+1; dPow++)
                    {
                        int dTarget = (1 << dPow);
                        Tuple<int, int, int, bool[]> T0;
                        if (dPow > 0)
                        {
                            T0 = Helpers.BestDepthReducingSetExcludeGreedy(DAG,dTarget);
                        }
                        else
                        {
                            T0 = Greedy.myGreedyDRSetsApxNoSorting(DAG, dTarget, radius, nNodes);
                        }
                        // Tuple<int, int, int, bool[]> T0 = Helpers.BestDepthReducingSetExcludeGreedy(P.Item1, dTarget);
                        int e = T0.Item2;
                        double minCost = P.Item3;
                        int d = T0.Item1;
                        int N = 1 << n;
                        int gMax = (int)Math.Min(N, 25* Math.Max(Math.Sqrt(N * 1.0 * d), e));
                        int gMin =  (int)Math.Max(dTarget, Math.Min(Math.Sqrt(N * 1.0 * d), e)/5);
                        int diff = gMax - gMin;
                        int[] g = new int[201];
                        for (int i = 0; i < g.Length; i++) g[i] = gMin + (int) ((i * diff) / 200.0);
                        //int[] g = new int[] { gMin, gMin + diff / 20, gMin + 2 * diff / 20, gMin + 3 * diff / 20, gMin + 4 * diff / 20, gMin + 5 * diff / 20, gMin + 6 * diff / 20, gMin + 7 * diff / 20, gMin + 8 * diff / 20, gMin + 9 * diff / 20, gMin + 10 * diff / 20, gMin + 11 * diff / 20, gMin + 12 * diff / 20, gMin + 13 * diff / 20, gMin + 14 * diff / 20, gMin + 15 * diff / 20, gMin + 16 * diff / 20, gMin + 17 * diff / 20, gMin + 17 * diff / 20, gMin + 18 * diff / 20, gMin + 19 * diff / 20, gMin + diff };
                        double cost = Helpers.EAttackCostUpperBound(N, DAG, T0.Item4, e, d, g, R);
                        //cost = Math.Min(Helpers.GreedyPebblingAttack(N, DAG, R),cost);
                        //Console.WriteLine("e = " + e);
                        Console.Write(".");
                        Console.WriteLine("Cost (log_2 dTgt =  " + dPow + "):" + cost);
                        Console.WriteLine("e = " + e);
                        Console.WriteLine("Actual Depth: " + Helpers.Depth(DAG, T0.Item4));
                        if (cost <= minCost)
                        {
                            minCost = cost;
                            P.Item3 = minCost;
                            P.Item4 = e;
                            P.Item5 = d;
                        }
                        minCost = Math.Min(cost, minCost);


                    }


          

                }
                Console.WriteLine();
                foreach (Pair<int[], string, double, int, int> P in DAGs)
                {
                    Console.WriteLine(P.Item2 + ": (" + n + "," + HonestECost(n) / P.Item3 + ")    " + "e=" + P.Item4 + "    d=" + P.Item5);
                    plots[P.Item2]+= "(" + n + "," + HonestECost(n) / P.Item3 + ") ";
                }

            }
            foreach(KeyValuePair<string,string> pair in plots)
            {
                Console.WriteLine("% "+ pair.Key);
                Console.WriteLine(pair.Value + "};");
            }

        }


        static double HonestECost(int n, double R=3000.0)
        {
            return (1 << n) * 1.0 * (1 << n) / 2 + (1 << n) * R;
        }

        static void Main(string[] args)
        {

            int seed = 123456;
            //depthExperiment3(21, 24, 16, seed);
            //attackExperiment(23, 23, seed);
            Console.ReadLine();
            //Console.WriteLine(Greedy.depth(DAG,S));

            /*
            


            
            Console.WriteLine("done");
            Console.ReadLine();
          */
        }
    }
}
