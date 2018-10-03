using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Globalization;
using System.Numerics;
using System.IO;

namespace ConsoleApplication1
{
    class Helpers
    {
        /*
         *  Constructors for several DAGs.
         * 
         */
        public static int[] Argon2iADAG(int n, int seed)
        {
            int[] DAG = new int[1 << n];
            DAG[0] = 0;
            Random r = new Random(seed);
            for (int i = 1; i < (1 << n); i++)
            {
                int j = r.Next() % i;
                DAG[i] = j;
            }
            return DAG;
        }

        public static int[] Argon2iADAGTwoPasses(int n, int seed)
        {
            int[] DAG = new int[1 << n];
            DAG[0] = 0; DAG[1] = 0;
            Random r = new Random(seed);
            for (int i = 2; i < (1 << n); i++)
            {
                int MaxDist = Math.Min(i-1, (1 << (n - 1)));
                int j = r.Next() % MaxDist;
                DAG[i] = i-j-1;
            }
            return DAG;
        }

        public static int[] Argon2iBDAG(int n, int seed)
        {
            int[] DAG = new int[1 << n];
            DAG[0] = 0;
            Random r = new Random(seed);
            for (int i = 1; i < (1 << n); i++)
            {
                double j = (r.Next() % (1 << 30)) * 1.0 / (1.0 * (1 << 30));
                int k = (int)Math.Floor(i * (1.0 - j * j));
                DAG[i] = k;

            }
            return DAG;
        }


        public static int[] Argon2iBDAGTwoPasses(int n, int seed)
        {
            int[] DAG = new int[1 << n];
            DAG[0] = 0;
            Random r = new Random(seed);
            for (int i = 1; i < (1 << n); i++)
            {
                int maxBackDist = Math.Min(i, (1 << (n - 1)));
                double j = (r.Next() % (1 << 30)) * 1.0 / (1.0 * (1 << 30));
                int k = (int)Math.Floor(maxBackDist * (1.0 - j * j));
                DAG[i] = k;

            }
            return DAG;
        }

        public static int[] Argon2iBDAGwithBug(int n, int seed)
        {
            int[] DAG = new int[1 << n];
            DAG[0] = 0;
            Random r = new Random(seed);
            for (int i = 1; i < (1 << n); i++)
            {
                double j = ((r.Next() % (1 << 29)) + (1<<29)) * 1.0 / (1.0 * (1 << 30));
                int k = (int)Math.Floor(i * (1.0 - j * j));
                DAG[i] = k;

            }
            return DAG;
        }
        /// <summary>
        /// DRSampleHybrid
        /// </summary>
        /// <param name="n"></param>
        /// <param name="seed"></param>
        /// <returns></returns>
        public static int[] CustomDAG(int n)
        {
            int[] DAG = new int[1 << n];
            DAG[0] = 0;
            for (int i = 1; i < 1 << n; i++)
            {
                int j = i % n;
                int k = i % (1 << j);
                Random r = new Random(k);
                int jj = Math.Min(i, 1 << (j + 1));
                int backDist = r.Next() % jj;
                DAG[i] = i - backDist - 1;
            }
            return DAG;
        }
        public static int[] DRSampleHybridDAG(int n, int seed)
        {
            Random r = new Random(seed);
            int[] DAG = new int[1 << n];
            DAG[0] = 0;
            for (int i = 1; i < 1 << n; i++)
            {
                int j = i % n;
                int k = i % (1 << j);

                int jj = Math.Min(i, 1 << (j + 1));
                int backDist = r.Next(jj >> 1, jj + 1);
                DAG[i] = i - backDist;
            }
            return DAG;
        }

        /// <summary>
        /// Deterministic
        /// </summary>
        /// <param name="n"></param>
        /// <param name="seed"></param>
        /// <returns></returns>
        public static int[] DeterministicDAG(int n)
        {
            int[] DAG = new int[1 << n];
            DAG[0] = 0;
            for (int i = 1; i < 1 << n; i++)
            {
                int j = i % n;


                DAG[i] = i - Math.Min(i, (1 << j));
            }
            return DAG;
        }
        /// <summary>
        /// DRSample
        /// </summary>
        /// <param name="n"></param>
        /// <param name="seed"></param>
        /// <returns></returns>
        public static int[] DRSampleDAG(int n, int seed)
        {
            Random r = new Random(seed);
            int[] DAG = new int[1 << n];
            DAG[0] = 0;
            DAG[1] = 1;
            for (int i = 2; i < 1 << n; i++)
            {
                int logi = (int)Math.Floor(Math.Log(i, 2.0));
                int j = r.Next() % logi;

                int jj = Math.Min(i, 1 << (j + 1));
                int backDist = r.Next(Math.Max(jj >> 1, 2), jj + 1);
                DAG[i] = i - backDist;
            }
            return DAG;
        }

        /// <summary>
        /// DRSample
        /// </summary>
        /// <param name="n"></param>
        /// <param name="seed"></param>
        /// <returns></returns>
        public static int[] DRSampleDAGTwoPasses(int n, int seed)
        {
            Random r = new Random(seed);
            int[] DAG = new int[1 << n];
            DAG[0] = 0;
            DAG[1] = 1;
            for (int i = 2; i < 1 << n; i++)
            {
                int maxDist = Math.Min(i, (1 << (n - 1)));
                int logi = (int)Math.Floor(Math.Log(maxDist, 2.0));
                int j = r.Next() % logi;

                int jj = Math.Min(i, 1 << (j + 1));
                int backDist = r.Next(Math.Max(jj >> 1, 2), jj + 1);
                DAG[i] = i - backDist;
            }
            return DAG;
        }
        /// <summary>
        /// DRSampleModified
        /// </summary>
        /// <param name="n"></param>
        /// <param name="seed"></param>
        /// <returns></returns>
        public static int[] DRSampleModifiedDAG(int n, int seed)
        {
            Random r = new Random(seed);
            int[] DAG = new int[1 << n];
            DAG[0] = 0;
            DAG[1] = 1;
            for (int i = 2; i < 1 << n; i++)
            {
                int logi = (int)Math.Floor(Math.Log(i, 2.0));
                
                int j = r.Next() % logi;
                

                int jj = (r.Next() % 2 ==1)? Math.Min(i, 1 << (j + 1)): i;
                int backDist = r.Next(Math.Max(jj >> 1, 2), jj + 1);
                DAG[i] = i - backDist;
            }
            return DAG;
        }


        /// <summary>
        /// DRSampleModified2Pass
        /// </summary>
        /// <param name="n"></param>
        /// <param name="seed"></param>
        /// <returns></returns>
        public static int[] DRSampleModifiedDAGWith2Passes(int n, int seed)
        {
            Random r = new Random(seed);
            int[] DAG = new int[1 << n];
            DAG[0] = 0;
            DAG[1] = 1;
            for (int i = 2; i < (1 << n); i++)
            {
                int maxDist = Math.Min(i, (1 << (n-1)));
                int logi = (int)Math.Floor(Math.Log(maxDist, 2.0));

                int j = r.Next() % logi;


                int jj = (r.Next() % 2 == 1) ? Math.Min(Math.Min(maxDist, (1<<(n-1))), 1 << (j + 1)) : maxDist;
                int backDist = r.Next(Math.Max(jj >> 1, 2), jj + 1);
                DAG[i] = i - backDist;
            }
            return DAG;
        }

        /// <summary>
        /// CPCSampleHybrid
        /// </summary>
        /// <param name="n"></param>
        /// <param name="seed"></param>
        /// <returns></returns>
        public static int[] CPCSampleHybridDAG(int n,int seed)
        {
            Random r = new Random(seed);
            int[] DAG = new int[1 << n];
            DAG[0] = 0;
            DAG[1] = 1;
            for (int i = 2; i < (1 << (n - 1)); i++)
            {
                int logi = (int)Math.Floor(Math.Log(i, 2.0));
                int j = r.Next() % logi;

                int jj = Math.Min(i, 1 << (j + 1));
                int backDist = r.Next(Math.Max(jj >> 1, 2), jj + 1);
                DAG[i] = i - backDist;
            }
  

            for (int i = (1 << (n - 1)); i < (1 << (n)); i++)
            {
                if (i % 2 == 1)
                {
                    DAG[i] = r.Next() % (1 << (n - 1));
                }
                else
                {
                    int logi = (int)Math.Floor(Math.Log(i, 2.0));
                    int j = r.Next() % logi;

                    int jj = Math.Min(i, 1 << (j + 1));
                    int backDist = r.Next(Math.Max(jj >> 1, 2), jj + 1);
                    DAG[i] = i - backDist;
                }
            }
            return DAG;
        }
        /// <summary>
        /// CPCSample
        /// </summary>
        /// <param name="n"></param>
        /// <param name="seed"></param>
        /// <returns></returns>
        public static int[] CPCSampleDAG(int n, double c, int seed)
        {
            Random r = new Random(seed);
            int[] DAG = new int[1 << n];
            DAG[0] = 0;
            DAG[1] = 1;
            for (int i = 2; i < (1 << (n - 1)); i++)
            {
                int logi = (int)Math.Floor(Math.Log(i, 2.0));
                int j = r.Next() % logi;

                int jj = Math.Min(i, 1 << (j + 1));
                int backDist = r.Next(Math.Max(jj >> 1, 2), jj + 1);
                DAG[i] = i - backDist;
            }
            int m = (int) Math.Floor(c *n);
            int nOverm = (int)Math.Floor((1<<n-1) * 1.0 / (1.0 * m));

            for (int i = (1 << (n - 1)); i < (1 << (n)); i++)
            {
                int b = (i - (1 << (n - 1))) % nOverm;
                DAG[i] = b*m;
               
            }
            return DAG;
        }


        public static int[] CustomDAG5(int n)
        {
            int[] DAG = new int[1 << n];
            DAG[0] = 0;
            for (int i = 1; i < 1 << n; i++)
            {
                int j = i % n;
                int k = i % (1 << j);
                Random r = new Random(k);
                int jj = Math.Min(i, 1 << (j + 1));
                int backDist = r.Next(jj >> 1, jj);
                DAG[i] = i - backDist - 1;
            }
            return DAG;
        }

        /// <summary>
        /// Constructs the DAG for the bit-swapping graph built on top of DRSample
        /// </summary>
        /// <param name="n">DAG Size</param>
        /// <param name="seed">PRNG Seed</param>
        public static int[] BitSwapDRSample(int n, int seed)
        {
            if(n > 31 || n < 2) throw new ArgumentException("n must be between 2 and 31");

            // Initialize
            Random r = new Random(seed);
            int[] DAG = new int[1 << n];
            DAG[0] = 0;
            DAG[1] = 1;

            // First half of DAG is DRSample
            for (int i = 2; i < 1 << (n-1); i++)
            {
                int logi = (int)Math.Floor(Math.Log(i, 2.0));
                int j = r.Next() % logi;

                int jj = Math.Min(i, 1 << (j + 1));
                int backDist = r.Next(Math.Max(jj >> 1, 2), jj + 1);
                DAG[i] = i - backDist;
            }

            // Construct the high order bit mask
            int mask = 0x00000001;
            for(int i = 1; i < (n-1); i++){
                mask <<= 1;
                mask |= 0x00000001;
            }

            // Construct the second half of the list

            // Swapping indices, temps, helpers
            int x;


            for(int k = (1 << (n-1)); k < 1 << n; k++){
                //Mask away all the high order bits
                int b = k & mask;

                // Swap remaining bits within their space
                int i = 0;
                int j = n-2;
                while(i < j)
                {
                    x = ((b >> i) ^ (b >> j)) & ((1 << 1) - 1);
                    b = b ^ ((x << i) | (x << j));
                    i++; j--;
                }

                DAG[k] = b;
            }
            
            return DAG;
        }

        /// <summary>
        /// Saves a copy of
        /// </summary>
        /// <param name="DRS">The depth-reducing set to save</param>
        /// <param name="filename">The filename where you want to save it</param>
        public static void saveDRS(bool[] DRS, string filename)
        {
            FileStream F = new FileStream(filename, FileMode.Open, FileAccess.Read, FileShare.Read);
            writeBooleans(F, DRS);
            F.Close();
        }

        /// <summary>
        /// Load a depth reducing set
        /// </summary>
        /// <param name="n">Size of array to allocate</param>
        /// <param name="filename">File where the DRS is stored</param>
        public static bool[] loadDRS(int n, string filename)
        {
            bool[] DRS = new bool[n];
            FileStream F = new FileStream(filename, FileMode.Open, FileAccess.Read, FileShare.Read);
            readBooleans(F, DRS);
            F.Close();
            return DRS;
        }


        public static void writeBooleans(System.IO.FileStream outfile, bool[] ar)
        {
            for (int i = 0; i < ar.Length; i += 8) {
                int b = 0;
                for (int j = Math.Min(i + 7, ar.Length-1); j >= i; j--) {
                    b = (b << 1) | (ar[j] ? 1 : 0);
                }
                outfile.WriteByte((byte)b);
            }
        }   

        private static void readBooleans(System.IO.FileStream infile, bool[] ar)
        {
            for (int i = 0; i < ar.Length; i += 8) {
                int b = infile.ReadByte();
                if (b < 0) throw new IndexOutOfRangeException();
                for (int j = i; j < i + 8 && j < ar.Length; j++)
                {
                    ar[j] = (b & 1) != 0;
                    b >>= 1;
                }
            }
        }

        /// <summary>
        /// Upper Bounds the cumulative cost of running the Alwen-Blocki Parallel Pebbling Attack on G given depth-reducing set S
        /// such that depth(G-S) <= d and |S| <= e.
        /// </summary>
        /// <param name="n">number of nodes in graph</param>
        /// <param name="DAG">input graph</param>
        /// <param name="S">depth reducing set</param>
        /// <param name="e">size of depth-reducing set |S|</param>
        /// <param name="d">upper bound on depth(G-S)</param>
        /// <param name="R">core memory area ratio (typically 3000.0). Size of core evaluating hash function divided by area to store one hash value in memory</param>
        /// <returns>double: Upper bound on aAT(G) </returns>
        public static double EAttackCostUpperBound(int n, int[] DAG, bool[] S, int e, int d, double R)
        {

            // Stores the depth of every node in DAG-S
            int[] depths = new int[n];
            int[] lastUsed = new int[n];
            int[] maxDepth = new int[n];

            double attackCost = 0.0;
            //Upper bound cost on parents of nodes in light phases + cost on S
            for (int i = 0; i < n; i++)
            {
                depths[i] = S[i] ? 0 : 1;
                maxDepth[i] = 0;

                // cost to pebble node i+1 for the first time
                attackCost += R;
                // Add Cost of pebbles on S
                if (S[i])
                {
                    // We keep a pebble on node i+1 in rounds i+1,...,n
                    attackCost += 1.0 * n - 1.0 * i;
                }


                // We start a balloon phase with <= e additional (not in S) pebbles on parents of next e nodes 
                // during a light phase we discard these pebbles as soon as we are done to during round j< g of a light phase
                // we have at most g-j pebbles on parents. The average number of pebbles per round is (g+(g-1)+...+1)/g = (g+1)/2.
                attackCost += e * 0.5 + 0.5;
            }

            // Count the number of nodes at depth d in G
            int[] depthCounts = new int[n];
            // Count the number of nodes whoses topologically last child has depth d
            int[] lastDepthCount = new int[n];

            for (int i = 1; i < DAG.Length; i++)
            {
                depthCounts[i] = 0;
                maxDepth[i] = 0;
                lastDepthCount[i] = 0;


            }
            depths[0] = (S[0]) ? 0 : 1;
            maxDepth[0] = depths[0];
            depthCounts[depths[0]] += 1;


            //Upper Bound costs of balloon phases
            for (int i = 1; i < DAG.Length; i++)
            {
                depths[i] = S[i] ? 0 : 1;
                if (S[i - 1] == false && S[i] == false) depths[i] = Math.Max(depths[i], depths[i - 1] + 1);
                if (S[DAG[i]] == false && S[i] == false) depths[i] = Math.Max(depths[i], depths[DAG[i]] + 1);
                depthCounts[depths[i]] += 1;
                lastUsed[i - 1] = depths[i];
                lastDepthCount[depths[i]]++;
                if (DAG[i] < i - 1 && depths[i] > lastUsed[DAG[i]])
                {
                    int theDepth = lastUsed[DAG[i]];
                    lastDepthCount[theDepth]--;
                    lastUsed[DAG[i]] = depths[i];
                    lastDepthCount[depths[i]]++;
                }

                maxDepth[i] = Math.Max(maxDepth[i - 1], depths[i]);
                if (0 == (i % e))
                {
                    int curDepth = maxDepth[i];
                    // Cost balloon phase <= curDepth*i - sum_{j=1}^{curDepth} \sum_{k =j+1}^{curdepth} depthCounts[k]
                    double savings = 0.0;
                    double savings2 = 0.0;
                    for (int j = 1; j <= curDepth; j++)
                    {
                        // We don't have pebbles on these nodes for j-1 of the balloon phase rounds
                        savings += (j - 1) * 1.0 * depthCounts[j];

                        //These pebbles can be removed (curDepth - j) rounds before the baloon phase finishes
                        savings2 += (curDepth - j) * 1.0 * lastDepthCount[j];
                        //Cost to place re-place pebble on node during balloon phase
                        attackCost += depthCounts[j] * R;
                    }

                    attackCost += curDepth * 1.0 * i - savings - Math.Max(0, savings2 - curDepth * e);
                }
            }


            return attackCost;
        }
        /// <summary>
        /// Simulates AB16 pebbling attack given depth-reducing set S. Searches for best attack parameter g in array.
        /// </summary>
        /// <param name="n">number of nodes in graph</param>
        /// <param name="DAG">input graph</param>
        /// <param name="S">depth reducing set</param>
        /// <param name="e">size of depth-reducing set |S|</param>
        /// <param name="d">upper bound on depth(G-S)</param>
        /// <param name="gs">attack parameters to try. Each g is length of light phase in AB16 attack.</param>
        /// <param name="R">core memory area ratio (typically 3000.0). Size of core evaluating hash function divided by area to store one hash value in memory</param>
        /// <returns>double: Upper bound on aAT(G) </returns>
        public static double EAttackCostUpperBound(int n, int[] DAG, bool[] S, int e, int d, int [] gs, double R)
        {
            int bestG = 0;
            double minCost = double.PositiveInfinity;
            for(int i=0;i<gs.Length;i++)
            {
                double newCost = EAttackCostUpperBound(n, DAG, S, e, d, gs[i], R);
                if (minCost > newCost)
                {
                    bestG = i;
                    minCost = newCost;
                }
            }
            //Console.WriteLine("   Optimal g = " + gs[bestG]);
            return minCost;
        }

        /// <summary>
        /// Simulates AB16 pebbling attack given depth-reducing set S. Searches for best attack parameter g in array.
        /// </summary>
        /// <param name="n">number of nodes in graph</param>
        /// <param name="DAG">input graph</param>
        /// <param name="S">depth reducing set</param>
        /// <param name="e">size of depth-reducing set |S|</param>
        /// <param name="d">upper bound on depth(G-S)</param>
        /// <param name="gs">attack parameters to try. Each g is length of light phase in AB16 attack.</param>
        /// <param name="R">core memory area ratio (typically 3000.0). Size of core evaluating hash function divided by area to store one hash value in memory</param>
        /// <returns>double: Upper bound on aAT(G) </returns>
        public static double EAttackCostUpperBoundTwoPass(int n, int[] DAG, bool[] S, int e, int d, int[] gs, double R)
        {
            int bestG = 0;
            double minCost = double.PositiveInfinity;
            for (int i = 0; i < gs.Length; i++)
            {
                double newCost = EAttackCostUpperBoundTwoPass(n, DAG, S, e, d, gs[i], R);
                if (minCost > newCost)
                {
                    bestG = i;
                    minCost = newCost;
                }
            }
            Console.WriteLine("   Optimal g = " + gs[bestG]);
            return minCost;
        }


        /// <summary>
        /// Similar to EAttackCostUpperBound, but includes the new GreedyPebblingAttack in addition to AB16 attacks. Returns best of both attacks.
        /// </summary>
        /// <param name="n">number of nodes in graph</param>
        /// <param name="DAG">input graph</param>
        /// <param name="S">depth reducing set</param>
        /// <param name="e">size of depth-reducing set |S|</param>
        /// <param name="d">upper bound on depth(G-S)</param>
        /// <param name="R">core memory area ratio (typically 3000.0). Size of core evaluating hash function divided by area to store one hash value in memory</param>
        /// <returns>double: Upper bound on aAT(G) </returns>
        public static double BestEAttackCostUpperBound(int n, int[] DAG, bool[] S, int e, int d, double R)
        {
            double minCost = GreedyPebblingAttack(n,DAG,R);
            minCost = Math.Min(minCost, EAttackCostUpperBound(n, DAG, S, e, d, R));

            return minCost;
        }


        /// <summary>
        /// Upper Bounds the cumulative cost of running the Alwen-Blocki Parallel Pebbling Attack on G given depth-reducing set S
        /// such that depth(G-S) <= d and |S| <= e.
        /// </summary>
        /// <param name="n">number of nodes in graph</param>
        /// <param name="DAG">input graph</param>
        /// <param name="S">depth reducing set</param>
        /// <param name="e">size of depth-reducing set |S|</param>
        /// <param name="d">upper bound on depth(G-S)</param>
        /// <param name="g">Tuneable attack Parameter for AB16 attack. Length of light phase.</param>
        /// <param name="R">core memory area ratio (typically 3000.0). Size of core evaluating hash function divided by area to store one hash value in memory</param>
        /// <returns>double: Upper bound on aAT(G) </returns>
        public static double EAttackCostUpperBound(int n, int[] DAG, bool[] S, int e, int d,int g, double R)
        {
            if (g <= d) return n * 0.5 * n + R * n;
            // Stores the depth of every node in DAG-S
            int[] depths = new int[n];
            int[] lastUsed = new int[n];
            int[] maxDepth = new int[n];

            double attackCost = 0.0;
            //Upper bound cost on parents of nodes in light phases + cost on S
            for (int i = 0; i < n; i++)
            {
                depths[i] = S[i] ? 0 : 1;
                maxDepth[i] = 0;

                // cost to pebble node i+1 for the first time
                attackCost += R;
                // Add Cost of pebbles on S
                if (S[i])
                {
                    // We keep a pebble on node i+1 in rounds i+1,...,n
                    attackCost += 1.0 * n - 1.0 * i;
                }


                // We start a balloon phase with <= g additional (not in S) pebbles on parents of next g nodes 
                // during a light phase we discard these pebbles as soon as we are done to during round j< g of a light phase
                // we have at most g-j pebbles on parents. The average number of pebbles per round is (g+(g-1)+...+1)/g = (g+1)/2.
                attackCost += g * 0.5 + 0.5;
            }

            // Count the number of nodes at depth d in G
            int[] depthCounts = new int[n];
            // Count the number of nodes whoses topologically last child has depth d
            int[] lastDepthCount = new int[n];

            for (int i = 1; i < DAG.Length; i++)
            {
                depthCounts[i] = 0;
                maxDepth[i] = 0;
                lastDepthCount[i] = 0;


            }
            depths[0] = (S[0]) ? 0 : 1;
            maxDepth[0] = depths[0];
            depthCounts[depths[0]] += 1;


            //Upper Bound costs of balloon phases
            for (int i = 1; i < DAG.Length; i++)
            {
                depths[i] = S[i] ? 0 : 1;
                if (S[i - 1] == false && S[i] == false) depths[i] = Math.Max(depths[i], depths[i - 1] + 1);
                if (S[DAG[i]] == false && S[i] == false) depths[i] = Math.Max(depths[i], depths[DAG[i]] + 1);
                depthCounts[depths[i]] += 1;
                lastUsed[i - 1] = depths[i];
                lastDepthCount[depths[i]]++;
                if (DAG[i] < i - 1 && depths[i] > lastUsed[DAG[i]])
                {
                    int theDepth = lastUsed[DAG[i]];
                    lastDepthCount[theDepth]--;
                    lastUsed[DAG[i]] = depths[i];
                    lastDepthCount[depths[i]]++;
                }

                maxDepth[i] = Math.Max(maxDepth[i - 1], depths[i]);
                if (0 == (i % g))
                {
                    int curDepth = maxDepth[i];
                    // Cost balloon phase <= curDepth*i - sum_{j=1}^{curDepth} \sum_{k =j+1}^{curdepth} depthCounts[k]
                    double savings = 0.0;
                    double savings2 = 0.0;
                    for (int j = 1; j <= curDepth; j++)
                    {
                        // We don't have pebbles on these nodes for j-1 of the balloon phase rounds
                        savings += (j - 1) * 1.0 * depthCounts[j];

                        //These pebbles can be removed (curDepth - j) rounds before the baloon phase finishes
                        savings2 += (curDepth - j) * 1.0 * lastDepthCount[j];
                        //Cost to place re-place pebble on node during balloon phase
                        attackCost += depthCounts[j] * R;
                    }

                    attackCost += curDepth * 1.0 * i - savings - Math.Max(0, savings2 - curDepth * g);
                }
            }


            return attackCost;
        }

        /// <summary>
        /// Upper Bounds the cumulative cost of running the Alwen-Blocki Parallel Pebbling Attack on G (two-pass) given depth-reducing set S
        /// such that depth(G-S) <= d and |S| <= e.
        /// </summary>
        /// <param name="n">number of nodes in graph</param>
        /// <param name="DAG">input graph</param>
        /// <param name="S">depth reducing set</param>
        /// <param name="e">size of depth-reducing set |S|</param>
        /// <param name="d">upper bound on depth(G-S)</param>
        /// <param name="g">Tuneable attack Parameter for AB16 attack. Length of light phase.</param>
        /// <param name="R">core memory area ratio (typically 3000.0). Size of core evaluating hash function divided by area to store one hash value in memory</param>
        /// <returns>double: Upper bound on aAT(G) </returns>
        public static double EAttackCostUpperBoundTwoPass(int n, int[] DAG, bool[] S, int e, int d, int g, double R)
        {
            if (g <= d) return (n/2) * 0.5 * (n/2) + (n/2)*1.0*(n/2) + R * n;
            // Stores the depth of every node in DAG-S
            int[] depths = new int[n];
            int[] lastUsed = new int[n];
            int[] maxDepth = new int[n];

            int lengthFirstPass = n / 2;
            double attackCost = 0.0;
            //Upper bound cost on parents of nodes in light phases + cost on S
            for (int i = 0; i < n; i++)
            {
                depths[i] = S[i] ? 0 : 1;
                maxDepth[i] = 0;

                // cost to pebble node i+1 for the first time
                attackCost += R;
                // Add Cost of pebbles on S
                if (S[i])
                {
                    // We keep a pebble on node i+1 in rounds i+1,...,n
                    attackCost += 1.0 * n - 1.0 * i;
                }


                // We start a balloon phase with <= g additional (not in S) pebbles on parents of next g nodes 
                // during a light phase we discard these pebbles as soon as we are done to during round j< g of a light phase
                // we have at most g-j pebbles on parents. The average number of pebbles per round is (g+(g-1)+...+1)/g = (g+1)/2.
                attackCost += g * 0.5 + 0.5;
            }

            // Count the number of nodes at depth d in G
            int[] depthCounts = new int[n];
            // Count the number of nodes whoses topologically last child has depth d
            int[] lastDepthCount = new int[n];

            for (int i = 1; i < DAG.Length; i++)
            {
                depthCounts[i] = 0;
                maxDepth[i] = 0;
                lastDepthCount[i] = 0;


            }
            depths[0] = (S[0]) ? 0 : 1;
            maxDepth[0] = depths[0];
            depthCounts[depths[0]] += 1;


            //Upper Bound costs of balloon phases
            for (int i = 1; i < DAG.Length; i++)
            {
                depths[i] = S[i] ? 0 : 1;
                if (S[i - 1] == false && S[i] == false) depths[i] = Math.Max(depths[i], depths[i - 1] + 1);
                if (S[DAG[i]] == false && S[i] == false) depths[i] = Math.Max(depths[i], depths[DAG[i]] + 1);
                if (i >= lengthFirstPass && S[i-lengthFirstPass]==false && S[i] == false) depths[i] = Math.Max(depths[i], depths[i - lengthFirstPass] + 1);
                depthCounts[depths[i]] += 1;
                lastUsed[i - 1] = depths[i];
                    lastDepthCount[depths[i]]++;
                if (DAG[i] < i - 1 && depths[i] > lastUsed[DAG[i]])
                {
                    int theDepth = lastUsed[DAG[i]];
                    lastDepthCount[theDepth]--;
                    lastUsed[DAG[i]] = depths[i];
                    lastDepthCount[depths[i]]++;
                }
                if (i >= lengthFirstPass && depths[i] > lastUsed[i - lengthFirstPass])
                {
                    int theDepth = lastUsed[i - lengthFirstPass];
                    lastDepthCount[theDepth]--;
                    lastUsed[i - lengthFirstPass] = Math.Max(lastUsed[i - lengthFirstPass], depths[i]);
                    lastDepthCount[depths[i]]++;
                }

                maxDepth[i] = Math.Max(maxDepth[i - 1], depths[i]);
                if (0 == (i % g))
                {
                    int curDepth = maxDepth[i];
                    // Cost balloon phase <= curDepth*i - sum_{j=1}^{curDepth} \sum_{k =j+1}^{curdepth} depthCounts[k]
                    double savings = 0.0;
                    double savings2 = 0.0;
                    for (int j = 1; j <= curDepth; j++)
                    {
                        // We don't have pebbles on these nodes for j-1 of the balloon phase rounds
                        savings += (j - 1) * 1.0 * depthCounts[j];

                        //These pebbles can be removed (curDepth - j) rounds before the baloon phase finishes
                        savings2 += (curDepth - j) * 1.0 * lastDepthCount[j];
                        //Cost to place re-place pebble on node during balloon phase
                        attackCost += depthCounts[j] * R;
                    }

                    attackCost += curDepth * 1.0 * i - savings - Math.Max(0, savings2 - curDepth * g);
                }
            }


            return attackCost;
        }


        /// <summary>
        /// Upper Bounds the cumulative cost of running the Alwen-Blocki Parallel Pebbling Attack on G given depth-reducing set S
        /// such that depth(G-S) <= d and |S| <= e.
        /// </summary>
        /// <param name="n">number of nodes in graph</param>
        /// <param name="DAG">input graph</param>
        /// <param name="S">depth reducing set</param>
        /// <param name="e">size of depth-reducing set |S|</param>
        /// <param name="d">upper bound on depth(G-S)</param>
        /// <returns>double: Upper bound on CC(G) </returns>
        public static double AttackCostUpperBound(int n, int[] DAG, bool[] S, int e, int d)
        {
            // Stores the depth of every node in DAG-S
            int[] depths = new int[n];
            double attackCost = 0.0;
            //Upper bound cost on parents of nodes in light phases + cost on S
            for (int i = 0; i < n; i++)
            {
                // Add Cost of pebbles on S
                if (S[i])
                {
                    // We keep a pebble on node i+1 in rounds i+1,...,n
                    attackCost += 1.0 * n - 1.0 * i;
                }
                // We start a balloon phase with <= e additional (not in S) pebbles on parents of next e nodes 
                // during a light phase we discard these pebbles as soon as we are done to during round j< g of a light phase
                // we have at most g-j pebbles on parents. The average number of pebbles per round is (g+(g-1)+...+1)/g = (g+1)/2.
                attackCost += e * 0.5 + 0.5;
            }
            // 
            int[] depthCounts = new int[n];
            int[] maxDepth = new int[n];

            for (int i = 1; i < DAG.Length; i++)
            {
                depthCounts[i] = 0;
                maxDepth[i] = 0;
            }
            depths[0] = (S[0]) ? 0 : 1;
            maxDepth[0] = depths[0];
            depthCounts[depths[0]] += 1;
            //Upper Bound costs of balloon phases
            for (int i = 1; i < DAG.Length; i++)
            {
                depths[i] = S[i] ? 0 : 1;
                if (S[i - 1] == false && S[i] == false) depths[i] = Math.Max(depths[i], depths[i - 1] + 1);
                if (S[DAG[i]] == false && S[i] == false) depths[i] = Math.Max(depths[i], depths[DAG[i]] + 1);
                depthCounts[depths[i]] += 1;

                maxDepth[i] = Math.Max(maxDepth[i - 1], depths[i]);
                if (0 == (i % e))
                {
                    int curDepth = maxDepth[i];
                    // Cost balloon phase <= curDepth*i - sum_{j=1}^{curDepth} \sum_{k =j+1}^{curdepth} depthCounts[k]
                    double savings = 0.0;
                    for (int j = 1; j <= curDepth; j++)
                    {
                        // We don't have pebbles on these nodes for j-1 of the balloon phase rounds
                        savings += (j - 1) * 1.0 * depthCounts[j];
                    }

                    attackCost += curDepth * 1.0 * i - savings;
                }
            }


            return attackCost;
        }

        /// <summary>
        /// Calculates energy cost of Greedy Pebbling Strategy where we 
        /// opportunistically discard a pebble as soon as we place a pebble on its 
        /// greatest child in the graph.
        /// </summary>
        /// <param name="n">number of nodes</param>
        /// <param name="DAG">array specifying edges in DAG</param>
        /// <param name="R">ratio of size of hash function on chip to size of RAM memory to store single output</param>
        /// <returns></returns>
        public static double GreedyPebblingAttack(int n, int[] DAG, double R)
        {
            int[] maxChild = new int[n];
            for (int i = 0; i <= n - 2; i++)
                maxChild[i] = i + 1;
            // special case: we remove pebble on last node at step n+1
            maxChild[n-1] = n+1;
            
            for (int i=1; i < n; i++)
            {
                int parent = DAG[i];
                maxChild[parent] = Math.Max(maxChild[parent], i);
            }


            // We never recompute/replace pebbles so this is cost of calling hash function n times
            double cost = R * n;

            // now maxChild[i] = max { v : (i,v) in E}
            // Greedy Pebble Attack 
            // 1) Places pebble on node i at time i
            //  2) Removes pebbles on node i at time maxChild[i]
            for (int i=0; i < n; i++)
            {
                // add total cost to keep pebble on node i
                cost+= maxChild[i] - i;
            }                
            return cost;
        }

        /// <summary>
        /// Calculates energy cost of Greedy Pebbling Strategy where we 
        /// opportunistically discard a pebble as soon as we place a pebble on its 
        /// greatest child in the graph.
        /// </summary>
        /// <param name="n">number of nodes</param>
        /// <param name="DAG">array specifying edges in DAG</param>
        /// <param name="R">ratio of size of hash function on chip to size of RAM memory to store single output</param>
        /// <returns></returns>
        public static double GreedyPebblingAttackOnTwoPassDAG(int n, int[] DAG, double R)
        {
            int[] maxChild = new int[n];
            for (int i = 0; i <= n - 2; i++)
                maxChild[i] = i + 1;
            // special case: we remove pebble on last node at step n+1
            maxChild[n - 1] = n + 1;

            for (int i = 1; i < n/2; i++)
            {
                int parent = DAG[i];
                maxChild[parent] = Math.Max(maxChild[parent], i);
            }
            // We never recompute/replace pebbles so this is cost of calling hash function n times
            // We pebble nodes 1...n/2 twice
            double cost = R * (n * 1.5);
            // aCcounts for memory costs for the second time we reppeble nodes 1...n/2
            for (int i = 0; i < n / 2; i++)
            {
                // add total cost to keep pebble on node i
                cost += maxChild[i] - i;
            }
            for (int i = n/2; i < n; i++)
            {
                int parent = DAG[i];
                maxChild[parent] = Math.Max(maxChild[parent], i);
            }


            // now maxChild[i] = max { v : (i,v) in E}
            // Greedy Pebble Attack 
            // 1) Places pebble on node i at time i
            //  2) Removes pebbles on node i at time maxChild[i]


            for (int i = 0; i < n; i++)
            {
                // add total cost to keep pebble on node i
                cost += maxChild[i] - i;
            }

            return cost;
        }

        /// <summary>
        /// One Round of Valiant Lemma
        /// </summary>
        /// <param name="DAG">The Input DAG</param>
        /// <returns>depth of new graph G-S and depth reducing set S</returns>
        public static Tuple<int, int, int, bool[]> BuildDepthReducingSetValiantLazy(int[] DAG, int targetDepth, int bbase = 2)
        {
            int[] depths = new int[DAG.Length];
            byte[] edgeLabels = new byte[DAG.Length];
            int[] EiCounts = new int[32];
            bool[] removed = new bool[32];
            int e = 0;
            for (int i = 0; i < 32; i++) removed[i] = false;
            bool[] S = new bool[DAG.Length];

            int currentDepth = DAG.Length;


            for (int i = 0; i < DAG.Length - 1; i++)
            {
                depths[i] = i + 1;
                EiCounts[MSDB(i, i + 1, bbase)]++;
                EiCounts[MSDB(i, DAG[i], bbase)]++;
            }
            EiCounts[MSDB(DAG.Length - 1, DAG[DAG.Length - 1], bbase)]++;



            while (currentDepth > targetDepth)
            {
                int j = -1;
                int min = 10 * DAG.Length;
                for (int i = 0; i < 30; i++)
                {
                    if (EiCounts[i] < min && EiCounts[i] > 0 && !removed[i])
                    {
                        min = EiCounts[i];
                        j = i;

                    }
                }

                removed[j] = true;

                int maxDepth = 0;
                int count = 0;
                for (int i = 1; i < DAG.Length; i++)
                {


                    if (MSDB(i - 1, i, bbase) == j || MSDB(i, DAG[i], bbase) == j)
                    {

                        S[i] = true;
                        depths[i] = 0;
                    }

                    if (S[i]) count++;
                }
                currentDepth = Depth(DAG, S); ;
                e = count;
            }


            return new Tuple<int, int, int, bool[]>(currentDepth, e, e, S);
        }

        /// <summary>
        /// Awsome attack implmentation
        /// </summary>
        /// <param name="DAG">cool graph</param>
        /// <param name="rounds"></param>
        /// <param name="bbase"></param>
        /// <returns>devastating results</returns>
        public static Tuple<int, int, int, bool[]> BuildDepthReducingSetValiant(int[] DAG, int rounds, int bbase = 2)
        {
            bool[] S = new bool[DAG.Length];
            for (int i = 0; i < DAG.Length; i++) S[i] = false;
            double honestCost = DAG.Length * 1.0 * DAG.Length / 2;
            double honestECost = DAG.Length * 3000.0 + DAG.Length * 1.0 * DAG.Length / 2;
            Tuple<int, int, int, bool[]> T = null;
            for (int i = 0; i < rounds; i++)
            {
                 

                T = BuildDepthReducingSetValiantHelper(DAG, S, bbase);
                S = T.Item4;
                int e = T.Item2;
                int d = T.Item1;
                int[] gs = new int[] { e / 4, e / 4 + 1 * e / 4, e / 4 + 2 * e / 4, e / 4 + 3 * e / 4, e / 4 + 4 * e / 4, e / 4 + 5 * e / 4, e / 4 + 6 * e / 4, e / 4 + 7 * e / 4, e / 4 + 8 * e / 4, e / 4 + 9 * e / 4, e / 4 + 10 * e / 4, e / 4 + 11 * e / 4, e / 4 + 12 * e / 4, e / 4 + 13 * e / 4, e / 4 + 14 * e / 4, e / 4 + 15 * e / 4 };

                Console.Out.WriteLine("----------------round = " + i + "------------");
                Console.Out.WriteLine("(d,d/n) = " + "(" + d + "," + d * 1.0 / (1.0 * DAG.Length) + ")");
                double attackCost = EAttackCostUpperBound(DAG.Length, DAG, S, e, d, gs, 3000.0);
                double attackCost1 = AttackCostUpperBound(DAG.Length, DAG, S, e, d);
                Console.Out.WriteLine("(e,e/n) = " + "(" + e + "," + e * 1.0 / (1.0 * DAG.Length) + ")");
                Console.WriteLine("attack quality= " + honestCost / (e * 1.0 * DAG.Length * 2.0 + DAG.Length * 1.0 * Math.Sqrt(DAG.Length * 1.0 * d * 1.0)));
                Console.WriteLine("attack quality= " + honestCost / attackCost1);
                Console.WriteLine("Energy attack quality= " + honestECost / attackCost);
            }
            return T;
        }

        /// <summary>
        /// Constructs depth-reducing sets S for G by iteratively applying Valiant's Lemma Attack
        /// </summary>
        /// <param name="DAG">Input DAG G</param>
        /// <param name="targetDepth">Goal reduce depth(G-S) <= d_{tgt} </param>
        /// <returns>Tuple.First component: depth of new graph G-S. Second/Third Component: e. Fourth Component::depth reducing set S</returns>
        public static Tuple<int, int, int, bool[]> DepthReducingSetValiant(int[] DAG, int targetDepth, int bbase = 2)
        {
            bool[] S = new bool[DAG.Length];

            Tuple<int, int, int, bool[]> T = null;

            int depth = DAG.Length;
            while (depth > targetDepth)
            {

                T = BuildDepthReducingSetValiantHelper(DAG, S, bbase);
                S = T.Item4;
                depth = T.Item1;
            }
            return T;

      
        }
        /// <summary>
        /// Constructs depth-reducing sets S for G in several different algorithms (Excluding the most computationally expensive approaches such as the greedy path cover algorithm). Returns the smallest set.
        /// </summary>
        /// <param name="DAG">Input DAG G</param>
        /// <param name="targetDepth">Goal reduce depth(G-S) <= d_{tgt} </param>
        /// <returns>Tuple.First component: depth of new graph G-S. Second/Third Component: e. Fourth Component::depth reducing set S</returns>
        public static Tuple<int, int, int, bool[]> BestDepthReducingSetExcludeGreedy(int[] DAG, int targetDepth)
        {
            Tuple<int, int, int, bool[]> T0 = Helpers.DepthReducingSetLayered(DAG, targetDepth , (int)Math.Floor(Math.Sqrt(targetDepth)));
            Tuple<int, int, int, bool[]> T1 = Helpers.DepthReducingSetHybrid(DAG, targetDepth, (int)Math.Floor(Math.Sqrt(targetDepth)));
            Tuple<int, int, int, bool[]> T2 = Helpers.DepthReducingSetValiant(DAG, targetDepth);
            Tuple<int, int, int, bool[]> T3 = Helpers.BuildDepthReducingSetValiantLazy(DAG, targetDepth);
            Tuple<int, int, int, bool[]> T4 = Helpers.BuildDepthReducingSetValiantLazy(DAG, targetDepth,3);
            int e0 = T0.Item2;
            int e1 = T1.Item2;
            int e2 = T2.Item2;
            int e3 = T3.Item2;
            int e4 = T4.Item2;
            if (e0 <= e1 && e0 <= e2 && e0 <= e3 && e0 <= e4) return T0;
            else if (e1 <= e2 && e1 <= e3 && e1 <= e4) return T1;
            else if (e2 <= e3 && e2 <= e4) return T2;
            else if (e3 <= e4) return T3;
            else return T4;
        }

        /// <summary>
        /// Outputs the most significant digit on which integers a and b (viewed in binary) differ.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="bbase">when b > 2 we can find the most significant differing digit when a and b viewed base bbase</param>
        /// <returns>index of the most significant differing digit (usually bit)</returns>
        static int MSDB(int a, int b, int bbase = 2)
        {
            int c = a ^ b;
            int r = -1;
            while (c > 0)
            {
                int rem = c % bbase;
                c = (c - rem) / bbase;
                r++;
            }
            if (r == -1) return 31;
            return r;
            // if (c == 0) return 31;
            //int pos = 0;
            //return  (int)(Math.Log(c, bbase));
            //return pos;
        }

        /// <summary>
        /// One Round of Valiant Lemma
        /// </summary>
        /// <param name="DAG">The Input DAG</param>
        /// <returns>depth of new graph G-S and depth reducing set S</returns>
        public static Tuple<int, int, int, bool[]> BuildDepthReducingSetValiantHelper(int[] DAG, bool[] currentS, int bbase = 2)
        {
            int[] depths = new int[DAG.Length];
            byte[] edgeLabels = new byte[DAG.Length];
            int[] EiCounts = new int[32];
            bool[] S = new bool[DAG.Length];
            depths[0] = (currentS[0]) ? 0 : 1;
            int currentDepth = 0;
            for (int i = 1; i < DAG.Length; i++)
            {
                depths[i] = currentS[i] ? 0:1;
                if (currentS[i - 1] == false) depths[i] = Math.Max(depths[i], depths[i - 1] + 1);
                if (currentS[DAG[i]] == false) depths[i] = Math.Max(depths[i], depths[DAG[i]] + 1);
                currentDepth = Math.Max(currentDepth, depths[i]);

            }
            for (int i = 0; i < DAG.Length - 1; i++)
            {
                if (!currentS[i] && !currentS[i + 1]) EiCounts[MSDB(depths[i], depths[i + 1], bbase)]++;
                if (!currentS[i] && !currentS[DAG[i]]) EiCounts[MSDB(depths[i], depths[DAG[i]], bbase)]++;
            }
            if (!currentS[DAG.Length - 1] && !currentS[DAG[DAG.Length - 1]]) EiCounts[MSDB(depths[DAG.Length - 1], depths[DAG[DAG.Length - 1]], bbase)]++;
            int j = -1;
            int min = 10 * DAG.Length;
            for (int i = 0; i < 30; i++)
            {
                if (EiCounts[i] < min && EiCounts[i] > 0)
                {
                    min = EiCounts[i];
                    j = i;
                }
            }
            int maxDepth = 0;
            int count = 0;
            for (int i = 1; i < DAG.Length; i++)
            {
                S[i] = currentS[i];

                if (MSDB(depths[i - 1], depths[i], bbase) == j || MSDB(depths[i], depths[DAG[i]], bbase) == j)
                {

                    S[i] = true;
                    depths[i] = 0;
                }

                if (S[i]) count++;
            }

            maxDepth = Depth(DAG, S);

            return new Tuple<int, int, int, bool[]>(maxDepth, count, count, S);
        }
        /// <summary>
        /// One Round of Valiant Lemma, output S which reduces depth to depth(G-S)<= 2^i whenver depth(G) <= 2^{i+1} 
        /// </summary>
        /// <param name="DAG">The Input DAG</param>
        /// <returns>Tuple.First component: depth of new graph G-S. Second/Third Component: e. Fourth Component::depth reducing set S</returns>
        public static Tuple<int, int, int, bool[]> BuildDepthReducingSetValiantHeuristic(int[] DAG, int bbase = 2)
        {
            int[] depths = new int[DAG.Length];
            byte[] edgeLabels = new byte[DAG.Length];
            int[] EiCounts = new int[32];
            bool[] S = new bool[DAG.Length];
            for (int i = 0; i < DAG.Length; i++)
            {

                depths[i] = i + 1;
            }
            for (int i = 0; i < DAG.Length - 1; i++)
            {
                EiCounts[MSDB(i, i + 1, bbase)]++;
                EiCounts[MSDB(i, DAG[i], bbase)]++;
            }
            EiCounts[MSDB(DAG.Length - 1, DAG[DAG.Length - 1], bbase)]++;
            int j = -1;
            int min = 10 * DAG.Length;
            for (int i = 0; i < 30; i++)
            {
                if (EiCounts[i] < min && EiCounts[i] > 0)
                {
                    min = EiCounts[i];
                    j = i;
                }
            }
            int maxDepth = 0;
            int count = 0;
            depths[0] = 1;
            for (int i = 1; i < DAG.Length; i++)
            {

                if (MSDB(i - 1, i, bbase) == j || MSDB(i, DAG[i], bbase) == j)
                {
                    count++;
                    S[i] = true;
                    depths[i] = 0;
                }
                else
                {
                    depths[i] = 1;
                    if (S[i - 1] == false) depths[i] = Math.Max(depths[i], depths[i - 1] + 1);
                    if (S[DAG[i]] == false) depths[i] = Math.Max(depths[i], depths[DAG[i]] + 1);
                    maxDepth = Math.Max(depths[i], maxDepth);
                }
            }



            return new Tuple<int, int, int, bool[]>(maxDepth, count, count, S);
        }

        /// <summary>
        /// Layered Attack from AB16 to construct depth-reducing set S s.t depth(G-S) <= d_{tgt} 
        /// </summary>
        /// <param name="DAG">input DAG G</param>
        /// <param name="targetDepth">target depth</param>
        /// <param name="numLayers">parameter of the attack</param>
        /// <returns> Tuple.First component: depth of new graph G-S. Second/Third Component: e. Fourth Component::depth reducing set S</returns>
        public static Tuple<int, int, int, bool[]> DepthReducingSetLayered(int[] DAG, int targetDepth, int numLayers)
        {
            int e = 0;
            bool[] S = new bool[DAG.Length];
            int gap = (int)Math.Floor(targetDepth * 1.0 / (1.0 * numLayers));
            int layerSize = (int)Math.Ceiling(DAG.Length / (1.0 * numLayers));

            for (int i = gap - 1; i < DAG.Length; i += gap)
            {
                S[i] = true;
                e++;
            }

            for (int j = 0; j < numLayers; j++)
            {
                int sizeOfCurrentLayer = (j < numLayers - 1) ? layerSize : DAG.Length - (layerSize * (numLayers - 1));
                for (int i = 0; i < sizeOfCurrentLayer; i++)
                {
                    int k = j * layerSize + i;
                    int placeInGap = (k + 1) % gap;
                    int parent = DAG[k];
                    int parentPlaceInGap = (k + 1) % gap;
                    int parentLayer = (int)Math.Floor(parent * 1.0 / (1.0 * layerSize));
                    if (parentLayer == j && parentPlaceInGap >= placeInGap && !S[k])
                    {

                        e++;
                        S[k] = true;

                    }

                }
            }
            int depth = Depth(DAG, S);

            Tuple<int, int, int, bool[]> T = new Tuple<int, int, int, bool[]>(depth, e, e, S);

            return T;


        }

        /// <summary>
        /// Builds a depth reducing set using the layered attack + CSP Solver.
        /// We divide the graph into layers and use CSP Solver to reduce the depth of each layer.
        /// </summary>
        /// <param name="DAG">Input Graph G</param>
        /// <param name="targetDepth">Goal to construct S s.t depth(G-S) <= d_{tgt}</param>
        /// <param name="numLayers">Parameter which specifies number of layers. Smaller numLayers will be computationally less efficient, but may result in better sets.</param>
        /// <returns>Tuple. First component is depth(G-S), Second/third component are e=|S| and final component is the depth reducing set S</returns>
        public static Tuple<int, int, int, bool[]> DepthReducingSetHybrid(int[] DAG, int targetDepth, int numLayers)
        {
            int e = 0;
            bool[] S = new bool[DAG.Length];
            int gap = (int)Math.Floor(targetDepth * 1.0 / (1.0 * numLayers));
            int layerSize = (int)Math.Ceiling(DAG.Length / (1.0 * numLayers));



            for (int j = 0; j < numLayers; j++)
            {
                int sizeOfCurrentLayer = (j < numLayers - 1) ? layerSize : DAG.Length - (layerSize * (numLayers - 1));
                int[] layerDag = new int[sizeOfCurrentLayer];

                layerDag[0] = 0;
                for (int i = 1; i < sizeOfCurrentLayer; i++)
                {

                    int k = j * layerSize + i;
                    int parent = DAG[k] - j * layerSize;
                    layerDag[i] = parent >= 0 ? parent : i - 1;

                }
                bool[] layerS = Helpers.DepthReducingSetValiant(layerDag, gap).Item4;
                for (int i = 0; i < sizeOfCurrentLayer; i++)
                {

                    int k = j * layerSize + i;
                    S[k] = layerS[i];
                    e++;
                }
            }

            for (int i = gap - 1; i < DAG.Length; i += gap)
            {
                if (!S[i])
                {
                    S[i] = true;
                    e++;
                }
            }
            int depth = Depth(DAG, S);

            Tuple<int, int, int, bool[]> T = new Tuple<int, int, int, bool[]>(depth, e, e, S);

            return T;


        }



        /// <summary>
        /// Outputs a DOT file to visualize a DAG
        /// </summary>
        /// <param name="DAG">A DAG with N = 2^n nodes</param>
        /// <param name="graphname">Name of the graph in DOT file</param>
        /// <returns>DOT file which can be rendered to visualize the DAG</returns>
        public static string PrintGraph(int[] DAG,string graphname)
        {
            string graph = "digraph " + graphname + " {\n rankdir=LR;\n";

            for(int i = 1; i < DAG.Length; i++)
            {
                
                string color = "black";

                
                graph += i - 1 + " -> " + i + "[penwidth=3,color=" + color + "];\n";
                int lengthColor = (int)Math.Log(i - DAG[i], 2.0);
                if (lengthColor < 1) color = "black";
                else if (lengthColor < 2) color = "purple";
                else if (lengthColor < 3) color = "violet";
                else if (lengthColor < 4) color = "indigo";
                else if (lengthColor < 5) color = "blue";
                else if (lengthColor < 6) color = "green";
                else if (lengthColor < 7) color = "yellow";
                else if (lengthColor < 8) color = "orange";
                else color = "red";
                if (DAG[i] < i - 1)  graph += DAG[i] + " -> " + i + "[penwidth=12,color=" + color + "];\n";
            }
            graph += "\n}";
            return graph;
        }
        /// <summary>
        /// Outputs a DOT file to visualize a DAG
        /// </summary>
        /// <param name="DAG">A DAG with N = 2^n nodes</param>
        /// <param name="graphname">Name of the graph in DOT file</param>
        /// <returns>DOT file which can be rendered to visualize the DAG</returns>
        public static string PrintAATGraph(int[] DAG, string graphname)
        {
            string graph = "digraph " + graphname + " {\n \n";
            graph += "subgraph cluster0 {\n fontsize =50; \n  rankdir=LR;\nstyle=filled; \n color=lightgrey;\n node [style=filled,color=white];\n";
            for (int i = 1; i < DAG.Length/2; i++)
            {

                string color = "black";


                graph += i - 1 + " -> " + i + "[penwidth=3,color=" + color + "];\n";
                int lengthColor = (int)Math.Log(i - DAG[i], 2.0);
                if (lengthColor < 1) color = "black";
                else if (lengthColor < 2) color = "purple";
                else if (lengthColor < 3) color = "violet";
                else if (lengthColor < 4) color = "indigo";
                else if (lengthColor < 5) color = "blue";
                else if (lengthColor < 6) color = "green";
                else if (lengthColor < 7) color = "yellow";
                else if (lengthColor < 8) color = "orange";
                else color = "red";
                if (DAG[i] < i - 1) graph += DAG[i] + " -> " + i + "[penwidth=12,color=" + color + "];\n";
            }
            graph += "\n label=\"DRSample\"\n}\n ";
            graph += "subgraph cluster1 {\n fontsize =50; \n rankdir=LR; \n node [style=filled];\n";
            for (int i = DAG.Length / 2+1; i < DAG.Length ; i++)
            {
                string color = "black";


                graph += i - 1 + " -> " + i + "[penwidth=3,color=" + color + "];\n";
                int lengthColor = (int)Math.Log(i - DAG[i], 2.0);
                if (lengthColor < 1) color = "black";
                else if (lengthColor < 2) color = "purple";
                else if (lengthColor < 3) color = "violet";
                else if (lengthColor < 4) color = "indigo";
                else if (lengthColor < 5) color = "blue";
                else if (lengthColor < 6) color = "green";
                else if (lengthColor < 7) color = "yellow";
                else if (lengthColor < 8) color = "orange";
                else color = "red";
                if (DAG[i] < i - 1 && DAG[i] >= DAG.Length/2) graph += DAG[i] + " -> " + i + "[penwidth=12,color=" + color + "];\n";
            }
            graph+= "color = blue; label=\"aAT Boost Gadget\"\n}\n";
            graph += DAG.Length/2-1 + " -> " + DAG.Length/2 + "[penwidth=3,color=black" +"];\n";

            for (int i =  1; i < DAG.Length; i++)
            {
                string color = "black";


                int lengthColor = (int)Math.Log(i - DAG[i], 2.0);
                if (i >= DAG.Length / 2 && DAG[i] < DAG.Length / 2)
                {
                    if (lengthColor < 1) color = "black";
                    else if (lengthColor < 2) color = "purple";
                    else if (lengthColor < 3) color = "violet";
                    else if (lengthColor < 4) color = "indigo";
                    else if (lengthColor < 5) color = "blue";
                    else if (lengthColor < 6) color = "green";
                    else if (lengthColor < 7) color = "yellow";
                    else if (lengthColor < 8) color = "orange";
                    else color = "red";
                    if (DAG[i] < i - 1) graph += DAG[i] + " -> " + i + "[penwidth=12,color=" + color + "];\n";
                }
            }

            graph += "\n}";
            return graph;
        }

        public static double LowDepthPebble(int [] DAG, bool [] S, double R=3000.0)
        {
            int depth = Depth(DAG, S);
            int[] pebbles = new int[DAG.Length];
            int[,] ancestorsOfVAtDepth = new int[DAG.Length, depth+1];

            ancestorsOfVAtDepth[0, 0] = (S[0] ? 0 : 1); for (int d = 1; d <= depth; d++) ancestorsOfVAtDepth[0, d] = 0;
            for (int v=1; v < DAG.Length; v++)
            {
                ancestorsOfVAtDepth[v, 0] = (S[v]?0: 1);
                //if (S[v]) continue;
                for (int d=1;d<= depth;d++)
                {
                    if (S[v]) ancestorsOfVAtDepth[v, d] = 0;
                    else
                    {
                        ancestorsOfVAtDepth[v, d] += ancestorsOfVAtDepth[v - 1, d - 1] + ancestorsOfVAtDepth[DAG[v], d - 1];
                    }
                }
            }


            double cost = 0;
            int numPebbles = 0;
            cost += R + 1;
            for (int v = 1; v < DAG.Length; v++)
            {
                for(int i=v; i < Math.Min(DAG.Length,v+depth); i++)
                {
                    numPebbles+=ancestorsOfVAtDepth[DAG[v], i-v];
                    cost+= (R)*ancestorsOfVAtDepth[DAG[v], i - v];
                }
                for (int i= Math.Max(0,v-1); i<= v-1; i++)
                {
                    numPebbles -= ancestorsOfVAtDepth[DAG[i], depth - (v - i)];
                }
                cost += numPebbles;
            }
            return cost;

        }



        /// <summary>
        /// Makes a deep copy of input array
        /// </summary>
        /// <param name="input"> input array</param>
        /// <returns>copied array</returns>
        public static BigInteger[] copyArray (BigInteger[] input)
        {
            BigInteger [] copy = new BigInteger[input.Length];
            for (int i = 0; i < input.Length; i++) copy[i] = input[i];
            return copy;
        }
        /// <summary>
        /// Makes a deep copy of input array
        /// </summary>
        /// <param name="input"> input array</param>
        /// <returns>copied array</returns>
        public static double[] copyArray(double[] input)
        {
            double[] copy = new double[input.Length];
            for (int i = 0; i < input.Length; i++) copy[i] = input[i];
            return copy;
        }
        

      

        /// <summary>
        /// Counts the number of length d paths in G-S, where G is a two pass DAG. This version only gives an approximation since we use the double type for dynamic programming.
        /// </summary>
        /// <param name="DAG">DAG G (two pass DAG includes edges (i,i+n/2))</param>
        /// <param name="d">path length</param>
        /// <param name="S">set of deleted nodes S</param>
        /// <returns>Tuple. Item 1 = Number of length d paths in G-s, Item 2 = max_v #depth d paths incident to v. Item 3= index of vertex v incident to most length d paths</returns>
        public static Tuple<double, double, int> CountPathsApxTwoPass(int[] DAG, int d, bool[] S)
        {
            int lengthFirstPass = DAG.Length / 2;
            int biggestDistanceInStack = 0;
            Stack<Tuple<double[], int>> DistDPathCountEndAtV = new Stack<Tuple<double[], int>>();
            double[] tailArray = new double[d + 1];
            double[] copyEndingAtNodeVOfLengthZero = new double[DAG.Length];
            double[] pathsEndingAtNodeVofLengthi = new double[DAG.Length];
            double[] pathsEndingAtNodeVofLengthiMinus1 = new double[DAG.Length];
            double[] pathsStartingAtNodeVofLengthi = new double[DAG.Length];
            double[] pathsStartingAtNodeVofLengthiMinus1 = new double[DAG.Length];
            double totalPathsOfLengthd = 0, maxVNumPathsOfLengthdIncidentToV = 0;
            double[] pathsIncidentToV = new double[DAG.Length];

            tailArray[0] = (DAG.Length >= d) ? 1 : 0;
            for (int j = 0; j < DAG.Length; j++)
            {
                copyEndingAtNodeVOfLengthZero[j] = (S[j] == false) ? 1 : 0;
                pathsEndingAtNodeVofLengthiMinus1[j] = (S[j] == false) ? 1 : 0;
                pathsStartingAtNodeVofLengthiMinus1[j] = (S[j] == false) ? 1 : 0;
            }
            DistDPathCountEndAtV.Push(new Tuple<double[], int>(copyArray(pathsEndingAtNodeVofLengthiMinus1), 0));
            for (int i = 1; i <= d; i++)
            {
                if ((i % 1000) == 0) Console.Write(".");
                pathsEndingAtNodeVofLengthi[0] = 0;
                for (int v = 1; v < DAG.Length; v++)
                {

                    if (S[v]) pathsEndingAtNodeVofLengthi[v] = 0;
                    else pathsEndingAtNodeVofLengthi[v] = 
                            (S[v - 1] ? 0 : pathsEndingAtNodeVofLengthiMinus1[v - 1])
                            + ((DAG[v] < v - 1 && S[DAG[v]] == false) ? pathsEndingAtNodeVofLengthiMinus1[DAG[v]] : 0)
                            + ((v >= lengthFirstPass && S[v-lengthFirstPass]==false)? pathsEndingAtNodeVofLengthiMinus1[v - lengthFirstPass] : 0) ;

                }
                //if ((i % 20) == 4) Console.WriteLine("Num Paths of Length" + i + "ending at  node" + (DAG.Length / 2-1) + " is: " + pathsEndingAtNodeVofLengthi[DAG.Length / 2  - 1]);
                if (i > 1 + d - (d - biggestDistanceInStack) / 2)
                {

                    DistDPathCountEndAtV.Push(new Tuple<double[], int>(copyArray(pathsEndingAtNodeVofLengthi), i));
                    biggestDistanceInStack = i;
                }
                for (int v = 0; v < DAG.Length; v++)
                {
                    pathsEndingAtNodeVofLengthiMinus1[v] = pathsEndingAtNodeVofLengthi[v];
                }
                tailArray[i] = pathsEndingAtNodeVofLengthi[DAG.Length - 1];
            }

            for (int v = 0; v < DAG.Length; v++)
            {
                pathsIncidentToV[v] += (S[v]) ? 0 : pathsEndingAtNodeVofLengthi[v];
                totalPathsOfLengthd += (S[v]) ? 0 : pathsEndingAtNodeVofLengthi[v];
            }
            //Console.WriteLine("\nHalfway There");
            double totalPathsOfLengthd2 = 0, totalPathsOfLengthd3 = 0;
            for (int i = 1; i <= d; i++)
            {
                if ((i % 1000) == 0) Console.Write(".");
                pathsStartingAtNodeVofLengthi[DAG.Length - 1] = 0;
                for (int v = DAG.Length - 1; v >= 0; v--)
                {
                    if (DAG.Length - v < i || S[v]) pathsStartingAtNodeVofLengthi[v] = 0;
                    else if (DAG.Length - v == i)
                    {
                        pathsStartingAtNodeVofLengthi[v] = 0;
                        if (DAG[v] < v - 1 && S[DAG[v]] == false) pathsStartingAtNodeVofLengthi[DAG[v]] += pathsStartingAtNodeVofLengthiMinus1[v];
                        if (v >= lengthFirstPass && S[v - lengthFirstPass] == false) pathsStartingAtNodeVofLengthi[v - lengthFirstPass] += pathsStartingAtNodeVofLengthiMinus1[v];
       
                    }
                    else
                    {
                        pathsStartingAtNodeVofLengthi[v] += S[v + 1] ? 0 : pathsStartingAtNodeVofLengthiMinus1[v + 1];
                        if (DAG[v] < v - 1 && S[DAG[v]] == false) pathsStartingAtNodeVofLengthi[DAG[v]] += pathsStartingAtNodeVofLengthiMinus1[v];
                        if (v >= lengthFirstPass && S[v-lengthFirstPass] == false) pathsStartingAtNodeVofLengthi[v - lengthFirstPass] += pathsStartingAtNodeVofLengthiMinus1[v];
                    }



                }


                while (biggestDistanceInStack > d - i) { if (DistDPathCountEndAtV.Count > 1) { DistDPathCountEndAtV.Pop(); biggestDistanceInStack = DistDPathCountEndAtV.Peek().Item2; } else biggestDistanceInStack = 0; }
                Tuple<double[], int> T = DistDPathCountEndAtV.Peek();
                if (T.Item2 >= d - i) { DistDPathCountEndAtV.Pop(); if (DistDPathCountEndAtV.Count > 0) biggestDistanceInStack = DistDPathCountEndAtV.Peek().Item2; else biggestDistanceInStack = 0; }
                pathsEndingAtNodeVofLengthi = copyArray(T.Item1);
                pathsEndingAtNodeVofLengthiMinus1 = copyArray(T.Item1);
                for (int y = T.Item2; y < d - i; y++)
                {
                    pathsEndingAtNodeVofLengthi[0] = 0;

                    for (int v = 1; v < DAG.Length; v++)
                    {

                        if (S[v]) pathsEndingAtNodeVofLengthi[v] = 0;// = pathsEndingAtNodeVofLengthiMinus1[v - 1];
                        else pathsEndingAtNodeVofLengthi[v] = (S[v - 1] ? 0 : pathsEndingAtNodeVofLengthiMinus1[v - 1])
                                + ((DAG[v] < v - 1 && S[DAG[v]] == false) ? pathsEndingAtNodeVofLengthiMinus1[DAG[v]] : 0)
                                +((v >= lengthFirstPass && S[v - lengthFirstPass] == false) ? pathsEndingAtNodeVofLengthiMinus1[v - lengthFirstPass] : 0);

                    }
                    if (y + 1 >= (d - i) - ((d - i) - biggestDistanceInStack) / 2 && y < d - i)
                    {
                        DistDPathCountEndAtV.Push(new Tuple<double[], int>(copyArray(pathsEndingAtNodeVofLengthi), y + 1));
                        biggestDistanceInStack = y + 1;
                    }
                    for (int v = 0; v < DAG.Length && y + 1 < d - i; v++)
                    {
                        pathsEndingAtNodeVofLengthiMinus1[v] = pathsEndingAtNodeVofLengthi[v];
                    }
                }
                for (int v = 0; v < DAG.Length; v++)
                {
                    pathsIncidentToV[v] += pathsEndingAtNodeVofLengthi[v] * pathsStartingAtNodeVofLengthi[v];

                }


                //    if (((d-i)%20)==4)    Console.WriteLine("Num Paths of Length" + (d-i) + " ending at  node " + (DAG.Length/2 - 1) + " is: " + pathsEndingAtNodeVofLengthi[DAG.Length / 2 - 1]);

                for (int v = 0; v < DAG.Length; v++)
                {
                    pathsStartingAtNodeVofLengthiMinus1[v] = pathsStartingAtNodeVofLengthi[v];
                    pathsStartingAtNodeVofLengthi[v] = 0;
                    pathsEndingAtNodeVofLengthi[v] = pathsEndingAtNodeVofLengthiMinus1[v];
                }
            }
            int vMax = 0;
            for (int v = 0; v < DAG.Length; v++)
            {
                totalPathsOfLengthd3 += pathsIncidentToV[v];
                pathsStartingAtNodeVofLengthi[v] = pathsStartingAtNodeVofLengthiMinus1[v];

                totalPathsOfLengthd2 += pathsStartingAtNodeVofLengthi[v];

                if (maxVNumPathsOfLengthdIncidentToV < pathsIncidentToV[v])
                {
                    maxVNumPathsOfLengthdIncidentToV = pathsIncidentToV[v];
                    vMax = v;
                }
                maxVNumPathsOfLengthdIncidentToV = Math.Max(maxVNumPathsOfLengthdIncidentToV, pathsIncidentToV[v]);
            }
            totalPathsOfLengthd3 = totalPathsOfLengthd3 / (d + 1);

            return new Tuple<double, double, int>(totalPathsOfLengthd, maxVNumPathsOfLengthdIncidentToV, vMax);
        }

        /// <summary>
        /// Computes the detph of G-S for a single pass DAG G and a set S of deleted nodes
        /// </summary>
        /// <param name="DAG">Input DAG G</param>
        /// <param name="S">set of deleted nodes</param>
        /// <returns></returns>
        public static int Depth(int[] DAG, bool[] S)
        {

            int[] depths = new int[DAG.Length];
            depths[0] = S[0] ? 0 : 1;
            int maxDepth = depths[0];
            for (int i = 1; i < DAG.Length; i++)
            {
                depths[i] = 0;
                int parent = DAG[i];
                int parDepth = S[parent] ? 0 : depths[parent];
                int parDepth2 = S[i - 1] ? 0 : depths[i - 1];
                depths[i] =  Math.Max(1+parDepth,1+ parDepth2);
                if (S[i]) depths[i] = 0;
                maxDepth = Math.Max(maxDepth, depths[i]);
            }
            return maxDepth;
        }

        /// <summary>
        /// Computes the detph of G-S for a two pass DAG G and a set S of deleted nodes
        /// </summary>
        /// <param name="DAG">Input DAG G</param>
        /// <param name="S">set of deleted nodes</param>
        /// <returns></returns>
        public static int DepthTwoPassDAG(int[] DAG, bool[] S)
        {

            int[] depths = new int[DAG.Length];
            depths[0] = S[0] ? 0 : 1;
            int maxDepth = depths[0];
            for (int i = 1; i < DAG.Length; i++)
            {
                depths[i] = 0;
                int parent = DAG[i];
                int parDepth = S[parent] ? 0 : depths[parent];
                int parDepth2 = S[i - 1] ? 0 : depths[i - 1];
                int parDepth3 = (i >= DAG.Length / 2 && S[i - DAG.Length / 2] == false) ? depths[i-DAG.Length/2] : 0;
                depths[i] = Math.Max(1 + parDepth, Math.Max(1 + parDepth2, 1+ parDepth3));
                if (S[i]) depths[i] = 0;
                maxDepth = Math.Max(maxDepth, depths[i]);
            }
            return maxDepth;
        }
    }
}
