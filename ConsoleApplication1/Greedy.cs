using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ConsoleApplication1
{
    class Greedy
    {
        public static double[] copyArray(double[] incidentSet)
        {
            double[] copy = new double[incidentSet.Length];
            for (int i = 0; i < incidentSet.Length; i++) copy[i] = incidentSet[i];
            return copy;
        }

        // Heap helper
        class FirstCompTuple<T, U> : IComparable where T : IComparable
        {
            public T first;
            public U second;

            public FirstCompTuple(T newfirst = default(T), U newsecond = default(U))
            {
                first = newfirst;
                second = newsecond;
            }

            int IComparable.CompareTo(object obj)
            {
                FirstCompTuple<T, U> toCompare = (FirstCompTuple<T, U>)obj;
                return this.first.CompareTo(toCompare.first);
            }
        }

        public static int[] GetKLargestUnvisited(int[] arr, int k, bool[] S)
        {
            // Stores index/value pairs for easier retrieval
            Heap<FirstCompTuple<int, int>> maxK = new Heap<FirstCompTuple<int, int>>();
            int i = 0;
            FirstCompTuple<int, int>[] firstCompToRet;

            int added = 0;
            // Try to fill the heap with at least k elements.
            while (i < arr.Length && added < k)
            {
                if (!S[i])
                {
                    maxK.insert(new FirstCompTuple<int, int>(arr[i], i));
                    added++;
                }
                i++;
            }
            if (added < k) // Probably wont happen, but just in case
                firstCompToRet = maxK.toSortedList();
            else
            {

                for (i = k; i < arr.Length; i++)
                {
                    if (!S[i] && arr[i] > maxK.minValue().first)
                    {
                        maxK.removeMin();
                        maxK.insert(new FirstCompTuple<int, int>(arr[i], i));

                    }
                }
                firstCompToRet = maxK.toSortedList();
            }
            // Turn this into a list of the node labels with highest num paths
            int[] toRet = new int[firstCompToRet.Length];
            for (i = 0; i < toRet.Length; i++)
                toRet[i] = firstCompToRet[i].second;
            return toRet;
        }

        // Use the more efficient GetKLargestUnvisited instead
        public static int[] myPickTopNodesNew(double[] arr, int topIndex, bool[] S)
        {
            // 

            // Pick the biggest "topIndex" nodes in array arr
            double mymax;
            int minIndex;

            int[] visited = new int[arr.Length];
            for (int k = 0; k < arr.Length; k++)
                visited[k] = 0;

            int[] topNodes = new int[topIndex];
            for (int i = 0; i < topIndex; i++) // obtain the first topIndex number of nodes in array
            {
                mymax = double.NegativeInfinity;
                minIndex = -1;
                for (int j = 0; j < arr.Length; j++)
                {
                    if (!S[j])
                    {
                        if (arr[j] > mymax && visited[j] == 0)
                        {
                            mymax = arr[j];
                            minIndex = j;
                        }
                    }

                }
                if (minIndex != -1)
                {
                    visited[minIndex] = 1;
                    topNodes[i] = minIndex;
                    //Console.Write("^");
                }
            }
            return topNodes;
        }

        public static int[] myPickTopNodes(double[] arr, int topIndex, bool[] S)
        {
            // Pick the biggest "topIndex" nodes in array arr
            double mymax;
            int minIndex;

            int[] visited = new int[arr.Length];
            for (int k = 0; k < arr.Length; k++)
                visited[k] = 0;

            int[] topNodes = new int[topIndex];
            for (int i = 0; i < topIndex; i++) // obtain the first topIndex number of nodes in array
            {
                mymax = -1;
                minIndex = -1;
                for (int j = 0; j < arr.Length; j++)
                {
                    if (!S[j])
                    {
                        if (arr[j] > mymax && visited[j] == 0)
                        {
                            mymax = arr[j];
                            minIndex = j;
                        }
                    }

                }
                if (minIndex != -1)
                {
                    visited[minIndex] = 1;
                    topNodes[i] = minIndex;
                }
            }
            return topNodes;
        }

        public static Tuple<int, int, int, bool[]> myGreedyDRSetsApxNoSortingNew(int[] DAG, int d, int radius, int topIndex, int dTarget = 16)
        {
            // Pick top nodes using brute force
            bool[] S = new bool[DAG.Length]; // default boolean is false
            List<int>[] DAGrev = findDAGrev(DAG);
            int[] incidentSet;
            int e = 0;
            bool[] inRadius = new bool[DAG.Length];
            int ii;

            // Array that stores whether or not a node is within the specified radius of a removed node
            for (ii = 0; ii < DAG.Length; ii++)
                inRadius[ii] = false;

            Tuple<int, int[]> T = CountPathsApxLimitDInt(DAG, d, S, dTarget);
            incidentSet = T.Item2;

            int[] topNodes = GetKLargestUnvisited(incidentSet, topIndex, S);

            // Basic Approach: Delete nodes incident to most paths, mark nearby nodes as inelligible and repeat with elligible nodes
            while (Helpers.Depth(DAG, S) > d)
            {
                S[topNodes[0]] = true;
                NodesWithinRadiusOfV(DAG, DAGrev, topNodes[0], S, radius, ref inRadius);

                for (int i = 1; i < topNodes.Length; i++)
                {
                    if (!inRadius[topNodes[i]]) // node i is not in the circle of nodes before it
                    {
                        S[topNodes[i]] = true;
                        NodesWithinRadiusOfV(DAG, DAGrev, topNodes[i], S, radius, ref inRadius);
                    }
                }
                T = CountPathsApxLimitDInt(DAG, d, S, dTarget);
                incidentSet = T.Item2;
                //  Console.Write("*" + CustomDouble.getContents(T.Item1) + "*" + CustomDouble.getContents(incidentSet[22]) + "&"); 
                topNodes = GetKLargestUnvisited(incidentSet, topIndex, S);
                S[topNodes[0]] = true; // check
                for (ii = 0; ii < DAG.Length; ii++)
                {
                    inRadius[ii] = false;
                }
            }
            e = 0;
            for (ii = 0; ii < DAG.Length; ii++)
                if (S[ii]) e++;

            int depth = Helpers.Depth(DAG, S);
            Tuple<int, int, int, bool[]> Tt = new Tuple<int, int, int, bool[]>(depth, e, e, S);
            return Tt;
        }


        public static Tuple<double, double[]> CountPathsApxLimitD(int[] DAG, int dTarget, bool[] S, int dLimit = 16)
        {
            dLimit = Math.Min(dTarget, dLimit);
            double[,] pathsEndingAtNodeV = new double[dLimit + 1, DAG.Length];
            double[,] pathsStartingAtNodeV = new double[dLimit + 1, DAG.Length];
            double[] pathsIncidentToV = new double[DAG.Length];
            double totalPathsOfLengthD = 0;

            for (int j = 0; j < DAG.Length; j++)
            {
                pathsEndingAtNodeV[0, j] = (S[j] == false) ? 1.0 : 0.0;
                pathsStartingAtNodeV[0, j] = (S[j] == false) ? 1.0 : 0.0;
                pathsIncidentToV[j] = 0.0;
            }

            for (int d = 1; d <= dLimit; d++)
            {
                pathsEndingAtNodeV[d, 0] = 0;
                for (int v = 1; v < DAG.Length; v++)
                {
                    if (S[v])
                    {
                        pathsEndingAtNodeV[d, v] = 0;
                        pathsStartingAtNodeV[d, v] = 0;
                    }
                    else
                    {
                        pathsEndingAtNodeV[d, v] = pathsEndingAtNodeV[d - 1, v - 1] + pathsEndingAtNodeV[d - 1, DAG[v]];
                        if (S[v - 1] == false) pathsStartingAtNodeV[d, v - 1] += pathsStartingAtNodeV[d - 1, v];
                        if (S[DAG[v]] == false) pathsStartingAtNodeV[d, DAG[v]] += pathsStartingAtNodeV[d - 1, v];
                    }
                }
            }

            for (int v = 0; v < DAG.Length; v++)
            {
                for (int d = 0; d <= dLimit; d++)
                {
                    pathsIncidentToV[v] += pathsStartingAtNodeV[d, v] * pathsEndingAtNodeV[dLimit - d, v];
                }
                totalPathsOfLengthD += pathsEndingAtNodeV[dLimit, v];
            }

            return new Tuple<double, double[]>(totalPathsOfLengthD, pathsIncidentToV);

        }

        public static Tuple<int, int[]> CountPathsApxLimitDInt(int[] DAG, int dTarget, bool[] S, int dLimit = 16)
        {
            dLimit = Math.Min(dTarget, dLimit);
            int[,] pathsEndingAtNodeV = new int[dLimit + 1, DAG.Length];
            int[,] pathsStartingAtNodeV = new int[dLimit + 1, DAG.Length];
            int[] pathsIncidentToV = new int[DAG.Length];
            int totalPathsOfLengthD = 0;

            for (int j = 0; j < DAG.Length; j++)
            {
                pathsEndingAtNodeV[0, j] = (S[j] == false) ? 1 : 0;
                pathsStartingAtNodeV[0, j] = (S[j] == false) ? 1 : 0;
                pathsIncidentToV[j] = 0;
            }

            for (int d = 1; d <= dLimit; d++)
            {
                pathsEndingAtNodeV[d, 0] = 0;
                for (int v = 1; v < DAG.Length; v++)
                {
                    if (S[v])
                    {
                        pathsEndingAtNodeV[d, v] = 0;
                        pathsStartingAtNodeV[d, v] = 0;
                    }
                    else
                    {
                        pathsEndingAtNodeV[d, v] = pathsEndingAtNodeV[d - 1, v - 1] + pathsEndingAtNodeV[d - 1, DAG[v]];
                        if (S[v - 1] == false) pathsStartingAtNodeV[d, v - 1] += pathsStartingAtNodeV[d - 1, v];
                        if (S[DAG[v]] == false) pathsStartingAtNodeV[d, DAG[v]] += pathsStartingAtNodeV[d - 1, v];
                    }
                }
            }

            for (int v = 0; v < DAG.Length; v++)
            {
                for (int d = 0; d <= dLimit; d++)
                {
                    pathsIncidentToV[v] += pathsStartingAtNodeV[d, v] * pathsEndingAtNodeV[dLimit - d, v];
                }
                totalPathsOfLengthD += pathsEndingAtNodeV[dLimit, v];
            }

            return new Tuple<int, int[]>(totalPathsOfLengthD, pathsIncidentToV);

        }


        public static List<int>[] findDAGrev(int[] DAG)
        {
            // based the DAG build a reversed graph storing all the successors
            List<int>[] DAGrev = new List<int>[DAG.Length];
            // Initialize lists
            for (int i = 0; i < DAG.Length; i++)
            {
                DAGrev[i] = new List<int>();
            }

            for (int i = 0; i < DAG.Length; i++)
                DAGrev[DAG[i]].Add(i);
            return DAGrev;
        }

        // Fixed BFS version 
        // double check DRS checking
        // add in our keyword
        public static bool[] NodesWithinRadiusOfV(int[] DAG, List<int>[] DAGrev, int v, bool[] S, int radius, ref bool[] inRadius)
        {
            // using G and Grev to do BFS starting from node v
            HashSet<int> BFSSet = new HashSet<int>();
            HashSet<int> BFSNext;

            // First expansion
            if (v > 0 && !S[v - 1])
            {
                inRadius[v - 1] = true;
                BFSSet.Add(v - 1);
            }
            if (v < DAG.Length - 1 && !S[v + 1])
            {
                inRadius[v + 1] = true;
                BFSSet.Add(v + 1);
            }
            if (!S[DAG[v]])
            {
                inRadius[DAG[v]] = true;
                BFSSet.Add(DAG[v]);
            }
            foreach (int x in DAGrev[v])
            {
                if (!S[x])
                {
                    inRadius[x] = true;
                    BFSSet.Add(x);
                }
            }

            // Expand "radius" - 1 more times
            for (int i = 1; i < radius; i++)
            {
                // Process all elements in the expanding set
                BFSNext = new HashSet<int>();
                foreach (int vertex in BFSSet)
                {
                    // Still need to continue BFS even if node is inradius of another node we deleted
                    //if(inRadius[vertex]) continue;
                    //else 
                    if (vertex < v && vertex > 0) // looking backward
                    {
                        if (!S[vertex - 1])
                        {
                            BFSNext.Add(vertex - 1);
                            inRadius[vertex - 1] = true;
                        }
                        if (!S[DAG[vertex]])
                        {
                            BFSNext.Add(DAG[vertex]);
                            inRadius[DAG[vertex]] = true;
                        }
                    }
                    else if (vertex > v && vertex < DAG.Length - 1) // looking forward
                    {
                        if (!S[vertex + 1])
                        {
                            BFSNext.Add(vertex + 1);
                            inRadius[vertex + 1] = true;
                        }
                        foreach (int nextNode in DAGrev[vertex])
                            if (!S[nextNode])
                            {
                                BFSNext.Add(nextNode);
                                inRadius[nextNode] = true;
                            }
                    }
                    BFSSet = BFSNext;
                }
            }

            return inRadius;

        }

    }
}


