using System;
using System.Collections.Generic;
// using System.Linq;
// using System.Text;
// using System.Globalization;


namespace rewrite
{
	class radiusBFS
	{
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

        public static bool[] myGreedyDRSetsApxNoSorting(int[] DAG, int d, int radius, int topIndex)
        {
            // Pick top nodes using brute force
            bool[] S = new bool[DAG.Length]; // default boolean is false
            List<int>[] DAGadj = findDAGadj(DAG);
            double[] incidentSet;
            int[] inRadius = new int[DAG.Length];
            int ii;
            for (ii = 0; ii < DAG.Length; ii++)
                inRadius[ii] = 0;
            
            Tuple<double, double[]> T = myCountPathsApx(DAG, d, S);
            incidentSet = T.Item2;
            
            int[] topNodes = myPickTopNodes(incidentSet, topIndex, S);
           
            while (T.Item1 > 0)
            {
                S[topNodes[0]] = true;
                inRadius = BFSfromVDAGadj(DAG, DAGadj, topNodes[0], S, radius, inRadius);
                
                for (int i = 1; i < topNodes.Length; i++)
                {
                    if (inRadius[topNodes[i]] != 1) // node i is not in the circle of nodes before it
                    {
                        S[topNodes[i]] = true;
                        inRadius = BFSfromVDAGadj(DAG, DAGadj, topNodes[i], S, radius, inRadius);
                    }
                }   
                
                T = myCountPathsApx(DAG, d, S); 
                incidentSet = T.Item2;
                topNodes = myPickTopNodes(incidentSet, topIndex, S);
                S[topNodes[0]] = true;
                for (ii = 0; ii < DAG.Length; ii++)
                    inRadius[ii] = 0;
            }
            return S;
        }

        public static Tuple<double, double[]> myCountPathsApx(int[] DAG, int d, bool[] S)
        {
            // Basically the same as the original method except the return
            int biggestDistanceInStack = 0;
            Stack<Tuple<double[], int>> DistDPathCountEndAtV = new Stack<Tuple<double[], int>>();
            // double[] tailArray = new double[d + 1];
            double[] copyEndingAtNodeVOfLengthZero = new double[DAG.Length];
            double[] pathsEndingAtNodeVofLengthi = new double[DAG.Length];
            double[] pathsEndingAtNodeVofLengthiMinus1 = new double[DAG.Length];
            double[] pathsStartingAtNodeVofLengthi = new double[DAG.Length];
            double[] pathsStartingAtNodeVofLengthiMinus1 = new double[DAG.Length];
            double totalPathsOfLengthd = 0;//, maxVNumPathsOfLengthdIncidentToV = 0;
            double[] pathsIncidentToV = new double[DAG.Length];

            // tailArray[0] = (DAG.Length >= d) ? 1 : 0;
            for (int j = 0; j < DAG.Length; j++)
            {
                copyEndingAtNodeVOfLengthZero[j] = (S[j] == false) ? 1 : 0;
                pathsEndingAtNodeVofLengthiMinus1[j] = (S[j] == false) ? 1 : 0;
                pathsStartingAtNodeVofLengthiMinus1[j] = (S[j] == false) ? 1 : 0;
            }
            DistDPathCountEndAtV.Push(new Tuple<double[], int>(ConsoleApplication1.Helpers.copyArray(pathsEndingAtNodeVofLengthiMinus1), 0));
            for (int i = 1; i <= d; i++)
            {
                // if ((i % 1000) == 0) Console.Write(".");
                pathsEndingAtNodeVofLengthi[0] = 0;
                for (int v = 1; v < DAG.Length; v++)
                {

                    if (S[v]) pathsEndingAtNodeVofLengthi[v] = 0;
                    else pathsEndingAtNodeVofLengthi[v] = (S[v - 1] ? 0 : pathsEndingAtNodeVofLengthiMinus1[v - 1])
                            + ((DAG[v] < v - 1 && S[DAG[v]] == false) ? pathsEndingAtNodeVofLengthiMinus1[DAG[v]] : 0);

                }
                //if ((i % 20) == 4) Console.WriteLine("Num Paths of Length" + i + "ending at  node" + (DAG.Length / 2-1) + " is: " + pathsEndingAtNodeVofLengthi[DAG.Length / 2  - 1]);
                if (i > 1 + d - (d - biggestDistanceInStack) / 2)
                {

                    DistDPathCountEndAtV.Push(new Tuple<double[], int>(ConsoleApplication1.Helpers.copyArray(pathsEndingAtNodeVofLengthi), i));
                    biggestDistanceInStack = i;
                }
                for (int v = 0; v < DAG.Length; v++)
                {
                    pathsEndingAtNodeVofLengthiMinus1[v] = pathsEndingAtNodeVofLengthi[v];
                }
                // tailArray[i] = pathsEndingAtNodeVofLengthi[DAG.Length - 1];
            }

            for (int v = 0; v < DAG.Length; v++)
            {
                // @Lu pathsEndingAtNodeVofLengthi stores the len d paths
                pathsIncidentToV[v] += (S[v]) ? 0 : pathsEndingAtNodeVofLengthi[v];
                totalPathsOfLengthd += (S[v]) ? 0 : pathsEndingAtNodeVofLengthi[v];
            }
            //Console.WriteLine("\nHalfway There");
            // double totalPathsOfLengthd2 = 0, totalPathsOfLengthd3 = 0;
            for (int i = 1; i <= d; i++)
            {
                // if ((i % 1000) == 0) Console.Write(".");
                pathsStartingAtNodeVofLengthi[DAG.Length - 1] = 0;
                for (int v = DAG.Length - 1; v >= 0; v--)
                {
                    if (DAG.Length - v < i || S[v]) pathsStartingAtNodeVofLengthi[v] = 0;
                    else if (DAG.Length - v == i)
                    {
                        pathsStartingAtNodeVofLengthi[v] = 0;
                        if (DAG[v] < v - 1 && S[DAG[v]] == false) pathsStartingAtNodeVofLengthi[DAG[v]] += pathsStartingAtNodeVofLengthiMinus1[v];
                    }
                    else
                    {
                        pathsStartingAtNodeVofLengthi[v] += S[v + 1] ? 0 : pathsStartingAtNodeVofLengthiMinus1[v + 1];
                        if (DAG[v] < v - 1 && S[DAG[v]] == false) pathsStartingAtNodeVofLengthi[DAG[v]] += pathsStartingAtNodeVofLengthiMinus1[v];
                    }

                }

                while (biggestDistanceInStack > d - i) { if (DistDPathCountEndAtV.Count > 1) { DistDPathCountEndAtV.Pop(); biggestDistanceInStack = DistDPathCountEndAtV.Peek().Item2; } else biggestDistanceInStack = 0; }
                Tuple<double[], int> T = DistDPathCountEndAtV.Peek();
                if (T.Item2 >= d - i) { DistDPathCountEndAtV.Pop(); if (DistDPathCountEndAtV.Count > 0) biggestDistanceInStack = DistDPathCountEndAtV.Peek().Item2; else biggestDistanceInStack = 0; }
                pathsEndingAtNodeVofLengthi = ConsoleApplication1.Helpers.copyArray(T.Item1);
                pathsEndingAtNodeVofLengthiMinus1 = ConsoleApplication1.Helpers.copyArray(T.Item1);
                // @Lu recompute the path number end at v
                for (int y = T.Item2; y < d - i; y++)
                {
                    pathsEndingAtNodeVofLengthi[0] = 0;

                    for (int v = 1; v < DAG.Length; v++)
                    {

                        if (S[v]) pathsEndingAtNodeVofLengthi[v] = 0;// = pathsEndingAtNodeVofLengthiMinus1[v - 1];
                        else pathsEndingAtNodeVofLengthi[v] = (S[v - 1] ? 0 : pathsEndingAtNodeVofLengthiMinus1[v - 1])
                                + ((DAG[v] < v - 1 && S[DAG[v]] == false) ? pathsEndingAtNodeVofLengthiMinus1[DAG[v]] : 0);
                    }
                    if (y + 1 >= (d - i) - ((d - i) - biggestDistanceInStack) / 2 && y < d - i)
                    {
                        DistDPathCountEndAtV.Push(new Tuple<double[], int>(ConsoleApplication1.Helpers.copyArray(pathsEndingAtNodeVofLengthi), y + 1));
                        biggestDistanceInStack = y + 1;
                    }
                    for (int v = 0; v < DAG.Length && y + 1 < d - i; v++)
                    {
                        pathsEndingAtNodeVofLengthiMinus1[v] = pathsEndingAtNodeVofLengthi[v];
                    }
                }
                for (int v = 0; v < DAG.Length; v++)
                {
                    //@Lu miss length d paths that start at v but this means v has no predecessors
                    pathsIncidentToV[v] += pathsEndingAtNodeVofLengthi[v] * pathsStartingAtNodeVofLengthi[v];

                }

                for (int v = 0; v < DAG.Length; v++)
                {
                    pathsStartingAtNodeVofLengthiMinus1[v] = pathsStartingAtNodeVofLengthi[v];
                    pathsStartingAtNodeVofLengthi[v] = 0;
                    pathsEndingAtNodeVofLengthi[v] = pathsEndingAtNodeVofLengthiMinus1[v];
                }
            }
            // int vMax = 0;
            // for (int v = 0; v < DAG.Length; v++)
            // {
            //     if (maxVNumPathsOfLengthdIncidentToV < pathsIncidentToV[v])
            //     {
            //         maxVNumPathsOfLengthdIncidentToV = pathsIncidentToV[v];
            //         vMax = v;
            //     }
            //     maxVNumPathsOfLengthdIncidentToV = Math.Max(maxVNumPathsOfLengthdIncidentToV, pathsIncidentToV[v]);
            // }
            // @Lu
            return new Tuple<double, double[]>(totalPathsOfLengthd, pathsIncidentToV);
            // return new Tuple<double, double, int>(totalPathsOfLengthd, maxVNumPathsOfLengthdIncidentToV, vMax);
        }

        public static List<int>[] findDAGadj(int[] DAG)
        {
            // based the DAG build a reversed graph storing all the successors
            List<int>[] dag2 = new List<int>[DAG.Length-1];
            int i;
            
            for (i = 0; i < DAG.Length - 1; i++)
            {
                dag2[i] = new List<int>();
                dag2[i].Add(i+1);
                if (DAG[i] != i-1)
                {
                    dag2[DAG[i]].Add(i);    
                }
            }   
            dag2[DAG[i]].Add(i);
            return dag2;
        }

        public static int[] BFSfromVDAGadj(int[] DAG, List<int>[] DAGadj, int v, bool[] S, int radius, int[] inRadius)
        {            
            // using G and Grev to do BFS starting from node v
            Queue<int> qHoldVertex = new Queue<int>();
            int currVertex = v, r;
            int[] distFromV = new int[DAG.Length];
            for (int i = 0; i < DAG.Length; i++)
                distFromV[i] = -1;
            
            inRadius[v] = 1;
            distFromV[v] = 0;
            qHoldVertex.Enqueue(v);
            r = 1;

            while (qHoldVertex.Count != 0)
            {
                currVertex = qHoldVertex.Dequeue();
                
                r = distFromV[currVertex];
                if (r >= radius)
                    continue;
                
                // v's predecessors using DAG
                if (currVertex > 0 && !S[currVertex - 1])
                {
                    if (distFromV[currVertex - 1] == -1 || distFromV[currVertex - 1] > r + 1)
                    {
                        distFromV[currVertex - 1] = r + 1;
                        inRadius[currVertex - 1] = 1;
                        qHoldVertex.Enqueue(currVertex - 1);
                    }
                }
                else if (!S[DAG[currVertex]])
                {
                    if (distFromV[DAG[currVertex]] == -1 || distFromV[DAG[currVertex]] > r + 1)
                    {
                        distFromV[DAG[currVertex]] = r + 1;
                        inRadius[DAG[currVertex]] = 1;
                        qHoldVertex.Enqueue(DAG[currVertex]);
                    }
                }
                
                // v's successors using DAGadj
                if (currVertex == DAG.Length - 1)// this node is the last node in DAG then it has no child
                    continue;
                List<int> vAdj = DAGadj[currVertex];
                foreach (int u in vAdj)
                {
                    if (!S[u] && (distFromV[u] == -1 || distFromV[u] > r + 1))
                    {
                        distFromV[u] = r + 1;
                        inRadius[u] = 1;
                        qHoldVertex.Enqueue(u);
                    }
                }
            }
            return inRadius;
        }

	}
}