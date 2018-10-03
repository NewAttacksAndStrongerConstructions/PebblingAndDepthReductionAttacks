using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ConsoleApplication1
{
    class CustomDouble
    {
        double contents = 0.0;

        CustomDouble(double n)
        {
            contents = n;// Math.Log10(n);
        }

        public static implicit operator CustomDouble(double n)
        {
            return new CustomDouble(n);
        }

        public static double getOri(CustomDouble a)
        {
            return Math.Pow(10, a.contents);
        }

        public static double getContents(CustomDouble a)
        {
            return a.contents;
        }

        public static CustomDouble operator +(CustomDouble a, CustomDouble b)
        {
            
            double temp;
            CustomDouble res;
            if (a.contents == Double.NegativeInfinity && b.contents == Double.NegativeInfinity)
            {
                res = Double.NegativeInfinity;
            }
            else if (a.contents > b.contents)
            {
                temp = (a.contents + Math.Log10(1 + Math.Pow(10, b.contents - a.contents)));
                //if (temp > Double.PositiveInfinity) Console.Write(">");
                res = temp;// Math.Pow(10,temp);
            }
            else
            {
                temp = b.contents + Math.Log10(1 + Math.Pow(10, a.contents - b.contents));
                //if (temp > Double.PositiveInfinity) Console.Write(">");
                res = temp;// Math.Pow(10, temp);
            }
            //Console.WriteLine(CustomDouble.getContents(res));
            return res;
        }

        public static CustomDouble operator *(CustomDouble a, CustomDouble b)
        {
            double temp;
            CustomDouble res;
            temp = a.contents + b.contents;
            // if (temp > Double.PositiveInfinity) Console.Write(">");
            res = temp;// Math.Pow(10, temp);
            return res;
        }

        public static Boolean operator >(CustomDouble a, CustomDouble b)
        {
            if (a.contents > b.contents) return true;
            else return false;
        }

        public static Boolean operator <(CustomDouble a, CustomDouble b)
        {
            if (a.contents < b.contents) return true;
            else return false;
        }
    }
}
