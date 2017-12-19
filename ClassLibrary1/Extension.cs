using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace Atkin_MoreinPrimeTest
{
    static class Extension
    {
        /// <summary>
        /// Символ Якоби(Для простых m - символ Лежандра)
        /// </summary>
        /// <param name="a">Вычет</param>
        /// <param name="m">Модуль</param>
        /// <returns></returns>
        public static int Lezhandr(BigInteger a, BigInteger m)
        {
            a = BigInteger.Remainder(a, m);
            int t = 1;
            while (a.CompareTo(0)!=0)
            {
                while(BigInteger.Remainder(a,2) == 0)
                {
                    a = BigInteger.Divide(a, 2);
                    if (BigInteger.Remainder(m, 8) == 3 || BigInteger.Remainder(m, 8) == 5)
                        t = -t;
                }
                BigInteger temp = a;
                a = m;
                m = temp;
                if (BigInteger.Remainder(a, 4) == 3 && BigInteger.Remainder(m, 4) == 3)
                    t = -t;
                a = BigInteger.Remainder(a, m);
            }
            if (m.CompareTo(1) == 0)
                return t;

            return 0;
        }
        /// <summary>
        /// Квадратный корень по простому модулю. x^2 mod p = a
        /// </summary>
        /// <param name="a">Целое число</param>
        /// <param name="p">Модуль</param>
        /// <returns>Если -1 то решения нет</returns>
        public static BigInteger SquareRootModPrime(BigInteger a,BigInteger p)
        {
            a = BigInteger.Remainder(a, p);
            while (a < 0) a += p; // не эффективно

            if (Lezhandr(a, p) != 1) return -1;

            
            //----------------Простейшие-случаи p mod 8 = 3,5,7--------------------------------------------
            if (BigInteger.Remainder(p, 8)==3 || BigInteger.Remainder(p, 8)==7)
            {
                BigInteger x = BigInteger.ModPow(a, (p + 1) / 4, p);
                return x;      
            }

            if (BigInteger.Remainder(p, 8) == 5)
            {
                BigInteger x = BigInteger.ModPow(a, (p + 3) / 8, p);
                BigInteger y = BigInteger.ModPow(x, 2, p);
                if (y != a) x = BigInteger.Remainder(x * BigInteger.ModPow(2, (p - 1) / 4, p), p);
                return x;             
            }
            //-------------------------Общий случай----------------------------------------------------------------------
            BigInteger b = 1;
            while (Lezhandr(b, p) != -1) b++;

            int s = 0; BigInteger t = p - 1;
            while (BigInteger.Remainder(t, 2) == 0)
            {
                t = BigInteger.Divide(t, 2);
                s++;
            } // разложение 2^s * t

            BigInteger A = Inverse(a, p);
            BigInteger c = BigInteger.ModPow(b, t, p);
            BigInteger r = BigInteger.ModPow(a,(t+1)/2,p);
            for(int i = 1; i< s;i++)
            {
                BigInteger exp = BigInteger.Pow(2, s - 1 - i);
                BigInteger d = BigInteger.ModPow(r * r * A, exp, p);
                if (d == p - 1) r = BigInteger.Remainder(r * c, p);
                c = BigInteger.ModPow(c, 2, p);
            }
            return r;
            //if (r + r > p) r = BigInteger.Abs(r - p);
            //return r + r > p ? BigInteger.Abs(r - p) : r;
        }
        /// <summary>
        /// Целочисленный квадратный корень Если -1 то ошибка
        /// </summary>
        /// <param name="n"></param>
        /// <param name="index"></param>
        /// <returns></returns>
        public static BigInteger Sqrt(BigInteger n, int index = 2)
        {
            if (n <= 0) return 0;
            BigInteger x = BigInteger.Divide(n, 3) + 1;
            bool flag = false;
            switch (index)
            {
                case 2:
                    while (true)
                    {
                        BigInteger nx = (x + BigInteger.Divide(n, x)) >> 1;
                        if (x == nx || nx > x & flag) break;
                        flag = nx < x;
                        x = nx;
                    }

                    return x;
                case 3:                     
                    while (true)
                    {
                        BigInteger nx = (BigInteger.Divide(2*x,3) + BigInteger.Divide(n, 3*BigInteger.Pow(x,2)));
                        if (x == nx || nx > x & flag) break;
                        flag = nx < x;
                        x = nx;
                    }
                    return x; // по тестам корень на 1 меньше получается 
            }
            return -1;          
        }
        /// <summary>
        ///  представление 4p = x^2 + |D|y^2
        /// </summary>
        /// <param name="p"></param>
        /// <param name="D"></param>
        /// <returns></returns>
        public static List<BigInteger> KornakiSmit(BigInteger p, int D)
        {
            var xy = new List<BigInteger>();
            if (p.CompareTo(2) == 0)
            {
                if(D+8 == BigInteger.Pow(Sqrt(D + 8),2))
                {
                    xy.Add(Sqrt(D + 8));
                    xy.Add(1);
                    return xy;
                }
                return xy;
            }
            BigInteger D1 = p + D;
            if (Lezhandr(D1, p) < 1) return xy; // D отрицательное
            BigInteger x0 = SquareRootModPrime(D1, p);
            if (BigInteger.Remainder(x0,2) != BigInteger.Remainder(Math.Abs(D),2)) x0 = p - x0;

            BigInteger a = 2 * p;
            BigInteger b = x0;
            BigInteger c = 2 * Sqrt(p);

            while(b>c)
            {
                BigInteger temp = a;
                a = b;
                b = BigInteger.Remainder(temp, b);
            }

            BigInteger t = 4 * p - BigInteger.Pow(b, 2);
            if (BigInteger.Remainder(t, Math.Abs(D)) != 0) return xy;
            if (t/ Math.Abs(D) != BigInteger.Pow(Sqrt(t/Math.Abs(D)), 2))
                return xy;

            xy.Add(b);
            xy.Add(Sqrt(BigInteger.Divide(t, Math.Abs(D))));
            return xy;
        }
        /// <summary>
        /// Гильбертов многочлен, число классов и множество привиденных форм (a,b,c) дискриминанта D
        /// </summary>
        /// <param name="D"></param>
        /// <returns></returns>
        public static ComplexPolynom GilbertPolynom(int D)
        {
            //-----Инициализация------------------------------------------

            ComplexPolynom T = new ComplexPolynom(new List<ComplexMonom> {new ComplexMonom(0,1) });
            int b = D & 0x01; // D mod 2
            int k = (int)Math.Sqrt(Math.Abs(D) / 3);
            int h = 0; // счетчик классов
            List<BigInteger[]> red = new List<BigInteger[]>(); // Множество примитивных приведенных форм (a,b,c)

            //----------Основной Цикл по b---------------------------------------

            while(b <= k)
            {
                int t = (b * b - D);
                if ((t & 0x03) != 0)
                {
                    b++;
                    continue;
                }
                int m = (b * b - D) / 4;
                for(int a = 1; a * a <= m;a++)
                {
                    if (m % a != 0) continue;
                    int c = m / a;
                    if (b > a) continue;

                    //---------Установки для многочлена-------------------------      
                    Complex r = new Complex(-b, Math.Sqrt(Math.Abs(D)));
                    r /= 2 * a;
                    Complex delt1 = Delta(Complex.Exp(4 * Math.PI * r * Complex.ImaginaryOne), 5000);
                    Complex delt2 = Delta(Complex.Exp(2 * Math.PI * r * Complex.ImaginaryOne), 5000);
                    var f = delt1 / delt2;
                    var j = Complex.Pow(256 * f + 1, 3);
                    j /= f;

                    //----------------------------------
                    if (b==a || c==a || b==0)
                    {
                        T = ComplexPolynom.MultiplyPol(T, j, "line");
                        h++;
                        red.Add(new BigInteger[] { a, b, c });
                    }
                    else
                    {
                        T = ComplexPolynom.MultiplyPol(T, j, "quadr");
                        h += 2;
                        red.Add(new BigInteger[] { a, b, c });
                        red.Add(new BigInteger[] { a, -b, c });
                    }
                }
                b++;
            }
            return T;

        } 

        private static Complex Delta(Complex q,int epsilon) // delta(q) epsilon - точность
        {
            Complex temp = new Complex(1,0);
            int i = 1;
            while(i<=epsilon)
            {
                temp += Math.Pow(-1, i) * (Complex.Pow(q, (i * (3 * i - 1)) >> 1) + Complex.Pow(q, (i * (3 * i + 1)) >> 1));
                i++;
            }
            temp = Complex.Pow(temp, 24);
            temp *= q;

            return temp;
        }

        public static BigInteger Inverse(BigInteger ch, BigInteger n)
        {
            if (ch == 1) return 1;
            BigInteger a = ch, y = 0;
            BigInteger b = n, x = BigInteger.Zero, d = BigInteger.One;
            while (a.CompareTo(BigInteger.Zero) == 1)//a>0
            {
                BigInteger q = BigInteger.Divide(b, a);
                y = a;
                a = BigInteger.Remainder(b, a);
                b = y;
                y = d;
                d = BigInteger.Subtract(x, BigInteger.Multiply(q, d));
                x = y;
            }
            x = BigInteger.Remainder(x, n);
            if (x.CompareTo(BigInteger.Zero) == -1)//x<0
            {
                x = BigInteger.Remainder(BigInteger.Add(x, n), n);
            }
            return x;
        } // ch^-1 mod n

        public static BigInteger Random (BigInteger A, BigInteger B)
        {
            Random rnd = new Random(DateTime.Now.Millisecond);
            int countByte = B.ToByteArray().Count();
            byte[] r = new byte[countByte + rnd.Next(-1, 2)];
            rnd.NextBytes(r);
            BigInteger randBigint = new BigInteger(r);
            return BigInteger.Abs(BigInteger.Remainder(randBigint, B - A)) + A;
        }

        
    }
}
