using System.Diagnostics;
using System.Collections.Generic;
using System.Threading.Tasks;
using System.Threading;
using System.Numerics;

namespace Atkin_MoreinPrimeTest
{
    class FactorOrders
    {
        Task sieveTask =  new Task(() => { Sieve(borderPrimeSieve); });
        private const int borderPrimeSieve = 50000000; // Граница для списка простых чисел
        private static List<ulong> primeNumbers;

        /// <summary>
        /// Список порядков к разложению
        /// </summary>
        List<BigInteger> N;
        /// <summary>
        /// Список границ разложения для порядков
        /// </summary>
        List<BigInteger> borderq;
        /// <summary>
        /// Разложение числа N
        /// </summary>
        public List<BigInteger> result
        {
            get;
            private set;
        }
        /// <summary>
        /// Порядковый номер числа, которое успешно разложилось
        /// </summary>
        public int NumberResult
        {
            get;
            private set;
        }

        public FactorOrders()
        {    
            primeNumbers = new List<ulong>();            
            sieveTask.Start();        
        }

        public void SetNewNumbers(List<BigInteger> Numbers)
        {
            result = null;
            NumberResult = new int();           
            N = new List<BigInteger>();
            borderq = new List<BigInteger>();
            foreach (var _n in Numbers)
            {
                N.Add(_n);
                borderq.Add(BigInteger.Pow(Extension.Sqrt(Extension.Sqrt(_n)) + 1, 2));
            }                  
        }

        public void StartFact() // Разложение по методу кто первый
        {
            if (N.Count == 0)
            {
                result = null;
                return;
            }
            Stopwatch st = new Stopwatch();
            Thread[] factTh = new Thread[N.Count];
            List<BigInteger>[] preresult = new List<BigInteger>[N.Count];
            sieveTask.Wait();
            for (int i = 0; i < factTh.Length; i++)
            {

                factTh[i] = new Thread(new ThreadStart(delegate () { preresult[i] = GetMnozh(N[i], borderq[i]); }));
                factTh[i].Start();
                Thread.Sleep(50); // Если убрать выпадывает исключение из getmnozh. Устал искать
            }
            st.Start();
            bool flagthreadrun = false;
            while (st.ElapsedMilliseconds < 30000) // 30 сек
            {
                flagthreadrun = false;
                for (int i = 0; i < factTh.Length; i++)
                {
                    if (factTh[i].ThreadState == System.Threading.ThreadState.Stopped)
                    {
                        if (preresult[i] != null)
                        {
                            if (preresult[i].Count != 1)
                            {
                                for (int j = 0; j < N.Count; j++)
                                {
                                    if (j != i)
                                        factTh[j].Abort();
                                }
                                result = preresult[i];
                                NumberResult = i;
                                st.Stop();
                                return;
                            }
                        }
                    }
                    else
                        flagthreadrun = true;
                }
                if (!flagthreadrun)
                {
                    st.Stop();
                    return;
                }
                Thread.Sleep(250);
            }
            // Кончилось время но ни один поток не решил задачу и есть работающие последняя проверка
            for (int i = 0; i < N.Count; i++)
            {
                if (factTh[i].ThreadState == System.Threading.ThreadState.Stopped)
                {
                    if (preresult[i] != null)
                        if (preresult[i].Count != 1)
                        {
                            for (int j = 0; j < N.Count; j++)
                            {
                                if (j != i)
                                    factTh[j].Abort();
                            }
                            result = preresult[i];
                            NumberResult = i;
                            st.Stop();
                            return;
                        }
                }
            }
            // Разложить не удалось, отменяем
            for (int j = 0; j < N.Count; j++)
            {
                factTh[j].Abort();
            }
        }



        /// <summary>
        /// разложение числа на множители
        /// </summary>
        /// <returns> </returns>
        static List<BigInteger> GetMnozh(BigInteger N, BigInteger BorderQ)
        {
            List<BigInteger> mnozh = new List<BigInteger>();
            Queue<BigInteger> numbers = new Queue<BigInteger>();
            numbers.Enqueue(N);
            while (numbers.Count != 0)
            {
                BigInteger n = numbers.Dequeue();
                if (MRTest(n, new List<int> { 2, 3, 5, 7, 11 }))
                {
                    mnozh.Add(n);
                }
                else
                {
                    BigInteger p = GetDivisor(n);
                    if (p == n || p == 1)
                        return null;
                    numbers.Enqueue(p);
                    numbers.Enqueue(n / p);
                }
            }       
            if (mnozh[mnozh.Count - 1] < BorderQ) return null;
            return mnozh;
        }


        static BigInteger GetDivisor(BigInteger n)
        {
            for (int i = 0; i < primeNumbers.Count; i++)
            {
                if (BigInteger.Remainder(n, primeNumbers[i]) == 0)
                {
                    return primeNumbers[i];
                }
            }        
            BigInteger p = PollardFact(n);
            return p;

        }
        static BigInteger PollardFact(BigInteger N, int _x0 = 1, int _y0 = 1)
        {
            BigInteger del = 1;
            BigInteger x0 = _x0;
            BigInteger y0 = _y0;
            BigInteger x1, y1;

            while (del.CompareTo(1) == 0 || del.CompareTo(N) == 0)
            {
                x1 = Func(x0, N);
                y1 = Func(Func(y0, N), N);
                del = BigInteger.GreatestCommonDivisor(BigInteger.Abs(BigInteger.Subtract(x1, y1)), N);
                x0 = x1;
                y0 = y1;
            }
            return del;
        }
        static BigInteger Func(BigInteger x, BigInteger N)
        {
            BigInteger q = BigInteger.Pow(x, 2) + x + 1;
            q = BigInteger.Remainder(q, N);
            while (q.CompareTo(0) == -1) q += N;
            return q;
        }


        /// <summary>
        /// Решето Эратосфена
        /// </summary>
        /// <param name="N">Верхняя граница</param>
        private static void Sieve(ulong N)
        {
            bool[] A = new bool[N];
            for (ulong i = 0; i < N; i++)
            {
                A[i] = true;
            }
            for (ulong i = 2; i * i < N; i++)
            {
                if (A[i] == true)
                {
                    for (ulong j = i * i; j < N; j = j + i)
                    {
                        A[j] = false;
                    }
                }
            }
            for (ulong i = 2; i < N; i++)
            {
                if (A[i] == true)
                {
                    primeNumbers.Add(i);
                }
            }
        }
        /// <summary>
        /// Тест Миллера-Рабина :: MRTest(1373653, new List<int> { 2,3 })
        /// </summary>
        /// <param name="n">число</param>
        /// <param name="a">список баз для проверки</param>
        /// <returns></returns>
        public static bool MRTest(BigInteger n, List<int> a)
        {
            if (n == 1 || n == 2 || n == 3 || n == 5 || n == 7 || n == 11) return true; // для корректной проверки чисел 2 3 5 7 по всем базам 2 3 5 7
            if ((n & 0x1) == 0) return false;
            BigInteger n1 = BigInteger.Subtract(n, BigInteger.One); // n-1
            BigInteger s = 0, t = n1;
            while ((t & 0x01) == 0)
            {
                t = t >> 1;
                s = BigInteger.Add(s, BigInteger.One);
            } // разложение 2^s * t
            foreach (int basemr in a)
            {

                BigInteger b = BigInteger.ModPow(basemr, t, n);


                if (b.CompareTo(BigInteger.One) == 0 || b.CompareTo(n1) == 0)
                {
                    continue;
                }
                for (int i = 0; i < s; i++)
                {
                    b = BigInteger.ModPow(b, 2, n);
                    if (b.CompareTo(BigInteger.One) == 0)
                    {
                        return false;
                    }
                    if (b.CompareTo(n1) == 0)
                    {
                        break;
                    }
                }
                if (b.CompareTo(n1) != 0)
                {
                    return false;
                }
            }
            return true;
        }
    }
}
