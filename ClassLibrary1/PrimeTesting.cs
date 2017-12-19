using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;

namespace Atkin_MoreinPrimeTest
{
    public class PrimeTesting
    {
        #region Свойства
        /// <summary>
        /// Дискриминанты с числом классов 1 и 2 по возрастанию
        /// </summary>
        private int[] D = new int [] {-3,-4,-7,-8,-11,-19,-43,-67,-163,-15,-20,-24,-35,-40,-51,-52,-88,-91,-115,-123,-148,-187,-232,-235,-267,-403,-427};
    
        /// <summary>
        /// Число для тестрования на простоту
        /// </summary>
        public static BigInteger Number
        {
            get;
            private set;
        }
        /// <summary>
        /// Класс разложения на множители
        /// </summary>
        private FactorOrders fct;
        private StringBuilder Cert;
        private Dictionary<int, ComplexPolynom> GilbertPolynoms;
        /// <summary>
        /// Сертификат числа
        /// </summary>
        public string Certificate
        {
            get
            {
                return Cert.ToString();
            }
            private set
            {

            }
        }
        #endregion

        #region Конструкторы
    
        public PrimeTesting()
        {
            Cert = new StringBuilder();
            fct = new FactorOrders();
            GilbertPolynoms = new Dictionary<int, ComplexPolynom>();
        }
        #endregion

        #region Методы
        public void StartTesting()
        {
            if (!FactorOrders.MRTest(Number, new List<int> { 2, 3, 5, 7, 11}))
            {
                Cert = new StringBuilder("Composite");
                return;
            }
            BigInteger res = Testing();  
            while(res!=-1 && res!=1)
            {
                if (res == -2) break;
                Number = res;
                res = Testing();
            }      
            if(res == -1)
            {
                Cert = new StringBuilder("Composite");
            } else if(res == -2)
            {
                Cert = new StringBuilder("Error");
            }      
            else if(res == 1)
            {
                Cert.Append(string.Format("{0} is Prime", Number));
            }    
        }
        public void SetNewNumber(BigInteger N)
        {
            Cert = new StringBuilder();
            Number = N;
        }

        BigInteger Testing(int _currentD = 0)
        {
            //--наименьшее псевдопростое по 4 базам
            if (Number < 3215031751)
                if (FactorOrders.MRTest(Number, new List<int> { 2, 3, 5, 7, 11 }))
                    return 1;
                else
                    return -1;

            if(_currentD == D.Length) //Кончились дискриминанты
            {
                //D = FundamentalDiscrim(9000,D);
                //return Testing(_currentD);
                return -2;
            }
            int currentD = _currentD;
            //--------------L(D,N) = 1------------------------------------
            while (Extension.Lezhandr(D[currentD], Number) != 1)
            {
                currentD++;
                if (currentD == D.Length)
                {
                    //D = FundamentalDiscrim(9000,D);
                    //return Testing(_currentD);
                    return -2;
                }// Кончились дискриминанты
            }
            //-----------Пытаемся-найти-представление-4p=u^2+|D|v^2-----------------------   
            List<BigInteger> uv = new List<BigInteger>();
            while (currentD < D.Length)
            {
                 uv = Extension.KornakiSmit(Number, D[currentD]);
                if(uv.Count==0)
                    currentD++;
                else       
                    break;              
            }
            if (currentD == D.Length)
            {
                //D = FundamentalDiscrim(9000,D);
                //return Testing(_currentD);
                return -2;
            }// Кончились дискриминанты
            //-----------------------------------------------------------------------
            //-------------------Получаем возможные порядки------------------------------
            List<BigInteger> ordersCurve = new List<BigInteger>();
            if (D[currentD] == -3) // 6 порядков
            {               
                ordersCurve.Add(Number + 1 + ((uv[0] + 3 * uv[1]) >> 1));
                ordersCurve.Add(Number + 1 - ((uv[0] + 3 * uv[1]) >> 1));
                ordersCurve.Add(Number + 1 + ((uv[0] - 3 * uv[1]) >> 1));
                ordersCurve.Add(Number + 1 - ((uv[0] - 3 * uv[1]) >> 1));
                ordersCurve.Add(Number + 1 + uv[0]);
                ordersCurve.Add(Number + 1 - uv[0]);
            }
            else if(D[currentD] == -4 )// 4 порядка
            {
                ordersCurve.Add(Number + 1 + 2 * uv[1]);
                ordersCurve.Add(Number + 1 - 2 * uv[1]);
                ordersCurve.Add(Number + 1 + uv[0]);
                ordersCurve.Add(Number + 1 - uv[0]);
            }
            else
            {
                ordersCurve.Add(Number + 1 + uv[0]);
                ordersCurve.Add(Number + 1 - uv[0]);
            }
            //-----------------------------------------------------------------------
            //-----------------Раскладываем порядки на множители---------------------

            fct.SetNewNumbers(ordersCurve);
            while (true)
            {
                fct.StartFact();
                if (fct.result != null)
                {
                    //---------------Если попал простой порядок-----------------------
                    if(fct.result.Count == 1)
                    {
                        ordersCurve.RemoveRange(0, fct.NumberResult + 1);
                        fct.SetNewNumbers(ordersCurve);
                        continue;
                    }
                    BigInteger orderCurve = ordersCurve[fct.NumberResult];
                    //-----------------Получаем параметры кривой---------------------
                    var curveParam = GetCurveComplexMultiply(D[currentD]);
                  //  var curveParam = GetClearCurveComplexMultiply(D[currentD]);

                    if (curveParam == null) return -2;// Проблемы с нахождением корней многочлена
                    var paramCurve = GetCurveParamForOrder(orderCurve, curveParam);
                    if (paramCurve == null)
                        return -1; // Составное
                    if (paramCurve.Length == 1)
                    {
                        ordersCurve.RemoveRange(0, fct.NumberResult + 1);
                        fct.SetNewNumbers(ordersCurve);
                        continue;
                    }

                    //---------------------------------------------------------------
                    //-------------Операции с точкой---------------------------------
                    EllipticCurvePoint P = new EllipticCurvePoint(Number, paramCurve[0], paramCurve[1]);
                    BigInteger k = orderCurve / fct.result[fct.result.Count - 1];
                    var U = EllipticCurvePoint.Multiply(k, P);
                    if (U.CheckPoint() == false)
                        return -1; // Составное
                    while (U.X == 0 && U.Z == 0)
                    {
                        P.NextPoint();
                        U = EllipticCurvePoint.Multiply(k, P);
                        if (U.CheckPoint() == false)
                            return -1; // Составное
                    }
                    var V = EllipticCurvePoint.Multiply(fct.result[fct.result.Count - 1], U);
                    if (V.X != 0 && V.Z != 0)
                        return -1; // Составное
                                   //--------------Формирование Сертификата---------------------

                    Cert.Append(string.Format("N = {0}",Number)+ "\r\n");
                    Cert.Append(string.Format("D = {0}  U = {1}  V = {2}",D[currentD],uv[0],uv[1])+ "\r\n");
                    Cert.Append("#E = ");
                    foreach(var del in fct.result)
                    {
                        Cert.Append(del + " * ");
                    }
                    Cert.Remove(Cert.Length - 3, 3);
                    Cert.Append("\r\n");
                    Cert.Append(string.Format("E({0}, {1}, {2}, {3}, {4}, {5} )", P.A, P.B, P.P, P.X, P.Y, P.Z)+ "\r\n");
                    Cert.Append("\r\n");
                    //---------------------------------------------------------------
                    //-----------Возврат числа q-------------------------------------

                    return fct.result[fct.result.Count - 1];

                    //---------------------------------------------------------------
                }
                else // если порядки разложить не удалось, то начинаем с начала со следующего дискриминанта
                {
                    return Testing(currentD + 1);
                }
            }   
        }
        #endregion

        #region privateMethods
        /// <summary>
        /// Параметры кривых полученных комплексным умножением
        /// </summary>
        /// <param name="D">Фундаментальный Дискриминант</param>
        /// <returns></returns>
        private  List<BigInteger[]> GetCurveComplexMultiply(int D)
        {
            List<BigInteger[]> paramCurve = new List<BigInteger[]>();
            BigInteger g = GetNonQuadr();
            if(D == -3)
            {
                for(int i = 0; i < 6; i++)
                {
                    paramCurve.Add(new BigInteger[] {0,-BigInteger.ModPow(g,i,Number)});
                }
                return paramCurve;            
            }
            else if (D == -4)
            {
                for (int i = 0; i < 4; i++)
                {
                    paramCurve.Add(new BigInteger[] { -BigInteger.ModPow(g, i, Number), 0 });
                }
                return paramCurve;
            }
            else
            {
                Polynom GilbPolFp;
                if (GilbertPolynoms.ContainsKey(D))
                {
                    ComplexPolynom GilbPol;
                    GilbertPolynoms.TryGetValue(D, out GilbPol);
                    GilbPolFp = GilbPol.GetPolynomReal(Number);
                }
                else
                {
                    ComplexPolynom GilbPol = Extension.GilbertPolynom(D);
                    GilbertPolynoms.Add(D, GilbPol);
                    GilbPolFp = GilbPol.GetPolynomReal(Number);
                }
                //---------------------------------------------------------------
                //-----------------Получение корня полинома----------------------   

                BigInteger root;
                if (GilbPolFp.Degree <= 2)
                {
                    var roots = GilbPolFp.GetRoots();
                    if (roots.Count != 0)
                        root = roots[0];
                    else
                    {
                        return null;
                    }
                }
                else
                {
                    #region Поиск нод 
                    //var nod = Polynom.GreatCommonDivisor(GilbPolFp, Polynom.Derivative(GilbPolFp));
                    //if(nod.Degree!=0 && nod.Degree<3)
                    //{
                    //    var roots = nod.GetRoots();
                    //    if (roots.Count != 0)
                    //        root = roots[0];
                    //    else
                    //    {
                    //        return null;
                    //    }
                    //}
                    //else
                    //{
                    //    using (StreamWriter sw = new StreamWriter(new FileStream("polynom.txt", FileMode.Append)))
                    //    {
                    //        sw.WriteLine("//------------------------------------------------------------------------");
                    //        sw.WriteLine(Number);
                    //        sw.WriteLine("D = " + D);
                    //        sw.WriteLine(GilbPolFp + " mod " + Number);
                    //        sw.WriteLine("//------------------------------------------------------------------------");
                    //        sw.Close();
                    //    }
                    //    return null;
                    //}   
                    #endregion
                    return null;
                }

                //---------------------------------------------------------------
                //------------------Получение параметров-------------------------
                BigInteger c = root * Extension.Inverse(root - 1728, Number);
                c = BigInteger.Remainder(c, Number);
                BigInteger r = BigInteger.Remainder(3 * BigInteger.Negate(c), Number);
                BigInteger s = BigInteger.Remainder(2 * c, Number);
                //---------------------------------------------------------------
                paramCurve.Add(new BigInteger[] { r, s });
                paramCurve.Add(new BigInteger[] { BigInteger.Remainder(r * BigInteger.ModPow(g, 2, Number), Number), BigInteger.Remainder(s * BigInteger.ModPow(g, 3, Number), Number) });
                return paramCurve;          
                //---------------------------------------------------------------
            }
        }
        /// <summary>
        /// Явные Параметры кривых полученных комплексным умножением
        /// </summary>
        /// <param name="D">Фундаментальный дискриминант</param>
        /// <returns></returns>
        private static List<BigInteger[]> GetClearCurveComplexMultiply(int D)
        {
            Task<BigInteger> t1 = Task<BigInteger>.Factory.StartNew(delegate() { return GetNonQuadr(); });
            List<BigInteger[]> paramCurve = new List<BigInteger[]>();
            switch (D)
            {
                case (-7):
                    BigInteger r = 125;
                    BigInteger s = 189;
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * r * BigInteger.ModPow(s, 3, Number), Number), BigInteger.Remainder(2 * r * BigInteger.ModPow(s, 5, Number), Number) });
                    t1.Wait();
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * r * BigInteger.ModPow(s, 3, Number)* BigInteger.ModPow(t1.Result, 2, Number), Number), BigInteger.Remainder(2 * r * BigInteger.ModPow(s, 5, Number) * BigInteger.ModPow(t1.Result, 3, Number), Number) });
                    break;
                case (-8):
                    BigInteger r0 = 125;
                    BigInteger s0 = 98;
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * r0 * BigInteger.ModPow(s0, 3, Number), Number), BigInteger.Remainder(2 * r0 * BigInteger.ModPow(s0, 5, Number), Number) });
                    t1.Wait();
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * r0 * BigInteger.ModPow(s0, 3, Number) * BigInteger.ModPow(t1.Result, 2, Number), Number), BigInteger.Remainder(2 * r0 * BigInteger.ModPow(s0, 5, Number) * BigInteger.ModPow(t1.Result, 3, Number), Number) });
                    break;
                case (-11):
                    BigInteger r1 = 512;
                    BigInteger s1 = 539;
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * r1 * BigInteger.ModPow(s1, 3, Number), Number), BigInteger.Remainder(2 * r1 * BigInteger.ModPow(s1, 5, Number), Number) });
                    t1.Wait();
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * r1 * BigInteger.ModPow(s1, 3, Number) * BigInteger.ModPow(t1.Result, 2, Number), Number), BigInteger.Remainder(2 * r1 * BigInteger.ModPow(s1, 5, Number) * BigInteger.ModPow(t1.Result, 3, Number), Number) });
                    break;
                case (-19):
                    BigInteger r2 = 512;
                    BigInteger s2 = 512;
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * r2 * BigInteger.ModPow(s2, 3, Number), Number), BigInteger.Remainder(2 * r2 * BigInteger.ModPow(s2, 5, Number), Number) });
                    t1.Wait();
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * r2 * BigInteger.ModPow(s2, 3, Number) * BigInteger.ModPow(t1.Result, 2, Number), Number), BigInteger.Remainder(2 * r2 * BigInteger.ModPow(s2, 5, Number) * BigInteger.ModPow(t1.Result, 3, Number), Number) });
                    break;
                case (-43):
                    BigInteger r3 = 512000;
                    BigInteger s3 = 512001;
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * r3 * BigInteger.ModPow(s3, 3, Number), Number), BigInteger.Remainder(2 * r3 * BigInteger.ModPow(s3, 5, Number), Number) });
                    t1.Wait();
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * r3 * BigInteger.ModPow(s3, 3, Number) * BigInteger.ModPow(t1.Result, 2, Number), Number), BigInteger.Remainder(2 * r3 * BigInteger.ModPow(s3, 5, Number) * BigInteger.ModPow(t1.Result, 3, Number), Number) });
                    break;
                case (-67):
                    BigInteger r4 = 85184000;
                    BigInteger s4 = 85184001;
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * r4 * BigInteger.ModPow(s4, 3, Number), Number), BigInteger.Remainder(2 * r4 * BigInteger.ModPow(s4, 5, Number), Number) });
                    t1.Wait();
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * r4 * BigInteger.ModPow(s4, 3, Number) * BigInteger.ModPow(t1.Result, 2, Number), Number), BigInteger.Remainder(2 * r4 * BigInteger.ModPow(s4, 5, Number) * BigInteger.ModPow(t1.Result, 3, Number), Number) });
                    break;
                case (-163):
                    BigInteger r5 = 151931373056000;
                    BigInteger s5 = 151931373056001;
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * r5 * BigInteger.ModPow(s5, 3, Number), Number), BigInteger.Remainder(2 * r5 * BigInteger.ModPow(s5, 5, Number), Number) });
                    t1.Wait();
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * r5 * BigInteger.ModPow(s5, 3, Number) * BigInteger.ModPow(t1.Result, 2, Number), Number), BigInteger.Remainder(2 * r5 * BigInteger.ModPow(s5, 5, Number) * BigInteger.ModPow(t1.Result, 3, Number), Number) });
                    break;
                case (-15):
                    BigInteger r6 = 1225 - 2080 * Extension.SquareRootModPrime(5,Number);
                    BigInteger s6 = 5929;
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * r6 * BigInteger.ModPow(s6, 3, Number), Number), BigInteger.Remainder(2 * r6 * BigInteger.ModPow(s6, 5, Number), Number) });
                    t1.Wait();
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * r6 * BigInteger.ModPow(s6, 3, Number) * BigInteger.ModPow(t1.Result, 2, Number), Number), BigInteger.Remainder(2 * r6 * BigInteger.ModPow(s6, 5, Number) * BigInteger.ModPow(t1.Result, 3, Number), Number) });
                    break;
                case (-20):
                    BigInteger r7 = 108250 + 29835 * Extension.SquareRootModPrime(5, Number);
                    BigInteger s7 = 174724;
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * r7 * BigInteger.ModPow(s7, 3, Number), Number), BigInteger.Remainder(2 * r7 * BigInteger.ModPow(s7, 5, Number), Number) });
                    t1.Wait();
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * r7 * BigInteger.ModPow(s7, 3, Number) * BigInteger.ModPow(t1.Result, 2, Number), Number), BigInteger.Remainder(2 * r7 * BigInteger.ModPow(s7, 5, Number) * BigInteger.ModPow(t1.Result, 3, Number), Number) });
                    break;
                case (-24):
                    BigInteger r8 = 1757 - 494 * Extension.SquareRootModPrime(2, Number);
                    BigInteger s8 = 1058;
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * r8 * BigInteger.ModPow(s8, 3, Number), Number), BigInteger.Remainder(2 * r8 * BigInteger.ModPow(s8, 5, Number), Number) });
                    t1.Wait();
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * r8 * BigInteger.ModPow(s8, 3, Number) * BigInteger.ModPow(t1.Result, 2, Number), Number), BigInteger.Remainder(2 * r8 * BigInteger.ModPow(s8, 5, Number) * BigInteger.ModPow(t1.Result, 3, Number), Number) });
                    break;
                case (-35):
                    BigInteger r9 = -1126400 - 1589760 * Extension.SquareRootModPrime(5, Number);
                    BigInteger s9 = 2428447;
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * r9 * BigInteger.ModPow(s9, 3, Number), Number), BigInteger.Remainder(2 * r9 * BigInteger.ModPow(s9, 5, Number), Number) });
                    t1.Wait();
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * r9 * BigInteger.ModPow(s9, 3, Number) * BigInteger.ModPow(t1.Result, 2, Number), Number), BigInteger.Remainder(2 * r9 * BigInteger.ModPow(s9, 5, Number) * BigInteger.ModPow(t1.Result, 3, Number), Number) });
                    break;
                case (-40):
                    BigInteger r11 = 54175 - 1020 * Extension.SquareRootModPrime(5, Number);
                    BigInteger s11 = 51894;
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * r11 * BigInteger.ModPow(s11, 3, Number), Number), BigInteger.Remainder(2 * r11 * BigInteger.ModPow(s11, 5, Number), Number) });
                    t1.Wait();
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * r11 * BigInteger.ModPow(s11, 3, Number) * BigInteger.ModPow(t1.Result, 2, Number), Number), BigInteger.Remainder(2 * r11 * BigInteger.ModPow(s11, 5, Number) * BigInteger.ModPow(t1.Result, 3, Number), Number) });
                    break;
                case (-51):
                    BigInteger r12 = 75520 - 7936 * Extension.SquareRootModPrime(17, Number);
                    BigInteger s12 = 108241;
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * r12 * BigInteger.ModPow(s12, 3, Number), Number), BigInteger.Remainder(2 * r12 * BigInteger.ModPow(s12, 5, Number), Number) });
                    t1.Wait();
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * r12 * BigInteger.ModPow(s12, 3, Number) * BigInteger.ModPow(t1.Result, 2, Number), Number), BigInteger.Remainder(2 * r12 * BigInteger.ModPow(s12, 5, Number) * BigInteger.ModPow(t1.Result, 3, Number), Number) });
                    break;
                case (-52):
                    BigInteger rq = 1778750 + 5125 * Extension.SquareRootModPrime(13, Number);
                    BigInteger sq = 1797228;
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * rq * BigInteger.ModPow(sq, 3, Number), Number), BigInteger.Remainder(2 * rq * BigInteger.ModPow(sq, 5, Number), Number) });
                    t1.Wait();
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * rq * BigInteger.ModPow(sq, 3, Number) * BigInteger.ModPow(t1.Result, 2, Number), Number), BigInteger.Remainder(2 * rq * BigInteger.ModPow(sq, 5, Number) * BigInteger.ModPow(t1.Result, 3, Number), Number) });
                    break;
                case (-88):
                    BigInteger rw = 181713125 - 44250 * Extension.SquareRootModPrime(2, Number);
                    BigInteger sw = 181650546;
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * rw * BigInteger.ModPow(sw, 3, Number), Number), BigInteger.Remainder(2 * rw * BigInteger.ModPow(sw, 5, Number), Number) });
                    t1.Wait();
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * rw * BigInteger.ModPow(sw, 3, Number) * BigInteger.ModPow(t1.Result, 2, Number), Number), BigInteger.Remainder(2 * rw * BigInteger.ModPow(sw, 5, Number) * BigInteger.ModPow(t1.Result, 3, Number), Number) });
                    break;
                case (-91):
                    BigInteger re = 74752 - 36352 * Extension.SquareRootModPrime(13, Number);
                    BigInteger se = 205821;
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * re * BigInteger.ModPow(se, 3, Number), Number), BigInteger.Remainder(2 * re * BigInteger.ModPow(se, 5, Number), Number) });
                    t1.Wait();
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * re * BigInteger.ModPow(se, 3, Number) * BigInteger.ModPow(t1.Result, 2, Number), Number), BigInteger.Remainder(2 * re * BigInteger.ModPow(se, 5, Number) * BigInteger.ModPow(t1.Result, 3, Number), Number) });
                    break;
                case (-115):
                    BigInteger ra = 269593600 - 89157120 * Extension.SquareRootModPrime(5, Number);
                    BigInteger sa = 468954981;
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * ra * BigInteger.ModPow(sa, 3, Number), Number), BigInteger.Remainder(2 * ra * BigInteger.ModPow(sa, 5, Number), Number) });
                    t1.Wait();
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * ra * BigInteger.ModPow(sa, 3, Number) * BigInteger.ModPow(t1.Result, 2, Number), Number), BigInteger.Remainder(2 * ra * BigInteger.ModPow(sa, 5, Number) * BigInteger.ModPow(t1.Result, 3, Number), Number) });
                    break;
                case (-123):
                    BigInteger rs = 1025058304000 - 1248832000 * Extension.SquareRootModPrime(41, Number);
                    BigInteger ss = 1033054730449;
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * rs * BigInteger.ModPow(ss, 3, Number), Number), BigInteger.Remainder(2 * rs * BigInteger.ModPow(ss, 5, Number), Number) });
                    t1.Wait();
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * rs * BigInteger.ModPow(ss, 3, Number) * BigInteger.ModPow(t1.Result, 2, Number), Number), BigInteger.Remainder(2 * rs * BigInteger.ModPow(ss, 5, Number) * BigInteger.ModPow(t1.Result, 3, Number), Number) });
                    break;
                case (-148):
                    BigInteger rd = 499833128054750 + 356500625 * Extension.SquareRootModPrime(37, Number);
                    BigInteger sd = 499835296563372;
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * rd * BigInteger.ModPow(sd, 3, Number), Number), BigInteger.Remainder(2 * rd * BigInteger.ModPow(sd, 5, Number), Number) });
                    t1.Wait();
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * rd * BigInteger.ModPow(sd, 3, Number) * BigInteger.ModPow(t1.Result, 2, Number), Number), BigInteger.Remainder(2 * rd * BigInteger.ModPow(sd, 5, Number) * BigInteger.ModPow(t1.Result, 3, Number), Number) });
                    break;
                case (-187):
                    BigInteger rz = 91878880000 - 1074017568000 * Extension.SquareRootModPrime(17, Number);
                    BigInteger sz = 4520166756633;
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * rz * BigInteger.ModPow(sz, 3, Number), Number), BigInteger.Remainder(2 * rz * BigInteger.ModPow(sz, 5, Number), Number) });
                    t1.Wait();
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * rz * BigInteger.ModPow(sz, 3, Number) * BigInteger.ModPow(t1.Result, 2, Number), Number), BigInteger.Remainder(2 * rz * BigInteger.ModPow(sz, 5, Number) * BigInteger.ModPow(t1.Result, 3, Number), Number) });
                    break;
                case (-232):
                    BigInteger rx = 1728371226151263375 - 11276414500 * Extension.SquareRootModPrime(29, Number);
                    BigInteger sx = 1728371165425912854;
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * rx * BigInteger.ModPow(sx, 3, Number), Number), BigInteger.Remainder(2 * rx * BigInteger.ModPow(sx, 5, Number), Number) });
                    t1.Wait();
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * rx * BigInteger.ModPow(sx, 3, Number) * BigInteger.ModPow(t1.Result, 2, Number), Number), BigInteger.Remainder(2 * rx * BigInteger.ModPow(sx, 5, Number) * BigInteger.ModPow(t1.Result, 3, Number), Number) });
                    break;
                case (-267):
                    BigInteger rc = BigInteger.Parse("3632253349307716000000") - 12320504793376000 * Extension.SquareRootModPrime(89, Number);
                    BigInteger sc = BigInteger.Parse("3632369580717474122449");
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * rc * BigInteger.ModPow(sc, 3, Number), Number), BigInteger.Remainder(2 * rc * BigInteger.ModPow(sc, 5, Number), Number) });
                    t1.Wait();
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * rc * BigInteger.ModPow(sc, 3, Number) * BigInteger.ModPow(t1.Result, 2, Number), Number), BigInteger.Remainder(2 * rc * BigInteger.ModPow(sc, 5, Number) * BigInteger.ModPow(t1.Result, 3, Number), Number) });
                    break;
                case (-403):
                    BigInteger rcq = BigInteger.Parse("16416107434811840000") - 4799513373120384000 * Extension.SquareRootModPrime(13, Number);
                    BigInteger scq = BigInteger.Parse("33720998998872514077");
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * rcq * BigInteger.ModPow(scq, 3, Number), Number), BigInteger.Remainder(2 * rcq * BigInteger.ModPow(scq, 5, Number), Number) });
                    t1.Wait();
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * rcq * BigInteger.ModPow(scq, 3, Number) * BigInteger.ModPow(t1.Result, 2, Number), Number), BigInteger.Remainder(2 * rcq * BigInteger.ModPow(scq, 5, Number) * BigInteger.ModPow(t1.Result, 3, Number), Number) });
                    break;
                case (-427):
                    BigInteger rcqa = BigInteger.Parse("564510997315289728000") - 5784785611102784000 * Extension.SquareRootModPrime(61, Number);
                    BigInteger scqa = BigInteger.Parse("609691617259594724421");
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * rcqa * BigInteger.ModPow(scqa, 3, Number), Number), BigInteger.Remainder(2 * rcqa * BigInteger.ModPow(scqa, 5, Number), Number) });
                    t1.Wait();
                    paramCurve.Add(new BigInteger[] { BigInteger.Remainder(-3 * rcqa * BigInteger.ModPow(scqa, 3, Number) * BigInteger.ModPow(t1.Result, 2, Number), Number), BigInteger.Remainder(2 * rcqa * BigInteger.ModPow(scqa, 5, Number) * BigInteger.ModPow(t1.Result, 3, Number), Number) });
                    break;
            }
            return paramCurve;
        }
        /// <summary>
        /// Случайный квадратичный невычет по модулю тестируемого числа
        /// </summary>
        /// <returns>Случайный квадратичный невычет</returns>
        private static BigInteger GetNonQuadr()
        {
            if (Number % 3 == 1)
            {
                BigInteger g = Extension.Random(1, Number);
                while ((Extension.Lezhandr(g, Number) != -1) || (BigInteger.ModPow(g,(Number - 1)/3, Number) == 1))
                {
                    g = Extension.Random(1, Number);
                }
                return g;
            }
            else
            {
                BigInteger g = Extension.Random(1, Number);
                while (Extension.Lezhandr(g, Number) != -1)
                {
                    g = Extension.Random(1, Number);
                }
                return g;
            }            
        }
        /// <summary>
        /// Поиск кривой с заданным порядком
        /// </summary>
        /// <param name="orderCurve">Необходимый порядок</param>
        /// <param name="paramCurve">Параметры кривых</param>
        /// <returns></returns>
        private static BigInteger[] GetCurveParamForOrder(BigInteger orderCurve, List<BigInteger[]> paramCurve )
        {
            int count = 0;
            List<int> curveNum = new List<int>(); ;
            int iteration = 0;
            while (iteration < 5)
            {
                curveNum = new List<int>();
                for (int i = 0; i < paramCurve.Count; i++)
                {
                    EllipticCurvePoint point = new EllipticCurvePoint(Number, paramCurve[i][0], paramCurve[i][1]);
                    if (!point.CheckPoint())
                        return null;
                    EllipticCurvePoint orderPoint = EllipticCurvePoint.Multiply(orderCurve, point);
                    if (!orderPoint.CheckPoint())
                        return null;
                    if (orderPoint.X == 0 && orderPoint.Z == 0)
                    {
                        count++;
                        curveNum.Add(i);
                    }
                }
                if (count == 1)
                    return paramCurve[curveNum[0]];
                count = 0;
                iteration++;
            }
            if (curveNum.Count == 0)
                return new BigInteger[] { 0 }; // Порядки не подошли
            return paramCurve[curveNum[1]];
        }
        /// <summary>
        /// Фундаментальные Дискриминанты
        /// </summary>
        /// <param name="border">Граница</param>
        /// <returns>Массив дискриминантов</returns>
        private static int[] FundamentalDiscrim(int border, int[] _oldD)
        {
            List<int> d = _oldD.ToList<int>();
            int[] prime = new int[] { 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271 };
            for (int i = 1; i < border; i++)
            {
                if (i % 16 == 3 || i % 16 == 4 || i % 16 == 7 || i % 16 == 8 || i % 16 == 11 || i % 16 == 15)
                {
                    int t = i;
                    while ((t & 0x1) == 0) t = t >> 1;
                    for (int j = 0; j < prime.Count() && prime[j] < (int)(Math.Sqrt(i)) + 1; j++)
                    {
                        if (t % (prime[j] * prime[j]) == 0)
                        {
                            t = 0;
                            break;
                        }
                    }
                    if (t != 0)
                    {
                        if(!d.Contains(-i))
                        {
                            d.Add(-i);
                        }
                    }
                }
            }
            return d.ToArray();
        }
        #endregion

    }
}
