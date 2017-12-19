using System;
using System.Numerics;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Atkin_MoreinPrimeTest
{
    class Polynom
    {
        #region свойства
        public BigInteger Degree
        {
            get;
            private set;
        }
        public BigInteger Fp // Поле Галуа ()
        {
            get;
            private set;
        }       
        List<Monom> coefficients;
        List<BigInteger> roots; // корни многочлена
        #endregion

        #region Конструкторы
        Polynom(Polynom a)
        {
            Degree = a.Degree;
            Fp = a.Fp;
            coefficients = a.coefficients;
        }
        public Polynom(List<Monom> coeff, BigInteger _Fp)
        {
            coefficients = coeff;
            Fp = _Fp;
            CheckPoly();
        }

        #endregion

        #region Математические Операции с Полиномами

        public static Polynom Derivative(Polynom x)
        {
            List<Monom> t = new List<Monom>();
            foreach(var m in x.coefficients)
            {
                t.Add(new Monom(m.Degree - 1, BigInteger.Remainder(m.Coefficient * m.Degree, x.Fp)));
            }
            return new Polynom(t, x.Fp);
        }
        /// <summary>
        /// Разность полиномов
        /// </summary>
        /// <param name="a">Уменьшаемое</param>
        /// <param name="b">Вычитаемое</param>
        /// <returns></returns>
        public static Polynom operator -(Polynom a, Polynom b)
        {
            a.NormPolynom();
            b.NormPolynom();

            List<Monom> c = new List<Monom>();
            int curnta = 0; // Текущий элемент в полиноме а
            int curntb = 0; // Текущий элемент в полиноме b

            while (curnta < a.coefficients.Count && curntb < b.coefficients.Count)
            {
                if (a.coefficients[curnta].Degree < b.coefficients[curntb].Degree)
                {
                    c.Add(Monom.Negate(b.coefficients[curntb]));
                    curntb++;
                } else if (a.coefficients[curnta].Degree > b.coefficients[curntb].Degree)
                {
                    c.Add(a.coefficients[curnta]);
                    curnta++;
                } else
                {
                    c.Add(a.coefficients[curnta].SubNumber(b.coefficients[curntb].Coefficient));
                    curnta++;
                    curntb++;
                }
            }

            if (curnta == a.coefficients.Count)
            {
                while (curntb < b.coefficients.Count)
                {
                    c.Add(Monom.Negate(b.coefficients[curntb]));
                    curntb++;
                }
            }
            else if (curntb == b.coefficients.Count)
            {
                while (curnta < a.coefficients.Count)
                {
                    c.Add(a.coefficients[curnta]);
                    curnta++;
                }
            }
            return new Polynom(c, a.Fp);
        }
        /// <summary>
        /// Сумма полиномов
        /// </summary>
        /// <param name="a">Слагаемое</param>
        /// <param name="b">Слагаемое</param>
        /// <returns></returns>
        public static Polynom operator +(Polynom a, Polynom b)
        {
            a.NormPolynom();
            b.NormPolynom();

            List<Monom> c = new List<Monom>();
            int curnta = 0; // Текущий элемент в полиноме а
            int curntb = 0; // Текущий элемент в полиноме b

            while (curnta < a.coefficients.Count && curntb < b.coefficients.Count)
            {
                if (a.coefficients[curnta].Degree < b.coefficients[curntb].Degree)
                {
                    c.Add(b.coefficients[curntb]);
                    curntb++;
                }
                else if (a.coefficients[curnta].Degree > b.coefficients[curntb].Degree)
                {
                    c.Add(a.coefficients[curnta]);
                    curnta++;
                }
                else
                {
                    c.Add(a.coefficients[curnta].AddNumber(b.coefficients[curntb].Coefficient));
                    curnta++;
                    curntb++;
                }
            }

            if (curnta == a.coefficients.Count)
            {
                while (curntb < b.coefficients.Count)
                {
                    c.Add(b.coefficients[curntb]);
                    curntb++;
                }
            }
            else if (curntb == b.coefficients.Count)
            {
                while (curnta < a.coefficients.Count)
                {
                    c.Add(a.coefficients[curnta]);
                    curnta++;
                }
            }
            return new Polynom(c, a.Fp);
        }
        /// <summary>
        /// Умножение полинома на число
        /// </summary>
        /// <param name="number">Множитель</param>
        public void MultiplyInt(BigInteger number)
        {
            for (int i = 0; i < coefficients.Count; i++)
            {
                coefficients[i] = coefficients[i].MultiplyInt(number);
            }
            NormPolynom();
        }
        /// <summary>
        /// Умножение полиномов
        /// </summary>
        /// <param name="A">Множитель</param>
        /// <param name="B">Множитель</param>
        /// <returns></returns>
        public static Polynom operator *(Polynom A, Polynom B)
        {
            List<Monom> c = new List<Monom>();
            for (int i = 0; i < A.coefficients.Count; i++)
            {
                for (int j = 0; j < B.coefficients.Count; j++)
                {
                    c.Add(A.coefficients[i] * B.coefficients[j]);
                }
            }
            return new Polynom(c, A.Fp);
        }
        /// <summary>
        /// Полином в степени
        /// </summary>
        /// <param name="A">Полином</param>
        /// <param name="exp">Степень</param>
        /// <returns></returns>
        public static Polynom ModPow(Polynom A, BigInteger exp, Polynom mod)
        {           
            Polynom t = A;
            // exp -> 100101010101......
            StringBuilder sb = new StringBuilder();
            while(exp.CompareTo(0)==1)
            {
                sb.Append(exp & 0x1);
                exp = exp >> 1;
            }
            var expstr = sb.ToString();
            //------------------------------------
            for(int i= expstr.Length-2; i>=0;i--)
            {
                t = t * t;
                t = Remainder(t, mod);
                if(expstr[i] == '1')
                {
                    t = t * A;
                    t = Remainder(t, mod);
                }
            }
            return t;
        }
        /// <summary>
        /// Остаток от деления.
        /// </summary>
        /// <param name="a">Делимое</param>
        /// <param name="b">Делитель</param>
        /// <returns>Остаток от деления</returns>
        public static Polynom Remainder(Polynom a, Polynom b)
        {
            if (a.Degree < b.Degree) return a;
            List<Monom> div = new List<Monom>(); // целая часть
            while (a.Degree >= b.Degree)
            {
                BigInteger current_deg = a.Degree - b.Degree;
                Monom c = new Monom(current_deg, 1); // целая часть от деления

                BigInteger currentInt = a.coefficients[0].Coefficient * Extension.Inverse(b.coefficients[0].Coefficient, b.Fp);
                c.MultiplyInt(currentInt);
                div.Add(c);
                Polynom temp = new Polynom(new List<Monom> {c}, a.Fp);
                a -= temp * b;
                a.NormPolynom();
            }
            return a;
        }
        /// <summary>
        /// Остаток от деления. С возвратом целой части
        /// </summary>
        /// <param name="a">Делимое</param>
        /// <param name="b">Делитель</param>
        /// <param name="_div">Целая цасть</param>
        /// <returns>Остаток</returns>
        public static Polynom Remainder(Polynom a, Polynom b, out Polynom _div)
        {
            if (a.Degree < b.Degree)
            {
                _div = null;
                return a;
            }
            List<Monom> div = new List<Monom>(); // целая часть
            while (a.Degree >= b.Degree)
            {
                BigInteger current_deg = a.Degree - b.Degree;
                Monom c = new Monom(current_deg, 1); // целая часть от деления

                BigInteger currentInt = a.coefficients[0].Coefficient * Extension.Inverse(b.coefficients[0].Coefficient, b.Fp);
                c.MultiplyInt(currentInt);
                div.Add(c);
                Polynom temp = new Polynom(new List<Monom> { c }, a.Fp);
                a -= temp * b;
                a.NormPolynom();
            }
            _div = new Polynom(div, a.Fp);
            return a;
        }
        /// <summary>
        /// Наибольший общий делитель
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static Polynom GreatCommonDivisor(Polynom a, Polynom b)
        {
            Polynom u = null;
            Polynom v = null;
            if (a.Degree >= b.Degree)
            {
                u = a;
                v = b;
            }
            else
            {
                u = b;
                v = a;
            }
            while(v.Degree!=-1) // пока v>0
            {
                Polynom utemp = u;
                u = v;               
                v = Polynom.QuickRemainder(utemp, v);
            }
            return u;
        }
        /// <summary>
        /// Корни многочлена 
        /// </summary>
        /// <param name="a">Многочлен</param>
        /// <returns>List с корнями</returns>
        public List<BigInteger> GetRoots()
        {
            roots = new List<BigInteger>();
            if (Degree <= 2)
            {
                Polynom x = new Polynom(coefficients, Fp);
                Roots(x);
            }
            else
            {
                roots = new List<BigInteger>();
                //-------------Многочлен x^p - x--------------------
                Polynom t = Parse(string.Format("+x^{0}-x", this.Fp), this.Fp);
                //----------------------------------------------------
                Polynom gcd = GreatCommonDivisor(t, this);
                if (gcd.Degree == 0) return null;
                if (gcd.GetValue(0) == 0)
                {
                    roots.Add(0);
                    Polynom rem = Polynom.Remainder(gcd, Parse("+x", this.Fp), out gcd);
                    if (rem.Degree != -1) throw new Exception("Остаток не равен нулю!");
                }
                Roots(gcd);
            }
            return roots;
        }
        void Roots(Polynom g)
        {            
            if(g.Degree <=2)
            {
                // Решение уравнений вида ax^2+bx+c = 0 mod p a,b,c>=0
                if(g.coefficients[0].Degree==1) // bx + c = 0
                {
                    if (g.coefficients.Count == 1) // bx = 0
                    {
                        roots.Add(0);
                        return;
                    }
                    else // bx + c = 0
                    {
                        BigInteger b = g.coefficients[0].Coefficient;
                        BigInteger c = g.coefficients[1].Coefficient;

                        BigInteger x = BigInteger.Remainder(BigInteger.Negate(c) * Extension.Inverse(b, g.Fp), g.Fp);
                        while (x < 0) x += g.Fp;
                        roots.Add(x);
                        return;
                    }
                    // ошибка
                }
                else // ax^2 + bx + c = 0
                {
                    if(g.coefficients[1].Degree == 0) // ax^2 + c = 0 
                    {
                        BigInteger a = g.coefficients[0].Coefficient;
                        BigInteger c = g.coefficients[1].Coefficient;
                        BigInteger ac = BigInteger.Remainder(BigInteger.Negate(c) * Extension.Inverse(a,g.Fp),g.Fp);
                        BigInteger r = Extension.SquareRootModPrime(ac, g.Fp);
                        roots.Add(r); // не проверенно
                        roots.Add(g.Fp - r);
                        return;
                    }
                    else // ax^2 + bx = 0 && ax^2 = 0 && ax^2 + bx + c = 0
                    {
                        if(g.coefficients.Count == 1)// ax ^ 2 = 0
                        {
                            roots.Add(0);
                            return;
                        }
                        else // ax^2 + bx = 0 && ax^2 + bx + c = 0
                        {
                            if (g.coefficients.Count == 2) // ax^2 + bx = 0
                            {
                                roots.Add(0);
                                Polynom div = null;
                                Polynom.Remainder(g, Parse("+x", this.Fp), out div);
                                Roots(div);
                            }
                            else // ax^2 + bx + c = 0
                            {
                                BigInteger a = g.coefficients[0].Coefficient;
                                BigInteger b = g.coefficients[1].Coefficient;
                                BigInteger c = g.coefficients[2].Coefficient;

                                BigInteger D = Extension.SquareRootModPrime(BigInteger.Remainder(b * b - 4 * a * c,g.Fp), g.Fp);

                                BigInteger x1 = BigInteger.Remainder((BigInteger.Negate(b) - D) * BigInteger.Remainder(Extension.Inverse(2 * a, g.Fp),g.Fp),g.Fp);
                                BigInteger x2 = BigInteger.Remainder((BigInteger.Negate(b) + D) * BigInteger.Remainder(Extension.Inverse(2 * a, g.Fp),g.Fp),g.Fp);
                                while (x1 < 0) x1 += g.Fp;
                                while (x2 < 0) x2 += g.Fp;
                                roots.Add(x1);
                                roots.Add(x2);
                            }
                        }
                    }
                }             
            }
            else
            {
                Polynom one = Parse("+1", g.Fp);
                Polynom h = one;
                while (h.Degree == 0 || h.Eguals(g) )
                {
                    BigInteger a = Extension.Random(1, g.Fp);
                    Polynom s = ModPow(Parse(string.Format("+x-{0}", a), g.Fp), (g.Fp - 1) / 2, g);
                    h = Polynom.GreatCommonDivisor(s-one, g);
                }
                Polynom div = null;
                Remainder(g, h, out div);
                Roots(h);
                Roots(div);                          
            }
        }
        /// <summary>
        /// Значение полинома в точке x
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public BigInteger GetValue(BigInteger x)
        {
            BigInteger t = 0;
            foreach(Monom m in coefficients)
            {
                t += m.Coefficient * BigInteger.ModPow(x, m.Degree, Fp);
            }
            return BigInteger.Remainder(t,Fp);
        }
        #endregion

        #region Побочные методы
        /// <summary>
        /// Нормировка полинома в заданном поле.
        /// </summary>
        private void NormPolynom()
        {
                for (int i = 0; i < coefficients.Count; i++)
                {
                    if (coefficients[i].Coefficient < 0)
                    {
                        BigInteger k = BigInteger.Abs(coefficients[i].Coefficient) / Fp + 1;
                        coefficients[i] = coefficients[i].AddNumber(Fp * k);
                    }
                    if (coefficients[i].Coefficient >= Fp)
                    {
                        coefficients[i] = coefficients[i].Remainder(Fp);
                    }
                }
        }
        /// <summary>
        /// Строковое представление
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            CheckPoly();

            if (Degree == 0 || Degree == -1)// Если нулевая степень полинома или нулевой
            {
                sb.Append(coefficients[0].Coefficient);
            }
            else
            {
                sb.Append(string.Format("{0}x^{1}", coefficients[0].Coefficient, coefficients[0].Degree)); // первый моном

                for (int i = 1; i < coefficients.Count - 1; i++) // от второго до предпоследнего монома
                {
                    if (coefficients[i].Coefficient > 0)
                    {
                        sb.Append(string.Format(" + {0}x^{1}", coefficients[i].Coefficient, coefficients[i].Degree));
                    }
                    else
                    {
                        sb.Append(string.Format(" - {0}x^{1}", BigInteger.Abs(coefficients[i].Coefficient), coefficients[i].Degree));
                    }
                }

                if (coefficients[coefficients.Count - 1].Degree == 0) // Если степень последнего нулевая (Свободный член)
                {
                    if (coefficients[coefficients.Count - 1].Coefficient > 0)
                    {
                        sb.Append(string.Format(" + {0}", coefficients[coefficients.Count - 1].Coefficient));
                    }
                    else
                    {
                        sb.Append(string.Format(" - {0}", BigInteger.Abs(coefficients[coefficients.Count - 1].Coefficient)));
                    }
                }
                else
                {
                    if (coefficients[coefficients.Count - 1].Coefficient > 0)
                    {
                        sb.Append(string.Format(" + {0}x^{1}", coefficients[coefficients.Count - 1].Coefficient, coefficients[coefficients.Count - 1].Degree));
                    }
                    else
                    {
                        sb.Append(string.Format(" - {0}x^{1}", BigInteger.Abs(coefficients[coefficients.Count - 1].Coefficient), coefficients[coefficients.Count - 1].Degree));
                    }
                }
            }
            return sb.ToString();
        }
        /// <summary>
        /// Создает объект типа Polynom из переданной строки
        /// </summary>
        /// <remarks>Перед первым членом обязательно указывается знак. Например: Polynom.Parse("+x^3+2x^2+3x^1+4,19) </remarks>
        /// <param name="_poly">Строка с полиномом</param>
        /// <param name="_Fp">Модуль(Поле Галуа)</param>
        public static Polynom Parse(string _poly, BigInteger _Fp)
        {
            if (_poly == "") throw new Exception("string is Empty"); // пустая строка

            StringBuilder poly = new StringBuilder(_poly);
            List<Monom> polycoeff = new List<Monom>(); // Коэфф. полинома в разреженном виде

            while(poly.Length!=0)
            {
                BigInteger coef = 1;
                BigInteger deg = 1;
                #region sign
                bool sign = true;
                if(poly[0]!= '-' && poly[0]!='+') throw new Exception("Input string was malformed");
                if (poly[0] =='-')
                {
                    sign = false;
                }
                #endregion

                if (poly[1] == 'x')
                {
                    if(poly.Length >2 && poly[2] =='^')
                    {
                        int countdeg = ParseNumber(poly, 3);
                        deg = BigInteger.Parse(poly.ToString(3, countdeg));
                        poly = poly.Remove(0, countdeg + 3);
                        if (sign)
                        {
                            polycoeff.Add(new Monom(deg, coef));
                        }
                        else
                        {
                            polycoeff.Add(new Monom(deg, -coef));
                        }
                    }
                    else
                    {
                        poly = poly.Remove(0, 2);
                        if (sign)
                        {
                            polycoeff.Add(new Monom(deg, coef));
                        }
                        else
                        {
                            polycoeff.Add(new Monom(deg, -coef));
                        }
                    }
                }
                else if (poly[1] >= '1' && poly[1] <= '9')
                {
                    int count = ParseNumber(poly, 1);
                    if (count + 1 < poly.Length) 
                    {
                        string c = poly.ToString(1, count);
                        if(sign)
                        {
                            coef = BigInteger.Parse(c);
                        }
                        else
                        {
                            coef = BigInteger.Negate(BigInteger.Parse(c));
                        }
                
                        if (poly[count+1] != 'x') throw new Exception("Input string was malformed");

                        if(count+2 < poly.Length)
                        {
                            if (poly[count + 2] == '^')
                            {
                                int countdeg = ParseNumber(poly, count + 3);
                                string d = poly.ToString(count + 3, countdeg);
                                deg = BigInteger.Parse(d);

                                polycoeff.Add(new Monom(deg, coef));
                                poly = poly.Remove(0, count + countdeg + 3);
                            }
                            else if (poly[count + 2] == '+' || poly[count + 2] == '-')
                            {
                                polycoeff.Add(new Monom(deg, coef));
                                poly = poly.Remove(0, count + 2);
                            }
                            else throw new Exception("Input string was malformed");

                        }
                        else
                        {
                            polycoeff.Add(new Monom(deg, coef));
                            poly = poly.Remove(0, count + 2);
                        }                                                
                    }
                    else // свободный член
                    {
                        string c = poly.ToString(1, count);
                        coef = BigInteger.Parse(c);
                        if(sign)
                         polycoeff.Add(new Monom(0, coef));
                        else
                            polycoeff.Add(new Monom(0, -coef));
                        poly = poly.Remove(0, count+1);
                    }                   
                }
                else throw new Exception("Input string was malformed");    
                          
            }
            return new Polynom(polycoeff,_Fp);
        }
        /// <summary>
        /// Выделяет целое число из строки с полиномом
        /// </summary>
        /// <param name="_poly">Полином в строковой записи</param>
        /// <param name="startindex">Начальный индекс числа</param>
        /// <returns></returns>
        static int ParseNumber(StringBuilder _poly, int startindex)
        {
            if (_poly[startindex] < '1' || _poly[startindex] > '9') return 0;
            int count = 1;
            while (count+ startindex < _poly.Length && (_poly[startindex + count] >= '0' && _poly[startindex + count] <= '9'))
            {
                count++;
            }
            return count;
        }
        /// <summary>
        /// Удаление нулевых коэффициентов
        /// </summary>
        void DeleteZeroCoef()
        {
            coefficients.RemoveAll(delegate (Monom a)
            {
                return a.Coefficient.CompareTo(0) == 0;
            });
        } 
        /// <summary>
        /// Сортировка по убыванию степеней
        /// </summary>
        void SortOfDegree()
        {
            coefficients.Sort(delegate (Monom a1, Monom a2)
            {
                return a2.Degree.CompareTo(a1.Degree);
            });
        }
        /// <summary>
        /// Суммирует мономы с одинаковыми степенями в пределах одного полинома
        /// </summary>
        void SummEqualsDegry()
        {
            SortOfDegree();
            for(int i = 0; i < coefficients.Count - 1;)
            {
                if (coefficients[i].Degree == coefficients[i + 1].Degree)
                {
                    coefficients[i] =  coefficients[i].AddNumber(coefficients[i + 1].Coefficient);
                    coefficients.RemoveAt(i + 1);
                }
                else
                    i++;
            }
        }
        /// <summary>
        /// Удаление нулевых коэффициентов, Сортировка,сложение одинаковых степеней, установка степени полинома
        /// </summary>
        void CheckPoly()
        {          
            SummEqualsDegry();
            NormPolynom();
            DeleteZeroCoef();
            if (coefficients.Count == 0)
            {
                Degree = -1;
                coefficients.Add(new Monom(0, 0));
            }
            else
                Degree = coefficients[0].Degree;
        }
        bool Eguals(Polynom b)
        {
            if(Fp == b.Fp) // поле одинаковое?
            {
                if(Degree == b.Degree) // степень одинаковая?
                {
                    if(coefficients.Count == b.coefficients.Count()) // Кол-во мономов одинаковое?
                    {
                        for(int i = 0; i<coefficients.Count(); i++)
                        {
                            if (coefficients[i].Coefficient != b.coefficients[i].Coefficient
                                || coefficients[i].Degree != b.coefficients[i].Degree) return false;
                        }
                    }
                    else
                    {
                        return false;
                    }
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
            return true;
        }

        #endregion


        #region Быстрое деление Многочленов
        /// <summary>
        /// Укорачивание до указанной степени 
        /// </summary>
        /// <param name="a">Полином</param>
        /// <param name="deg">Степень укорачивания</param>
        /// <returns>Укороченный полином степени меньшей deg</returns>
        static Polynom Shortening(Polynom a,BigInteger deg)
        {
            List<Monom> coef = new List<Monom>();
            foreach(var m in a.coefficients)
            {
                coef.Add(m);
            }             
            for(int i = coef.Count-1; i >= 0; i--)
            {
                if(coef[i].Degree >= deg)
                {
                    coef.RemoveRange(0, i + 1);
                    break;
                }
            }
            return new Polynom(coef, a.Fp);
        }
        /// <summary>
        /// Быстрое обращение многочлена
        /// </summary>
        /// <param name="x">Многочлен</param>
        /// <param name="deg">Степень укорачивания</param>
        /// <returns></returns>
        static Polynom R(Polynom x, BigInteger deg)
        {
            Polynom g = Polynom.Parse("+1", x.Fp);
            Polynom two = Polynom.Parse("+2", x.Fp);
            BigInteger n = 1;
            while(n < deg+1)
            {
                n = n << 1;
                if (n > deg + 1) n = deg + 1;
                Polynom h = Shortening(x, n);
                h = Shortening(h * g, n);
                g = Shortening(g * (two - h), n);
            }
            return g;
        }
        /// <summary>
        /// Померанц 575 страница
        /// </summary>
        /// <param name="x"></param>
        /// <param name="deg"></param>
        /// <returns></returns>
        static Polynom rev(Polynom x, BigInteger deg)
        {
            //Формируем список со степенями полинома x
            List<BigInteger> degree = new List<BigInteger>();
            foreach(var a in x.coefficients)
            {
                degree.Add(a.Degree);
            }
            List<Monom> temp = new List<Monom>();
            int i = 0;
            while(i<degree.Count)
            {
                BigInteger d = degree[i];
                BigInteger j = deg - d;
                if (j < 0)
                {
                    i++;
                    continue;
                }
                BigInteger cf = SearchCoefWithDeg(x, deg - j);
                if (cf != 0)
                {
                    temp.Add(new Monom(j, cf));
                }
                i++;
            }            
            return new Polynom(temp, x.Fp);
        }
        /// <summary>
        /// Померанц 575 страница
        /// </summary>
        /// <param name="x"></param>
        /// <param name="deg"></param>
        /// <returns></returns>
        static BigInteger ind(Polynom x, BigInteger deg)
        {
            int i = x.coefficients.Count() - 1;
            while (x.coefficients[i].Degree < deg)
            {
                i--;
                if (i < 0) return 0;
            }
            return x.coefficients[i].Degree;
        }
        /// <summary>
        /// Коэффициент при степени
        /// </summary>
        /// <param name="x">Полином</param>
        /// <param name="deg">Необходимая степень</param>
        /// <returns></returns>
        static BigInteger SearchCoefWithDeg(Polynom x, BigInteger deg)
        {
            int i = 0;
            while(x.coefficients[i].Degree>deg)
            {
                i++;
            }
            if (x.coefficients[i].Degree == deg)
                return x.coefficients[i].Coefficient;
            else
                return 0;
        }


        public static Polynom QuickRemainder(Polynom x, Polynom y)
        {
            if (y.Degree == 0) return new Polynom(new List<Monom> { }, x.Fp);
            BigInteger d = x.Degree - y.Degree;
            if (d < 0) return x;

            Polynom X = rev(x, x.Degree);
            Polynom Y = rev(y, y.Degree);

            Polynom q = R(Y, d);

            q = Shortening(q * X, d + 1);
            Polynom r = X - q * Y;
            BigInteger i = ind(r, d + 1);
            Remainder(r, Parse(string.Format("+x^{0}", i), x.Fp), out r);
            return rev(r, x.Degree - i);
        }
        #endregion


    }
    struct Monom
    {
        public BigInteger Degree
        {
            get;
            private set;
        }
        public BigInteger Coefficient
        {
            get;
            private set;
        }
        public Monom(BigInteger _degre, BigInteger _coeff)
        {
            Degree = _degre;
            Coefficient = _coeff;
        }
        /// <summary>
        /// Прибавляет к коэффициенту число
        /// </summary>
        /// <param name="number">Число</param>
        public Monom AddNumber(BigInteger number)
        {
            Coefficient += number;
            return this;
        }
        /// <summary>
        /// Отнимает от коэффициента число
        /// </summary>
        /// <param name="number">Число</param>
        public Monom SubNumber(BigInteger number)
        {
            Coefficient -= number;
            return this;
        }
        /// <summary>
        /// Остаток по модулю
        /// </summary>
        /// <param name="number">Модуль</param>
        public Monom Remainder(BigInteger number)
        {
            Coefficient = BigInteger.Remainder(Coefficient,number);
            return this;
        }
        /// <summary>
        /// Изменяет знак коэффициента
        /// </summary>
        /// <returns></returns>
        static public Monom Negate(Monom a)
        {
            return new Monom(a.Degree, BigInteger.Negate(a.Coefficient));
        }
        /// <summary>
        /// Умножение монома на число
        /// </summary>
        /// <param name="number">Множитель</param>
        /// <returns></returns>
        public Monom MultiplyInt(BigInteger number)
        {
            Coefficient = BigInteger.Multiply(Coefficient, number);
            return this;
        }
        /// <summary>
        /// Умножение монома на моном
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static Monom operator *(Monom a,Monom b)
        {
            Monom c = new Monom();
            c.Coefficient = a.Coefficient * b.Coefficient;
            c.Degree = a.Degree + b.Degree;
            return c;
        }
    }
}
