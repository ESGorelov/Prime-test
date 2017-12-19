using System;
using System.Numerics;
using System.Collections.Generic;


namespace Atkin_MoreinPrimeTest
{
    class ComplexPolynom
    {
        #region свойства
        public BigInteger Degree
        {
            get;
            private set;
        }  
        List<ComplexMonom> coefficients;
        #endregion
        #region Конструкторы
        public ComplexPolynom(List<ComplexMonom> coeff)
        {
            coefficients = coeff;
            CheckPoly();
        }

        #endregion

        #region Математические Операции с Полиномами       
        public static ComplexPolynom MultiplyPol(ComplexPolynom poly, Complex j, string type = "line")
        {
            if (type == "line") // (x-j)
            {
                ComplexPolynom t = new ComplexPolynom(new List<ComplexMonom> { new ComplexMonom(1, 1), new ComplexMonom(0, -j) });
                t = t * poly;
                return t;
            }
            else if (type == "quadr") // X^2 - 2Re(j)X + |j|^2
            {
                ComplexPolynom t = new ComplexPolynom(new List<ComplexMonom> { new ComplexMonom(2, 1), new ComplexMonom(1, 2*j.Real),new ComplexMonom(0,j.Magnitude*j.Magnitude) });
                t = t * poly;
                return t;
            }
            else return null;
        }

        /// <summary>
        /// Умножение полинома на число
        /// </summary>
        /// <param name="number">Множитель</param>
        public void MultiplyInt(Complex number)
        {
            for (int i = 0; i < coefficients.Count; i++)
            {
                coefficients[i] = coefficients[i].MultiplyInt(number);
            }
        }
        /// <summary>
        /// Умножение полиномов
        /// </summary>
        /// <param name="A">Множитель</param>
        /// <param name="B">Множитель</param>
        /// <returns></returns>
        public static ComplexPolynom operator *(ComplexPolynom A, ComplexPolynom B)
        {
            List<ComplexMonom> c = new List<ComplexMonom>();
            for (int i = 0; i < A.coefficients.Count; i++)
            {
                for (int j = 0; j < B.coefficients.Count; j++)
                {
                    c.Add(A.coefficients[i] * B.coefficients[j]);
                }
            }
            return new ComplexPolynom(c);
        }
        #endregion
        #region Побочные методы                    
        /// <summary>
        /// Удаление нулевых коэффициентов
        /// </summary>
        void DeleteZeroCoef()
        {
            coefficients.RemoveAll(delegate (ComplexMonom a)
            {
                return a.Coefficient.Equals(new Complex(0,0));
            });
        } 
        /// <summary>
        /// Сортировка по убыванию степеней
        /// </summary>
        void SortOfDegree()
        {
            coefficients.Sort(delegate (ComplexMonom a1, ComplexMonom a2)
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
            DeleteZeroCoef();
            if (coefficients.Count == 0)
            {
                Degree = -1;
                coefficients.Add(new ComplexMonom(0, 0));
            }
            else
                Degree = coefficients[0].Degree;
        }
        /// <summary>
        /// Получение полинома в кольце целых чисел
        /// </summary>
        /// <returns></returns>
        public Polynom GetPolynomReal(BigInteger Fp)
        {
            List<Monom> t = new List<Monom>();
            foreach (var a in coefficients)
            {
                t.Add(new Monom(a.Degree, new BigInteger(Math.Round(a.Coefficient.Real))));
            }
            return new Polynom(t, Fp);
        }
        #endregion
    }
    struct ComplexMonom
    {
        public BigInteger Degree
        {
            get;
            private set;
        }
        public Complex Coefficient
        {
            get;
            private set;
        }
        public ComplexMonom(BigInteger _degre, Complex _coeff)
        {
            Degree = _degre;
            Coefficient = _coeff;
        }
        /// <summary>
        /// Прибавляет к коэффициенту число
        /// </summary>
        /// <param name="number">Число</param>
        public ComplexMonom AddNumber(Complex number)
        {
            Coefficient += number;
            return this;
        }
        /// <summary>
        /// Отнимает от коэффициента число
        /// </summary>
        /// <param name="number">Число</param>
        public ComplexMonom SubNumber(Complex number)
        {
            Coefficient -= number;
            return this;
        }
        /// <summary>
        /// Изменяет знак коэффициента
        /// </summary>
        /// <returns></returns>
        static public ComplexMonom Negate(ComplexMonom a)
        {
            return new ComplexMonom(a.Degree, Complex.Negate(a.Coefficient));
        }
        /// <summary>
        /// Умножение монома на число
        /// </summary>
        /// <param name="number">Множитель</param>
        /// <returns></returns>
        public ComplexMonom MultiplyInt(Complex number)
        {
            Coefficient = Complex.Multiply(Coefficient, number);
            return this;
        }
        /// <summary>
        /// Умножение монома на моном
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static ComplexMonom operator *(ComplexMonom a, ComplexMonom b)
        {
            ComplexMonom c = new ComplexMonom();
            c.Coefficient = a.Coefficient * b.Coefficient;
            c.Degree = a.Degree + b.Degree;
            return c;
        }
    }
}
