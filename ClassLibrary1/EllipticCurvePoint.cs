using System.Numerics;


namespace Atkin_MoreinPrimeTest
{
    class EllipticCurvePoint
    {
        #region Поля

        private BigInteger p; // GF(p)
        public BigInteger P
        {
            get { return p; }
        }
        private BigInteger a; // коэфф прямой
        public BigInteger A
        {
            get { return a; }
        }
        private BigInteger b; // коэфф прямой
        public BigInteger B
        {
            get { return b; }
        }
        private BigInteger x; // координата х 
        public BigInteger X
        {
            get { return x; }
            set { x = value; }
        }
        private BigInteger y; // координата y 
        public BigInteger Y
        {
            get { return y; }
            set { y = value; }
        }
        private BigInteger z; // координата z 
        public BigInteger Z
        {
            get { return z; }
            set { z = value; }
        }
        #endregion

        #region Building
        /// <summary>
        /// Создание эллиптической кривой с заданными параметрами и случайной точкой на ней.
        /// </summary>
        /// <param name="P">Модуль</param>
        /// <param name="A">Параметр А</param>
        /// <param name="B">Параметр B</param>
        public EllipticCurvePoint(BigInteger P, BigInteger A, BigInteger B)
        {
            p = P;
            a = A;
            b = B;
            var temp = GetPoint();
            x = temp[0];
            y = temp[1];
            z = 1;
        }
        /// <summary>
        /// Создание кривой с пользовательскими параметрами
        /// </summary>
        /// <param name="P">Модуль</param>
        /// <param name="A">Коэф. А</param>
        /// <param name="B">Коэф. B</param>
        /// <param name="X">Координата X</param>
        /// <param name="Y">Координата Y</param>
        /// <param name="Z">Координата Z</param>
        EllipticCurvePoint(BigInteger P, BigInteger A, BigInteger B, BigInteger X, BigInteger Y, BigInteger Z)
        {
            p = P;
            a = A;
            b = B;
            x = X;
            y = Y;
            z = Z;
        }
        EllipticCurvePoint(EllipticCurvePoint otherPoint)
        {
            p = otherPoint.P;
            a = otherPoint.A;
            b = otherPoint.B;
            x = otherPoint.X;
            y = otherPoint.Y;
            z = otherPoint.Z;
        }
        #endregion

        #region Сложение точек
        public static EllipticCurvePoint AddPoint(EllipticCurvePoint pointA, EllipticCurvePoint pointB)
        {
            if ((pointA.X == pointA.Z) && (pointA.X == 0)) { return pointB; } 
            if ((pointB.X == pointB.Z) && (pointB.X == 0)) { return pointA; }       

            BigInteger A = pointB.Y * pointA.Z - pointA.Y * pointB.Z;
            A = BigInteger.Remainder(A, pointA.P);
            if (A < 0) { while (A < 0) { A += pointA.P; } }

            BigInteger B = pointB.X * pointA.Z - pointA.X * pointB.Z;
            B = BigInteger.Remainder(B, pointA.P);
            if (B < 0) { while (B < 0) { B += pointA.P; } }

            BigInteger C = pointB.X * pointA.Z + pointA.X * pointB.Z;
            C = BigInteger.Remainder(C, pointA.P);
            if (C < 0) { while (C < 0) { C += pointA.P; } }

            BigInteger D = pointB.X * pointA.Z + 2 * pointA.X * pointB.Z;
            D = BigInteger.Remainder(D, pointA.P);
            if (D < 0) { while (D < 0) { D += pointA.P; } }

            BigInteger E = pointB.Z * pointA.Z;
            E = BigInteger.Remainder(E, pointA.P);
            if (E < 0) { while (E < 0) { E += pointA.P; } }

            BigInteger X = B * (E * A * A - C * B * B);
            X = BigInteger.Remainder(X, pointA.P);
            if (X < 0) { while (X < 0) { X += pointA.P; } }

            BigInteger Y = A * (B * B * D - A * A * E) - pointA.Y * pointB.Z * B * B * B;
            Y = BigInteger.Remainder(Y, pointA.P);
            if (Y < 0) { while (Y < 0) { Y += pointA.P; } }

            BigInteger Z = B * B * B * E;
            Z = BigInteger.Remainder(Z, pointA.P);
            if (Z < 0) { while (Z < 0) { Z += pointA.P; } }


            if (X == 0 && Z == 0)
            {
                return new EllipticCurvePoint(pointA.P, pointA.A, pointA.B, 0, 1, 0);
            }

            return new EllipticCurvePoint(pointA.P, pointA.A, pointA.B, X, Y, Z);
        }
        private static EllipticCurvePoint X2(EllipticCurvePoint point)
        {
            BigInteger A = 3 * point.X * point.X + point.A * point.Z * point.Z;
            A = BigInteger.Remainder(A, point.P);
            if (A < 0) { while (A < 0) { A += point.P; } }


            BigInteger B = 2 * point.Y * point.Z;
            B = BigInteger.Remainder(B, point.P);
            if (B < 0) { while (B < 0) { B += point.P; } }


            BigInteger X = B * (A * A - 4 * point.X * point.Y * B);
            X = BigInteger.Remainder(X, point.P);
            if (X < 0) { while (X < 0) { X += point.P; } }


            BigInteger Y = A * (6 * point.Y * point.X * B - A * A) - 2 * point.Y * point.Y * B * B;
            Y = BigInteger.Remainder(Y, point.P);
            if (Y < 0) { while (Y < 0) { Y += point.P; } }


            BigInteger Z = B * B * B;
            Z = BigInteger.Remainder(Z, point.P);
            if (Z < 0) { while (Z < 0) { Z += point.P; } }
            return new EllipticCurvePoint(point.P, point.A, point.B, X, Y, Z);
        }
        public static EllipticCurvePoint Multiply(BigInteger k, EllipticCurvePoint point)
        {
            EllipticCurvePoint temp = point;
            k--;
            while (k != 0)
            {
                if (BigInteger.Remainder(k, 2) != 0)
                {
                    if ((temp.X == point.X) && (temp.Y == point.Y))
                        temp = X2(temp);
                    else
                        temp = AddPoint(temp, point);
                    k--;
                }
                k = k / 2;
                point = X2(point);
            }
            return temp;
        }
        #endregion

        
              
        /// <summary>
        /// Привидение точки к афинным координатам
        /// </summary>
        /// <param name="p">Точка на кривой</param>
        /// <returns></returns>
        public static EllipticCurvePoint AffineCoords(EllipticCurvePoint p)
        {
            if (p.Z == 0) return p;
            p.X = BigInteger.Remainder(p.X * Extension.Inverse(p.Z, p.P), p.P);
            p.Y = BigInteger.Remainder(p.Y * Extension.Inverse(p.Z, p.P), p.P);
            p.Z = 1;
            return p;
            
        }
        /// <summary>
        /// Случайная точка на Данной кривой.
        /// </summary>    
        /// <returns></returns>
        BigInteger[] GetPoint()
        {
            BigInteger rnd = Extension.Random(1, P);
            while (Extension.Lezhandr(BigInteger.Remainder((rnd * rnd * rnd + A * rnd + B), P), P) != 1)
            {
                rnd = Extension.Random(1, P);
            }
            BigInteger xcoord = BigInteger.Remainder((rnd * rnd * rnd + A * rnd + B), P); // Правая часть уравнения кривой
            return new BigInteger[] {rnd, Extension.SquareRootModPrime(xcoord,P) };
        }
        /// <summary>
        /// Проверка принадлежности точки кривой
        /// </summary>
        /// <returns></returns>
        public bool CheckPoint()
        {
            if (this.X == 0 && this.Z == 0)
                return true;
            EllipticCurvePoint z = new EllipticCurvePoint(this);
            z = AffineCoords(z);
            BigInteger left = BigInteger.Remainder(z.Y * z.Y, z.P);
            BigInteger right = BigInteger.Remainder((z.X * z.X * z.X + z.A * z.X + z.B), z.P);
            return right.Equals(left);
        }
        /// <summary>
        /// Следующая случайная точка на кривой
        /// </summary>
        public void NextPoint()
        {
            var point = this.GetPoint();
            X = point[0];
            Y = point[1];
            Z = 1;
        }
        public override string ToString()
        {
            string temp = "X= " + X +" Y= " + Y + " Z= " + Z;
            return temp;
        }      
    }
}
