using MathNet.Numerics.LinearAlgebra;

using Vector = MathNet.Numerics.LinearAlgebra.Vector<double>;
using Matrix = MathNet.Numerics.LinearAlgebra.Matrix<double>;

namespace FEM2D
{
    internal class FEM
    {
        public Model m { get; set; }

        private Matrix KK;
        private Vector F;

        public FEM(Model m)
        {
            this.m = m;
            KK = Matrix.Build.Dense(m.n, m.n, 0);
            F = Vector.Build.Dense(m.n, 0);
        }

        internal void AssembleGlobalStiffness()
        {
            for (int i = 0; i < m.nElements; i++)
            {
                // Local stiffness matrix of the element i
                var kl = TrussKl(i);
                // Transformation matrix for this element to global coordinate system
                var C = TrussC(i);

                // Transform element stiffness matrix from local to global coordniate system
                var kg = C * kl * C.Transpose();

                // Element steering vector
                var g = TrussG(i);

                // Add to global matrix
                FormKK(kg, g);
            }
        }

        internal void FormTrussF()
        {
            for (int i = 0; i < m.nNodes; i++)
                for (int j = 0; j < m.DOFPerNode; j++)
                    if (Math.Abs(m.NF[i, j]) > double.Epsilon)
                        F[(int)m.NF[i, j] - 1] = m.Load[i, j];
        }

        /// <summary>
        /// Returns the global stiffness matrix
        /// </summary>
        /// <param name="i"></param>
        /// <returns></returns>
        private void FormKK(Matrix kg, Vector g)
        {
            for (int i = 0; i < m.DOFPerElement; i++)
                if (Math.Abs(g[i]) > double.Epsilon)
                    for (int j = 0; j < m.DOFPerElement; j++)
                        if (Math.Abs(g[j]) > double.Epsilon)
                            KK[(int)g[i] - 1, (int)g[j] - 1] += kg[i, j];
        }

        /// <summary>
        /// This function forms the steering vector for element i
        /// </summary>
        /// <param name="i"></param>
        /// <returns></returns>
        private Vector TrussG(int i)
        {
            var node1 = m.Topology[i, 0];
            var node2 = m.Topology[i, 1];
            var g = Vector.Build.Dense(4);
            g[0] = m.NF[(int)node1, 0];
            g[1] = m.NF[(int)node1, 1];
            g[2] = m.NF[(int)node2, 0];
            g[3] = m.NF[(int)node2, 1];
            return g;
        }

        /// <summary>
        /// This function forms the transformation between local and global coordinates
        /// </summary>
        /// <param name="i"></param>
        /// <returns></returns>
        private Matrix TrussC(int i)
        {
            // Retrieve the nodes of element i
            var node1 = m.Topology[i, 0];
            var node2 = m.Topology[i, 1];

            // Retrieve the x and y coordinates of nodes 1 and 2
            var x1 = m.Geometry[(int)node1, 0]; var y1 = m.Geometry[(int)node1, 1];
            var x2 = m.Geometry[(int)node2, 0]; var y2 = m.Geometry[(int)node2, 1];

            // Evaluate the angle that the member makes with the global axis X
            double theta;
            if (Math.Abs(x2 - x1) < double.Epsilon)
            { // Its orthogonal to it.
                if (y2 > y1)
                {
                    theta = 2 * Math.Atan(1);
                }
                else
                {
                    theta = -2 * Math.Atan(1);
                }
            }
            else
            { // Slope of the line
                theta = Math.Atan((y2 - y1) / (x2 - x1));
            }

            // Construct the transformation matrix
            return Matrix.Build.DenseOfArray(new double[,] {
                {Math.Cos(theta), -Math.Sin(theta), 0, 0 },
                {Math.Sin(theta),  Math.Cos(theta), 0, 0 },
                { 0, 0, Math.Cos(theta),-Math.Sin(theta) },
                { 0, 0, Math.Sin(theta), Math.Cos(theta) }
            });
        }

        /// <summary>
        /// Returns element stiffness matrix KL in local coordinates
        /// </summary>
        /// <param name="i"></param>
        /// <returns></returns>
        private Matrix TrussKl(int i)
        {
            var node1 = m.Topology[i, 0];
            var node2 = m.Topology[i, 1];

            // Retrieve the x and y coordinates of nodes 1 and 2
            var x1 = m.Geometry[(int)node1, 0]; var y1 = m.Geometry[(int)node1, 1];
            var x2 = m.Geometry[(int)node2, 0]; var y2 = m.Geometry[(int)node2, 1];

            // Evaluate length of element i
            double L = Math.Sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));

            // Retrieve section properties of element i
            double E = m.Properties[i, 0];
            double A = m.Properties[i, 1];

            // Calculate element stiffness matrix in its local coordinates
            return Matrix.Build.DenseOfArray(new double[,] {
                { E * A / L, 0, - E * A / L, 0 },
                {         0, 0,           0, 0 },
                {-E * A / L, 0,   E * A / L, 0 },
                {         0, 0,           0, 0 }
            });
        }

        internal Vector Solve()
        {
            return KK.Solve(F);
        }

        internal Matrix ExtractNodalDisplacements(Vector delta)
        {
            var nodeDisp = Matrix.Build.Dense(m.nNodes, m.DOFPerNode, 0);
            for (int i = 0; i < m.nNodes; i++)
            {
                for (int j = 0; j < m.DOFPerNode; j++)
                {
                    nodeDisp[i, j] = 0;
                    if (Math.Abs(m.NF[i, j]) > double.Epsilon)
                    {
                        nodeDisp[i, j] = delta[(int)m.NF[i, j] - 1];
                    }
                }
            }
            return nodeDisp;
        }

        internal Vector BuildLocalForces(Vector delta)
        {
            var force = Vector.Build.Dense(m.nElements);
            // Calculate the forces acting on each element
            // in local coordinates, and store them in the
            // vector force().
            for (int i = 0; i < m.nElements; i++)
            {
                var kl = TrussKl(i);
                var C = TrussC(i);
                var kg = C * kl * C.Transpose();
                var g = TrussG(i);


                var edg = Vector.Build.Dense(m.DOFPerElement);
                for (int j = 0; j < m.DOFPerElement; j++)
                {
                    if (Math.Abs(g[j]) < double.Epsilon)
                    {
                        // Displacement of zero for restrained freedom
                        edg[j] = 0;
                    }
                    else
                    {
                        edg[j] = delta[(int)g[j] - 1];
                    }
                }

                // Element force vector in global XY
                var fg = kg * edg;
                // Element force vector in local xy
                var fl = C.Transpose() * fg;
                force[i] = fl[2];
            }
            return force;
        }
    }
}
