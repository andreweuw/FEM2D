using MathNet.Numerics.LinearAlgebra;

using Matrix = MathNet.Numerics.LinearAlgebra.Matrix<double>;

namespace FEM2D
{
    internal class Model
    {
        internal enum Scenario
        {
            Simple,
            Bridge
        }

        public int nNodes { get; set; }
        public int nElements { get; set; }
        public int nNodesPerElement { get; set; }
        public int DOFPerNode { get; set; }
        public int n { get; set; }

        public int DOFPerElement { get; set; }

        public Matrix Geometry { get; set; }
        public Matrix Topology { get; set; }
        public Matrix Properties { get; set; }
        public Matrix NF { get; set; }
        public Matrix Load { get; set; }

        public Model(Scenario scenario)
        {
            Console.WriteLine("Initializing the model with the " + scenario + " scenario...");
            switch (scenario)
            {
                case Scenario.Simple: InitSimple(); break;
                case Scenario.Bridge: InitBridge(); break;
            }
            
        }

        private void InitBridge()
        {
            nNodes = 9;
            nElements = 15;
            nNodesPerElement = 2;
            DOFPerNode = 2;

            DOFPerElement = nNodesPerElement * DOFPerNode;

            Geometry = Matrix.Build.DenseOfArray(new double[,] {
                { 0, 0 }, // node1
                { 1, 2 }, // node2
                { 2, 0 }, // node3
                { 3, 2 }, // node4
                { 4, 0 }, // node5
                { 5, 2 }, // node6
                { 6, 0 }, // node7
                { 7, 2 }, // node8
                { 8, 0 }  // node9
            });
            Topology = Matrix.Build.DenseOfArray(new double[,] {
                { 0, 1 }, // el1
                { 0, 2 }, // el2
                { 1, 2 }, // el3
                { 1, 3 }, // el4
                { 2, 3 }, // el5
                { 2, 4 }, // el6
                { 3, 4 }, // el7
                { 3, 5 }, // el8
                { 4, 5 }, // el9
                { 4, 6 }, // el10
                { 5, 6 }, // el11
                { 5, 7 }, // el12
                { 6, 7 }, // el13
                { 6, 8 }, // el14
                { 7, 8 }  // el15
            });

            double E = 30e6;   // Yungs modulus
            double A1 = 0.02;   // Cross section of an diagonal edge
            double A2 = 0.045;   // Cross section of an horizontal edge
            // Edge == element, there is a pair of E, A for each element
            Properties = Matrix.Build.DenseOfArray(new double[,] {
                { E, A1 }, // el1   .... diagonals have lower cross section
                { E, A2 }, // el2
                { E, A1 }, // el3
                { E, A2 }, // el4
                { E, A1 }, // el5
                { E, A2 }, // el6
                { E, A1 }, // el7
                { E, A2 }, // el8
                { E, A1 }, // el9
                { E, A2 }, // el10
                { E, A1 }, // el11
                { E, A2 }, // el12
                { E, A1 }, // el13
                { E, A2 }, // el14
                { E, A1 }  // el15
            });
            // Boundary conditions
            // Nodal freedom matrix (number of DOFs at each node)
            NF = Matrix.Build.Dense(nNodes, DOFPerNode, 1);
            NF[0, 0] = 0; NF[0, 1] = 0; // Node 1
            NF[8, 1] = 0; // Node 3

            // Count the free degrees of freedom
            CountDegreesOfFreedom();

            // Forces
            Load = Matrix.Build.Dense(nNodes, DOFPerNode, 0);
            // Forces in X and Y direction at node 3
            Load[1, 0] = 15; // force in X dir at node 2
            Load[2, 1] = -5; // force in Y dir at node 3
            Load[3, 1] = -7; // force in Y dir at node 4
            Load[6, 1] = -10;// force in Y dir at node 7
        }

        private void InitSimple()
        {
            nNodes = 3;
            nElements = 3;
            nNodesPerElement = 2;
            DOFPerNode = 2;

            DOFPerElement = nNodesPerElement * DOFPerNode;

            Geometry = Matrix.Build.DenseOfArray(new double[,] {
                { 0, 0 },
                { 4000, 0 },
                { 4000, 6000 }
            });
            Topology = Matrix.Build.DenseOfArray(new double[,] {
                { 0, 1 },
                { 1, 2 },
                { 0, 2 }
            });

            double E = 200000; // Yungs modulus
            double A = 2300;   // Cross section of an edge
            // Edge == element, there is a pair of E, A for each element
            Properties = Matrix.Build.DenseOfArray(new double[,] {
                { E, A },
                { E, A },
                { E, A }
            });
            // Boundary conditions
            // Nodal freedom matrix (number of DOFs at each node)
            NF = Matrix.Build.Dense(nNodes, DOFPerNode, 1);
            NF[0, 0] = 0; NF[0, 1] = 0; // Node 1
            NF[1, 1] = 0; // Node 2

            // Count the free degrees of freedom
            CountDegreesOfFreedom();
            

            // Forces
            Load = Matrix.Build.Dense(nNodes, DOFPerNode, 0);
            // Forces in X and Y direction at node 3
            Load[2, 0] = 1200;
            Load[2, 1] = 0;
        }

        private void CountDegreesOfFreedom()
        {
            // Count the free degrees of freedom
            n = 0;
            for (int i = 0; i < nNodes; i++)
                for (int j = 0; j < DOFPerNode; j++)
                    if (Math.Abs(NF[i, j]) > double.Epsilon)
                    {
                        n++;
                        NF[i, j] = n;
                    }
        }
    }
}
