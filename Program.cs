namespace FEM2D
{
    internal class Program
    {

        static void Main(string[] args)
        {
            RunFemBasic(Model.Scenario.Simple);
            RunFemBasic(Model.Scenario.Bridge);
        }

        private static void RunFemBasic(Model.Scenario scenario)
        {
            Model m = new(scenario);

            FEM fem = new(m);

            // Global stiffness matrix
            fem.AssembleGlobalStiffness();

            // Form global force vector
            // Global forces
            fem.FormTrussF();

            // solve for unknown displacements
            var delta = fem.Solve();

            // Extract nodal displacements
            var nodeDisp = fem.ExtractNodalDisplacements(delta);
            Console.WriteLine("Nodal Displacements:" + nodeDisp);

            var force = fem.BuildLocalForces(delta);
            Console.WriteLine("Element forces in local coordinate system \n (positive - Tension; negative - Compression): " + force);
        }
    }
}