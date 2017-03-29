class Solvers_Class():

    def __init__(self):
        pass

    def ProperDistanceTabulate(self, Input_Param, z_max):
        from Distance_Solver import Proper_Distance_Tabulate
        Proper_Distance_Tabulate(Input_Param, z_max)

    def LxTxSolver(self, Halos):
        from LxTx_Solver import LxTx_Solver
        solver = LxTx_Solver()
        solver.solve(Halos)


