import sympy

from data_parser import Parser
from solver import Solver
from plotter import Plotter

if __name__ == '__main__':
    sympy.init_printing(use_unicode=True)
    t = sympy.symbols('t')

    parser = Parser(t)
    data = parser.parse_input_data()

    solver = Solver(data, t, data['x'])
    t, x = solver.solve()

    plotter = Plotter(t, x)
    plotter.plot_x()
