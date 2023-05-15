#!/usr/bin/python3

import sympy

from data_parser import Parser
from solver import Solver
from plotter import Plotter


def main():
    sympy.init_printing(use_unicode=True)

    parser = Parser()
    data = parser.parse_input_data()

    solver = Solver(data)
    t, x = solver.solve()

    plotter = Plotter(t, x)
    plotter.plot_x()


if __name__ == '__main__':
    main()
