#!/usr/bin/python3

import argparse
import sympy

from data_parser import Parser, is_number
from solver import Solver
from plotter import Plotter


def main():
    sympy.init_printing(use_unicode=True)

    parser = argparse.ArgumentParser(description='Решение кревой задачи методом продолжения по параметру')
    parser.add_argument('--problem_name', type=str, default='', help='Вы можете выбрать одну из известных задач (по умолчанию данные вводятся в полном объёме с клавиатуры)')
    args = parser.parse_args()

    parser = Parser(args.problem_name)
    data = parser.parse_input_data()

    solver = Solver(data)
    t, x = solver.solve()
    print("t =", t)
    print("x =", x)
    
    plotter = Plotter(t, x)
    plotter.plot_x()


if __name__ == '__main__':
    main()
