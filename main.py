#!/usr/bin/python3

import argparse

from data_parser import Parser
from solver import Solver
from plotter import Plotter


def main():
    help_msg = 'Вы можете выбрать одну из известных задач (по умолчанию данные вводятся в полном объёме с клавиатуры)'
    parser = argparse.ArgumentParser(description='Решение краевой задачи методом продолжения по параметру')
    parser.add_argument('--problem_name', type=str, default='', help=help_msg)
    args = parser.parse_args()

    parser = Parser(args.problem_name)
    data = parser.parse_input_data()

    solver = Solver(data)
    t, x = solver.solve()
    # print("t =", t)
    # print("x =", x)

    plotter = Plotter(t, x)
    plotter.plot_x()


if __name__ == '__main__':
    main()
