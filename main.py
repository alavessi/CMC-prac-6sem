#!/usr/bin/python3

import argparse

from gui import GUI


def main():
    parser = argparse.ArgumentParser(description='Решение краевой задачи методом продолжения по параметру')
    parser.add_argument('--nmax', type=int, default=6, help='Максимальныая размерность задачи')
    args = parser.parse_args()
    nmax = args.nmax
    if nmax < 1:
        return
    app = GUI(nmax)
    app.execute()


if __name__ == '__main__':
    main()
