import numpy as np
from scipy.integrate import solve_ivp
from sympy import diff, symbols, Function, Matrix


class Solver:
    def __init__(self, data: dict, t, x):
        self.n_ = data['n']  # нужно ли это поле? это ведь len(self.f_) (или len(self.R_))
        self.x_ = x
        self.t_ = t
        self.p_ = symbols('p_:%d' % self.n_)
        self.a_ = data['a']
        self.b_ = data['b']
        self.t_star_ = data['t*']
        self.p0_ = data['p0']
        self.f_ = data['f']
        self.R_ = data['R']
        # print(type(self.n_))
        # print(type(self.x_[0]))
        # print(type(self.t_))
        # print(type(self.p_[0]))
        # print(type(self.a_))
        # print(type(self.b_))
        # print(type(self.t_star_))
        # print(type(self.p0_[0]))
        # print(type(self.f_[0]))
        # print(type(self.R_[0]))

    def __f(self, t, x):  # f(x,t) returns [f_0(x,t), f_1(x,t), ...]
        fun = []
        for i in range(self.n_):
            f = self.f_[i].subs('t', t)
            for j in range(self.n_):
                f = f.subs(self.x_[j], x[j])
            fun.append(f.evalf())
        return fun

    def __find_x(self, tau):
        t_grid = np.linspace(self.a_, self.b_, int((self.b_ - self.a_) / tau))
        sol = solve_ivp(self.__f, [self.a_, self.b_], self.p0_, t_eval=t_grid, dense_output=True)
        return (sol.t, sol.y)

    def __find_X(self):
        return

    def __solve_inner(self, J):
        tau = 0.1
        t, x = self.__find_x(tau)
        A = []
        for i in range(len(t)):
            A.append(J.subs('t', t[i]))
            for j in range(self.n_):
                A[i] = A[i].subs(self.x_[j], x[j][i])
        # print(A)
        return (t, x)
        self.__find_X()

    def solve(self):
        J = Matrix(self.f_).jacobian(Matrix(self.x_))
        # print(J)
        return self.__solve_inner(J)