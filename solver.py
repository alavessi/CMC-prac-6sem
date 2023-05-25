import numpy as np
from scipy.integrate import solve_ivp
from sympy import symbols, Matrix
from odeintw import odeintw


class Solver:
    def __init__(self, data: dict):
        self.n_ = data['n']  # нужно ли это поле? это ведь len(self.f_) (или len(self.R_))
        self.x_ = data['x']
        self.xa_ = data['x(a)']
        self.xb_ = data['x(b)']
        self.p_ = symbols('p_:%d' % self.n_)
        self.a_ = data['a']
        self.b_ = data['b']
        self.t_star_ = data['t*']
        self.p0_ = data['p0']
        self.f_ = data['f']
        self.R_ = data['R']

    def __f(self, t, x):   # f(t, x) returns [f_0(t, x), f_1(t, x), ...]
        fun = []
        for i in range(self.n_):
            f = self.f_[i].subs('t', t)
            for j in range(self.n_):
                f = f.subs(self.x_[j], x[j])
            fun.append(f.evalf())
        return fun

    def __find_x(self):
        sol_left = solve_ivp(fun=self.__f, t_span=[self.t_star_, self.a_], y0=self.p0_, method='Radau', dense_output=True)
        sol_right = solve_ivp(fun=self.__f, t_span=[self.t_star_, self.b_], y0=self.p0_, method='Radau', dense_output=True)
        t = np.hstack((sol_left.t[:0:-1], sol_right.t))
        x = np.hstack((sol_left.y[:, :0:-1], sol_right.y))
        return t, x

    def __matprod(self, X, t, A):  # returns A * X
        idx = max(0, min(len(A) - 1, round(((t - self.a_) / (self.b_ - self.a_) * (len(A) - 1)))))
        return A[idx]@X

    def __find_X(self, A):
        sol_left = odeintw(func=self.__matprod, y0=np.eye(self.n_), t=np.linspace(self.t_star_, self.a_), args=(A,))
        sol_right = odeintw(func=self.__matprod, y0=np.eye(self.n_), t=np.linspace(self.t_star_, self.b_), args=(A,))
        return np.vstack((sol_left[:0:-1], sol_right))

    def __solve_inner(self, J):
        t, x = self.__find_x()
        A = []
        for i in range(len(t)):
            A.append(J.subs('t', t[i]))
            for j in range(self.n_):
                A[i] = A[i].subs(self.x_[j], x[j][i])
        A = np.array(A, dtype='float64')
        X = self.__find_X(A)
        return (t, x)

    def solve(self):
        J = Matrix(self.f_).jacobian(Matrix(self.x_))
        return self.__solve_inner(J)
