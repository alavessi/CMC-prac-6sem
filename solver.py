import numpy as np
from scipy.integrate import solve_ivp
from sympy import symbols, Matrix
from odeintw import odeintw


class Solver:
    def __init__(self, data: dict):
        self.n_ = data['n']
        self.x_ = data['x']
        self.xa_ = data['x(a)']
        self.xb_ = data['x(b)']
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

    def __find_x(self, p):
        sol_left = solve_ivp(fun=self.__f, t_span=[self.t_star_, self.a_], y0=p)
        sol_right = solve_ivp(fun=self.__f, t_span=[self.t_star_, self.b_], y0=p)
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

    def __Phi(self, x):
        Phi = np.ndarray(self.n_)
        for i in range(self.n_):
            R = self.R_[i]
            for j in range(self.n_):
                R = R.subs({self.xa_[j]: x[j][0], self.xb_[j]: x[j][-1]})
            Phi[i] = R.evalf()
        return Phi

    def __dRdx(self, x, Rdx):  # подстановка x(a,p), x(b,p) в матрицы частных производных R'x, R'y
        A = Rdx
        for i in range(self.n_):
            A = A.subs({self.xa_[i]: x[i][0], self.xb_[i]: x[i][-1]})
        return A

    def __rhs_ext(self, mu, p, J, Phi0):
        dPhidp = self.__solve_inner(J, p)
        # print(dPhidp)
        return (-1) * np.linalg.inv(dPhidp)@Phi0

    def __solve_external(self, J, Phi0):
        sol = solve_ivp(fun=self.__rhs_ext, t_span=[0, 1], y0=self.p0_, args=(J, Phi0))
        return sol.y

    def __solve_inner(self, J, p):
        t, x = self.__find_x(p)
        A = []
        for i in range(len(t)):
            A.append(J.subs('t', t[i]))
            for j in range(self.n_):
                A[i] = A[i].subs(self.x_[j], x[j][i])
        A = np.array(A, dtype='float64')
        X = self.__find_X(A)
        Rdx = Matrix(self.R_).jacobian(Matrix(self.xa_))
        Rdy = Matrix(self.R_).jacobian(Matrix(self.xb_))
        dRdx = self.__dRdx(x, Rdx)
        dRdy = self.__dRdx(x, Rdy)
        dPhidp = dRdx@X[0] + dRdy@X[-1]
        return np.array(dPhidp, dtype=float)

    def solve(self):
        J = Matrix(self.f_).jacobian(Matrix(self.x_))
        x_0 = self.__find_x(self.p0_)[1]
        print(x_0)
        Phi0 = self.__Phi(x_0)
        p = self.__solve_external(J, Phi0)
        print(p)
        sol_left = solve_ivp(fun=self.__f, t_span=[self.t_star_, self.a_], y0=p[:, -1])
        sol_right = solve_ivp(fun=self.__f, t_span=[self.t_star_, self.b_], y0=p[:, -1])
        t = np.hstack((sol_left.t[:0:-1], sol_right.t))
        ans = np.hstack((sol_left.y[:, :0:-1], sol_right.y))
        return t, ans
