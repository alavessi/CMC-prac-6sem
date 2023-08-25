import numpy as np
from scipy.integrate import solve_ivp
from sympy import symbols, Matrix
from odeintw import odeintw


class Solver:
    def __init__(self, data: dict) -> None:
        self.__n = data['n']
        self.__x = tuple(list(symbols('x_:%d' % (self.__n + 1)))[1:])
        self.__xa = tuple(list(symbols('x_:%da' % (self.__n + 1)))[1:])
        self.__xb = tuple(list(symbols('x_:%db' % (self.__n + 1)))[1:])
        self.__a = data['a']
        self.__b = data['b']
        self.__t_star = data['t*']
        self.__p0 = data['p0']
        self.__f = data['f']
        self.__R = data['R']
        self.__inner_method = data['inner method']
        self.__external_method = data['external method']

    def __fun(self, t, x):   # returns [f_0(t, x), f_1(t, x), ...]
        fun = []
        for i in range(self.__n):
            f = self.__f[i].subs('t', t)
            for j in range(self.__n):
                f = f.subs(self.__x[j], x[j])
            fun.append(f.evalf())
        return fun

    def __find_x(self, p):  # решение векторной внутренней задачи
        if self.__t_star == self.__a:
            sol = solve_ivp(fun=self.__fun, t_span=[self.__a, self.__b], y0=p, method=self.__inner_method, t_eval=np.linspace(self.__a, self.__b, 500))
            return sol.t, sol.y
        if self.__t_star == self.__b:
            sol = solve_ivp(fun=self.__fun, t_span=[self.__b, self.__a], y0=p, method=self.__inner_method, t_eval=np.linspace(self.__b, self.__a, 500))
            return sol.t[::-1], sol.y[:, ::-1]
        sol_left = solve_ivp(fun=self.__fun, t_span=[self.__t_star, self.__a], y0=p, method=self.__inner_method, t_eval=np.linspace(self.__t_star, self.__a, 500))
        sol_right = solve_ivp(fun=self.__fun, t_span=[self.__t_star, self.__b], y0=p, method=self.__inner_method, t_eval=np.linspace(self.__t_star, self.__b, 500))
        t = np.hstack((sol_left.t[:0:-1], sol_right.t))
        x = np.hstack((sol_left.y[:, :0:-1], sol_right.y))
        return t, x

    def __rxs_matrix(self, X, t, A):  # returns A * X
        idx = max(0, min(len(A) - 1, round((t - self.__a) / (self.__b - self.__a) * len(A)) - 1))
        return np.dot(A[idx], X)

    def __find_X(self, A):  # решение матричной внутренней задачи
        sol_left = odeintw(func=self.__rxs_matrix, y0=np.eye(self.__n), t=np.linspace(self.__t_star, self.__a), args=(A,))
        sol_right = odeintw(func=self.__rxs_matrix, y0=np.eye(self.__n), t=np.linspace(self.__t_star, self.__b), args=(A,))
        return sol_left[-1], sol_right[-1]

    def __Phi(self, x):  # returns Φ(p) ≡ R(x(a,p), x(b,p))
        Phi = np.ndarray(self.__n)
        for i in range(self.__n):
            R = self.__R[i]
            for j in range(self.__n):
                R = R.subs({self.__xa[j]: x[j][0], self.__xb[j]: x[j][-1]})
            Phi[i] = R.evalf()
        return Phi

    def __dRdx(self, x, Rdx):  # подстановка x(a,p), x(b,p) в матрицы частных производных R'x, R'y
        A = Rdx
        for i in range(self.__n):
            A = A.subs({self.__xa[i]: x[i][0], self.__xb[i]: x[i][-1]})
        return np.array(A, dtype=float)

    def __solve_inner(self, J, p):  # решение внутренней задачи
        t, x = self.__find_x(p)
        A = []
        for i in range(len(t)):
            A.append(J.subs('t', t[i]))
            for j in range(self.__n):
                A[i] = A[i].subs(self.__x[j], x[j][i])
        Xa, Xb = self.__find_X(np.array(A, dtype=float))
        Rdx = Matrix(self.__R).jacobian(Matrix(self.__xa))
        Rdy = Matrix(self.__R).jacobian(Matrix(self.__xb))
        dRdx = self.__dRdx(x, Rdx)
        dRdy = self.__dRdx(x, Rdy)
        dPhidp = np.dot(dRdx, Xa) + np.dot(dRdy, Xb)
        return dPhidp

    def __rhs_ext(self, mu, p, J, Phi0, window, progressbar):  # returns −[Φ′(p)]^(−1) * Φ(p0)
        progressbar["value"] = mu * 100
        window.update()
        dPhidp = self.__solve_inner(J, p)
        return (-1) * np.dot(np.linalg.inv(dPhidp), Phi0)

    def __solve_external(self, J, Phi0, window, progressbar):  # решение внешней задачи
        sol = solve_ivp(fun=self.__rhs_ext, t_span=[0, 1], y0=self.__p0, method=self.__external_method, t_eval=[1], args=(J, Phi0, window, progressbar))
        return sol.y.flatten()

    def solve(self, window, progressbar):
        J = Matrix(self.__f).jacobian(Matrix(self.__x))
        x0 = self.__find_x(self.__p0)[1]
        Phi0 = self.__Phi(x0)
        p = self.__solve_external(J, Phi0, window, progressbar)
        print(f"ans = {p}")
        if self.__t_star == self.__a:
            sol = solve_ivp(fun=self.__fun, t_span=[self.__a, self.__b], y0=p, t_eval=np.linspace(self.__a, self.__b, 100))
            return sol.t, sol.y
        if self.__t_star == self.__b:
            sol = solve_ivp(fun=self.__fun, t_span=[self.__b, self.__a], y0=p, t_eval=np.linspace(self.__b, self.__a, 100))
            return sol.t[::-1], sol.y[:, ::-1]
        sol_left = solve_ivp(fun=self.__fun, t_span=[self.__t_star, self.__a], y0=p, t_eval=np.linspace(self.__t_star, self.__a, 100))
        sol_right = solve_ivp(fun=self.__fun, t_span=[self.__t_star, self.__b], y0=p, t_eval=np.linspace(self.__t_star, self.__b, 100))
        t = np.hstack((sol_left.t[:0:-1], sol_right.t))
        x = np.hstack((sol_left.y[:, :0:-1], sol_right.y))
        return t, x

    def calculate_functional(self, fun, t_res, x_res):
        fun_res = []
        for i in range(len(t_res)):
            fun_res.append(fun.subs('t', t_res[i]))
            for j in range(self.__n):
                fun_res[i] = fun_res[i].subs(self.__x[j], x_res[j][i])
        integral = 0
        for i in range(len(t_res) - 1):
            integral += (fun_res[i] + fun_res[i + 1]) * (t_res[i + 1] - t_res[i]) / 2
        return integral
