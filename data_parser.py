from sympy import symbols
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application, function_exponentiation, convert_xor


def is_number(str) -> bool:
    try:
        float(str)
        return True
    except ValueError:
        return False


class Parser:
    def __init__(self):
        self.n_ = 0
        self.a_, self.b_ = 0, 0
        self.t_star_ = 0
        self.x_ = symbols('x:%d' % self.n_)
        self.p0_ = []
        self.f_ = []
        self.R_ = []
        self.__transformations = standard_transformations + (implicit_multiplication_application, convert_xor, function_exponentiation)

    def parse_input_data(self) -> dict:
        self.__parse_dimension()
        self.__parse_time_interval()
        self.__parse_time_star()
        self.__parse_parameter()
        self.__parse_function_f()
        self.__parse_function_R()
        return {'n': self.n_, 'a': self.a_, 'b': self.b_, 't*': self.t_star_, 'x': self.x_, 'p0': self.p0_, 'f': self.f_, 'R': self.R_}

    def __parse_dimension(self):
        s = input("Задайте размерность задачи: n = ")
        while not (s.isdigit() and int(s) != 0):
            print("ЭТО ДОЛЖНО БЫТЬ НАТУРАЛЬНОЕ ЧИСЛО")
            s = input("n = ")
        self.n_ = int(s)
        self.x_ = tuple(list(symbols('x_:%d' % (self.n_ + 1)))[1:])

    def __parse_time_interval(self):
        print("Задайте начальный момент времени a и конечный момент времени b")
        str_a = input("a = ")
        str_b = input("b = ")
        while not (is_number(str_a) and is_number(str_b) and float(str_a) >= 0 and float(str_a) < float(str_b)):
            print("a И b ДОЛЖНЫ БЫТЬ ЧИСЛАМИ, ПРИЧЕМ 0 <= a < b")
            str_a = input("a = ")
            str_b = input("b = ")
        self.a_, self.b_ = float(str_a), float(str_b)

    def __parse_time_star(self):
        str_t = (input(f"Задайте момент времени t* из отрезка [{self.a_}, {self.b_}]: t* = "))
        while not (is_number(str_t) and self.a_ <= float(str_t) and float(str_t) <= self.b_):
            str_t = (input(f"ЭТО ДОЛЖНО БЫТЬ ЧИСЛО ИЗ ОТРЕЗКА [{self.a_}, {self.b_}]: t* = "))
        self.t_star_ = float(str_t)

    def __parse_parameter(self):
        for i in range(self.n_):
            str_p = input(f"Задайте начальное приближение параметра: p0_{i + 1} = ")
            while not is_number(str_p):
                str_p = (input(f"ЭТО ДОЛЖНО БЫТЬ ЧИСЛО: p0_{i + 1} = "))
            self.p0_.append(float(str_p))

    def __parse_function_f(self):
        print("Введите функции f(x, t) из правых частей ОДУ")
        for i in range(self.n_):
            self.f_.append(parse_expr(input(f"dx_{i + 1}/dt = "), transformations=self.__transformations))

    def __parse_function_R(self):
        print("Введите функции R из краевых условий R(x(a), x(b)) = 0. Введем обозначения: y = x(a), z = x(b)")
        y, z = symbols('y z')
        for i in range(self.n_):
            self.R_.append(parse_expr(input(f"R_{i + 1}(y, z) = "), transformations=self.__transformations))
