from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application, function_exponentiation, convert_xor


def is_number(str) -> bool:
    try:
        float(str)
        return True
    except ValueError:
        return False


class Parser:
    def __init__(self, problem_name=''):
        self.problem_name_ = problem_name
        self.__transformations = standard_transformations + (implicit_multiplication_application, convert_xor, function_exponentiation)

    def parse_input_data(self) -> dict:
        self.n_ = self.__parse_dimension()
        self.a_, self.b_ = self.__parse_time_interval()
        self.t_star_ = self.__parse_time_star()
        self.p0_ = self.__parse_parameter()
        self.f_ = self.__parse_function_f()
        self.R_ = self.__parse_function_R()
        return {'n': self.n_,
                'a': self.a_,
                'b': self.b_,
                't*': self.t_star_,
                'p0': self.p0_,
                'f': self.f_,
                'R': self.R_}

    def __parse_dimension(self):
        match self.problem_name_:
            case "Краевая задача двух тел":
                return 4
            case "Предельные циклы в системе Эквейлера":
                return 4
            case "Функционал типа 'энергия' для трехкратного интегратора":
                return 6
            case _:
                s = input("Задайте размерность задачи: n = ")
                while not (s.isdigit() and int(s) != 0):
                    print("ЭТО ДОЛЖНО БЫТЬ НАТУРАЛЬНОЕ ЧИСЛО")
                    s = input("n = ")
                return int(s)

    def __parse_time_interval(self):
        match self.problem_name_:
            case "Предельные циклы в системе Эквейлера":
                return 0, 1
            case "Краевая задача двух тел" | "Функционал типа 'энергия' для трехкратного интегратора":
                str_T = input("Задайте конечный момент времени T = ")
                while not (is_number(str_T) and float(str_T) > 0):
                    print("ЭТО ДОЛЖНО БЫТЬ ПОЛОЖИТЕЛЬНОЕ ЧИСЛО")
                    str_T = input("T = ")
                return 0, float(str_T)
            case _:
                print("Задайте начальный момент времени a и конечный момент времени b")
                str_a = input("a = ")
                str_b = input("b = ")
                while not (is_number(str_a) and is_number(str_b) and float(str_a) >= 0 and float(str_a) < float(str_b)):
                    print("a И b ДОЛЖНЫ БЫТЬ ЧИСЛАМИ, ПРИЧЕМ 0 <= a < b")
                    str_a = input("a = ")
                    str_b = input("b = ")
                return float(str_a), float(str_b)

    def __parse_time_star(self):
        str_t = (input(f"Задайте момент времени t* из отрезка [{self.a_}, {self.b_}]: t* = "))
        while not (is_number(str_t) and self.a_ <= float(str_t) and float(str_t) <= self.b_):
            str_t = (input(f"ЭТО ДОЛЖНО БЫТЬ ЧИСЛО ИЗ ОТРЕЗКА [{self.a_}, {self.b_}]: t* = "))
        return float(str_t)

    def __parse_parameter(self):
        p0 = []
        for i in range(self.n_):
            str_p = input(f"Задайте начальное приближение параметра: p0_{i + 1} = ")
            while not is_number(str_p):
                str_p = (input(f"ЭТО ДОЛЖНО БЫТЬ ЧИСЛО: p0_{i + 1} = "))
            p0.append(float(str_p))
        return p0

    def __parse_function_f(self):
        match self.problem_name_:
            case "Краевая задача двух тел":
                funcs = ['x_3', 'x_4', '-x_1/(x_1**2 + x_2**2)**(3/2)', '-x_2/(x_1**2 + x_2**2)**(3/2)']
            case "Предельные циклы в системе Эквейлера":
                funcs = ['x_3*x_2', 'x_3*(-x_1+sin(x_2))', '0', '0']
            case "Функционал типа 'энергия' для трехкратного интегратора":
                v = float(input("Задайте малый параметр v = "))
                funcs = ['x_2', 'x_3', f'1/2*(({v}+(x_6+1)**2)+({v}+(x_6-1)**2))', '0', '-x_4', '-x_5']
            case _:
                print("Введите функции f(x, t) из правых частей ОДУ")
                funcs = [input(f"dx_{i + 1}/dt = ") for i in range(self.n_)]
        return [parse_expr(f, transformations=self.__transformations) for f in funcs]

    def __parse_function_R(self):
        match self.problem_name_:
            case "Краевая задача двух тел":
                a1 = float(input("a1 = "))
                a2 = float(input("a2 = "))
                b1 = float(input("b1 = "))
                b2 = float(input("b2 = "))
                R = [f'x_1a-{a1}', f'x_2a-{a2}', f'x_1b-{b1}', f'x_2b-{b2}']
            case "Предельные циклы в системе Эквейлера":
                R = ['x_1a-x_4a', 'x_1b-x_4b', 'x_2a', 'x_2b']
            case "Функционал типа 'энергия' для трехкратного интегратора":
                R = ['x_1a-1', 'x_2a', 'x_3a', 'x_1b', 'x_2b', 'x_3b']
            case _:
                print("Введите функции R из краевых условий R(x(a), x(b)) = 0")
                R = [input(f"R_{i + 1}(x(a), x(b)) = ") for i in range(self.n_)]
        return [parse_expr(r, transformations=self.__transformations) for r in R]
