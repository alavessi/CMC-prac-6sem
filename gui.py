from tkinter import Tk, Menu, Frame, Button, Label, Text, LabelFrame
from tkinter.ttk import Combobox, Progressbar
from tkinter.filedialog import asksaveasfile, askopenfile
from tkinter.constants import LEFT, RIGHT, TOP, DISABLED, NORMAL, END, WORD, BOTH, CENTER

from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application, function_exponentiation, convert_xor
from PIL import Image, ImageTk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from pandas import DataFrame

from solver import Solver


METHODS = ["RK45", "RK23", "DOP853", "Radau", "BDF", "LSODA"]
TRANSFORMATIONS = standard_transformations + (implicit_multiplication_application, convert_xor, function_exponentiation)


def get_frame(frame: Frame, text: str, side: str = LEFT, padx: int = None):
    frame = Frame(frame)
    frame.pack(side=TOP)
    Label(frame, text=text, font=("Inter", 15)).pack(side=side, expand=True, pady=5, padx=padx)
    return frame


def frame_time_bounds(frame: Frame):
    frame = get_frame(frame=frame, text="Начальный и конечный моменты времени:")
    init_a = Text(frame, width=5, height=1, font=("Arial", 15))
    init_a.pack(side=LEFT, expand=True, pady=5, padx=10)
    init_b = Text(frame, width=5, height=1, font=("Arial", 15))
    init_b.pack(side=LEFT, expand=True, pady=5)
    return init_a, init_b


def frame_time_star(frame: Frame) -> Text:
    frame = get_frame(frame=frame, text="Момент времени t* = ")
    init_time_star = Text(frame, width=20, height=1, font=("Arial", 15))
    init_time_star.pack(side=LEFT, expand=True, pady=5)
    return init_time_star


def frame_inner_method(frame: Frame) -> Combobox:
    frame = get_frame(frame=frame, text="Метод решения внутренней задачи:", padx=15)
    inner_method = Combobox(frame, width=5, state="readonly", values=METHODS)
    inner_method.pack(pady=5, side=LEFT)
    return inner_method


def frame_external_method(frame: Frame) -> Combobox:
    frame = get_frame(frame=frame, text="Метод решения внешней задачи:", padx=15)
    external_method = Combobox(frame, width=5, state="readonly", values=METHODS)
    external_method.pack(pady=5, side=LEFT)
    return external_method


def frame_init_parameter(frame: Frame) -> Text:
    frame = get_frame(frame=frame, text="Вектор начального приближения параметра:")
    init_parameter = Text(frame, width=20, height=1, font=("Arial", 15))
    init_parameter.pack(side=LEFT, expand=True, pady=5, padx=10)
    return init_parameter


def about_program():
    about_program_window = Tk()
    about_program_window.title('О программе')
    about_program_window.geometry('800x300')
    documentation = Text(about_program_window, width=150, height=100, wrap=WORD)
    documentation.pack(side=TOP)
    with open("README.md", encoding='utf-8') as f:
        documentation.insert(END, f.read())
    about_program_window.mainloop()


def about_author():
    about_author_window = Tk()
    about_author_window.title('Об авторе')
    about_author_window.geometry('700x500')
    image = Image.open("my_photo.jpg").resize((300, 400))
    photo = ImageTk.PhotoImage(image, master=about_author_window)
    label = Label(about_author_window, image=photo, anchor=CENTER)
    label.pack(fill=BOTH, expand=True)
    author_info = Text(about_author_window, wrap=WORD)
    author_info.pack(side=TOP)
    about_me = "Автор программы: Эмиров Самир, студент 313 группы ВМК МГУ, кафедры оптимального управления.\n" + \
               "email: samir.emirov.2001@mail.ru\nВесна 2023"
    author_info.insert(END, about_me)
    about_author_window.mainloop()


def help_program():
    help_window = Tk()
    help_window.title('Помощь')
    help_window.geometry('800x350')
    help_info = Text(help_window, width=600, height=350, wrap=WORD)
    help_info.pack(side=TOP)
    with open("help.txt", encoding="utf-8") as f:
        help_info.insert(END, f.read())
    help_window.mainloop()


class GUI:
    def __init__(self, maximal_dimension: int):
        self.window = Tk()
        self.window.title("Метод продолжения по параметру")
        self.window.geometry("1200x600")
        mainmenu = Menu(self.window)
        self.window.config(menu=mainmenu)

        file_menu = Menu(mainmenu, tearoff=0)
        file_menu.add_command(label="Открыть задачу...", command=self.open_from_file)
        file_menu.add_command(label="Сохранить задачу...", command=self.save_problem)
        file_menu.add_command(label="Сохранить решение...", command=self.save_solution)
        file_menu.add_command(label="Выход", command=self.window.destroy)
        mainmenu.add_cascade(label="Файл", menu=file_menu)

        help_menu = Menu(mainmenu, tearoff=0)
        help_menu.add_command(label="О программе", command=about_program)
        help_menu.add_command(label="Об авторе", command=about_author)
        help_menu.add_command(label="Помощь", command=help_program)
        mainmenu.add_cascade(label="Информация", menu=help_menu)

        self.init_odes = []
        self.init_boundary_conditions = []
        self.maximal_dimension = maximal_dimension

    def __get_data_from_gui(self):
        f = []
        for i in range(self.maximal_dimension):
            expr = self.init_odes[i].get("1.0", "end-1c")
            if expr:
                f.append(parse_expr(expr, transformations=TRANSFORMATIONS, evaluate=True))
        R = [None] * len(f)
        for i in range(len(f)):
            expr = self.init_boundary_conditions[i].get("1.0", "end-1c")
            R[i] = parse_expr(expr, transformations=TRANSFORMATIONS, evaluate=True)
        a = float(self.init_a.get("1.0", "end-1c"))
        b = float(self.init_b.get("1.0", "end-1c"))
        t_star = float(self.init_time_star.get("1.0", "end-1c"))
        p0 = list(map(float, self.init_parameter.get("1.0", "end-1c").split(",")))
        inner_method = self.inner_method.get()
        external_method = self.external_method.get()
        return {'n': len(f), 'a': a, 'b': b, 't*': t_star, 'p0': p0, 'f': f, 'R': R,
                'inner method': inner_method, 'external method': external_method}

    def __exec(self):
        self.window.update()
        self.button_draw["state"] = DISABLED
        self.button_functional["state"] = DISABLED
        self.window.update()
        data = self.__get_data_from_gui()
        self.solver = Solver(data)
        self.t_res, self.x_res = self.solver.solve(self.window, self.progressbar)
        self.button_draw["state"] = NORMAL
        self.button_functional["state"] = NORMAL
        self.window.update()

    def open_from_file(self):
        for value in *self.init_odes, *self.init_boundary_conditions, self.init_a, self.init_b, self.init_time_star, self.init_parameter:
            value.delete("1.0", "end-1c")
        filetypes = [('text files', '*.txt'), ('All files', '*.*')]
        f = askopenfile(title='Open a file', initialdir='./problems', filetypes=filetypes)
        s = f.readline()
        i = 0
        while s != '\n':
            self.init_odes[i].insert("1.0", s.strip())
            i += 1
            s = f.readline()
        i = 0
        s = f.readline()
        while s != '\n':
            self.init_boundary_conditions[i].insert("1.0", s.strip())
            i += 1
            s = f.readline()
        for value in [self.init_a, self.init_b, self.init_time_star, self.init_parameter]:
            value.insert("1.0", f.readline().strip())

    def save_problem(self):
        filetypes = [("All Files", "*.*"), ("Text Documents", "*.txt")]
        f = asksaveasfile(initialfile='untitled.txt', defaultextension=".txt", filetypes=filetypes)
        for i in range(self.maximal_dimension):
            ode_rhs = self.init_odes[i].get("1.0", "end-1c")
            if ode_rhs == "":
                break
            f.write(ode_rhs)
            f.write("\n")
        f.write("\n")
        for i in range(self.maximal_dimension):
            boundary_lhs = self.init_boundary_conditions[i].get("1.0", "end-1c")
            if boundary_lhs == "":
                break
            f.write(boundary_lhs)
            f.write("\n")
        f.write("\n")
        f.write(self.init_a.get("1.0", "end-1c"))
        f.write("\n")
        f.write(self.init_b.get("1.0", "end-1c"))
        f.write("\n")
        f.write(self.init_time_star.get("1.0", "end-1c"))
        f.write("\n")
        f.write(self.init_parameter.get("1.0", "end-1c"))


    def save_solution(self):
        solution_dict = {"t": self.t_res}
        i = 1
        for x in self.x_res:
            solution_dict[f"x_{i}"] = x
            i += 1
        df = DataFrame(solution_dict)
        filetypes = [("All Files", "*.*"), ("Text Documents", "*.txt"), ("Table documents", "*.csv")]
        f = asksaveasfile(initialfile="untitled.csv", defaultextension=".csv", filetypes=filetypes)
        df.to_csv(path_or_buf=f)

    def calc_functional(self):
        expr = self.integral_txt.get("1.0", "end-1c")
        if not expr:
            return
        fun = parse_expr(expr, transformations=TRANSFORMATIONS, evaluate=True)
        functional = self.solver.calculate_functional(fun, self.t_res, self.x_res)
        self.integral_value.config(text=str(functional))

    def __frame_odes(self, frame: Frame):
        get_frame(frame=frame, text="Введите систему дифференциальных уравнений:", side=TOP)
        for i in range(self.maximal_dimension):
            fr = Frame(frame)
            fr.pack(side=TOP)
            self.init_odes.append(Text(fr, width=30, height=1, font=("Arial", 15)))
            Label(fr, text=f"dx_{i + 1}/dt =", font=("Inter", 15)).pack(side=LEFT, expand=True, pady=5, padx=40)
            self.init_odes[i].pack(side=LEFT, fill="x", expand=True)

    def __frame_boundary_conditions(self, frame: Frame):
        frame = get_frame(frame=frame, text="Введите систему краевых условий:", side=TOP)
        for i in range(self.maximal_dimension):
            fr = Frame(frame)
            fr.pack(side=TOP)
            self.init_boundary_conditions.append(Text(fr, width=30, height=1, font=("Arial", 15)))
            self.init_boundary_conditions[i].pack(side=LEFT, fill="x", expand=True)
            Label(fr, text="= 0", font=("Inter", 15)).pack(side=LEFT, expand=True, pady=5, padx=40)

    def __draw(self, dr1, dr2):
        conv = {f"x_{i + 1}": self.x_res[i] for i in range(len(self.x_res))}
        conv["t"] = self.t_res
        plot = Tk()
        plot.title("График решения задачи")
        fig = Figure(figsize=(5, 5), dpi=100)
        plot1 = fig.add_subplot(111)
        plot1.plot(conv[dr1.get()], conv[dr2.get()])
        plot1.grid()
        canvas = FigureCanvasTkAgg(fig, master=plot)
        canvas.draw()
        canvas.get_tk_widget().pack()
        toolbar = NavigationToolbar2Tk(canvas, plot)
        toolbar.update()
        canvas.get_tk_widget().pack()
        plot.mainloop()

    def __frame_plot(self, frame: Frame):
        frame = get_frame(frame=frame, text="Оси графиков:", padx=15)
        axes = ["t"] + [f"x_{i + 1}" for i in range(self.maximal_dimension)]
        dr1 = Combobox(frame, width=5, state="readonly", values=axes)
        dr2 = Combobox(frame, width=5, state="readonly", values=axes)
        dr1.pack(pady=5, side=LEFT)
        dr2.pack(pady=5, side=LEFT)
        self.button_draw = Button(frame, text="Нарисовать", font=("Arial", 15), command=lambda: self.__draw(dr1, dr2), state=DISABLED)
        self.button_draw.pack(side=LEFT, fill='both', pady=5, padx=15)

    def create_interface(self):
        frame_left = Frame(self.window)
        frame_left.pack(side=LEFT)
        self.__frame_odes(frame_left)
        self.init_a, self.init_b = frame_time_bounds(frame_left)
        self.init_time_star = frame_time_star(frame_left)
        self.inner_method = frame_inner_method(frame_left)
        self.external_method = frame_external_method(frame_left)
        frame_right = Frame(self.window)
        frame_right.pack(side=RIGHT)
        self.__frame_boundary_conditions(frame_right)
        self.init_parameter = frame_init_parameter(frame_right)
        return frame_right

    def execute(self):
        frame = self.create_interface()
        self.progressbar = Progressbar(frame, orient="horizontal", length=200, value=0)
        Button(frame, text="Решить задачу", font=("Arial", 15), command=self.__exec).pack(side=TOP, fill="y", pady=5)
        Label(frame, text="Алгоритм выполняется...", font=("Inter", 15)).pack(side=TOP, expand=True, pady=5, padx=40)
        self.progressbar.pack(side=TOP)

        self.__frame_plot(frame)

        label_functional = LabelFrame(frame, text="Интегральный функционал от решения задачи:", font=("Inter", 15))
        label_functional.pack(side="left", expand=True, pady=5)
        Label(label_functional, text="J(x) = ", font=("Inter", 15)).pack(side="left", expand=True, pady=5)
        self.integral_txt = Text(label_functional, width=10, height=1, font=("Times", 15))
        self.integral_txt.pack(side="left", expand=True, pady=5, padx=10)
        self.integral_value = Label(label_functional, text="", font=("Times", 15), state=DISABLED)
        self.button_functional = Button(label_functional, text="Вычислить", font=("Arial", 15), width=10, command=self.calc_functional, state=DISABLED)
        self.button_functional.pack(side='left', fill=BOTH, padx=5)
        self.integral_value.pack(side='left', expand=True, pady=5, padx=10)

        self.window.mainloop()
