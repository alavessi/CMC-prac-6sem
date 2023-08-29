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


def frame_time_interval(frame: Frame):
    frame = get_frame(frame=frame, text="Временной промежуток:")
    init_a = Text(frame, width=5, height=1, font=("Arial", 15))
    init_a.pack(side=LEFT, expand=True, pady=5, padx=10)
    init_b = Text(frame, width=5, height=1, font=("Arial", 15))
    init_b.pack(side=LEFT, expand=True, pady=5)
    return init_a, init_b


def frame_time_star(frame: Frame) -> Text:
    frame = get_frame(frame=frame, text="Момент времени t* = ")
    init_time = Text(frame, width=20, height=1, font=("Arial", 15))
    init_time.pack(side=LEFT, expand=True, pady=5)
    return init_time


def frame_inner_method(frame: Frame) -> Combobox:
    frame = get_frame(frame=frame, text="Метод для решения внутренней задачи:", padx=15)
    inner_method = Combobox(frame, width=5, state="readonly", values=METHODS)
    inner_method.pack(pady=5, side=LEFT)
    return inner_method


def frame_external_method(frame: Frame) -> Combobox:
    frame = get_frame(frame=frame, text="Метод для решения внешней задачи:", padx=15)
    external_method = Combobox(frame, width=5, state="readonly", values=METHODS)
    external_method.pack(pady=5, side=LEFT)
    return external_method


def frame_initial_value(frame: Frame) -> Text:
    frame = get_frame(frame=frame, text="Вектор начального приближения параметра:")
    init_value = Text(frame, width=20, height=1, font=("Arial", 15))
    init_value.pack(side=LEFT, expand=True, pady=5, padx=10)
    return init_value


def about_program():
    about_program_window = Tk()
    about_program_window.title('О программе')
    about_program_window.geometry('800x300')
    text = Text(about_program_window, width=150, height=100, wrap=WORD)
    text.pack(side=TOP)
    with open("README.md", encoding='utf-8') as f:
        documentation = f.read()
    text.insert(END, documentation)
    about_program_window.mainloop()


def about_author():
    about_author_window = Tk()
    about_author_window.title('Об авторе')
    about_author_window.geometry('700x500')
    im = Image.open("my_photo.jpg")
    im = im.resize((300, 400))
    ph = ImageTk.PhotoImage(im, master=about_author_window)
    label = Label(about_author_window, image=ph, anchor=CENTER)
    label.pack(fill=BOTH, expand=True)
    text = Text(about_author_window, wrap=WORD)
    text.pack(side=TOP)
    about_me = "Автор программы: Эмиров Самир, студент 313 группы ВМК МГУ, кафедры оптимального управления.\n" + \
               "email: samir.emirov.2001@mail.ru\nВесна 2023"
    text.insert(END, about_me)
    about_author_window.mainloop()


def help_program():
    help_window = Tk()
    help_window.title('Помощь')
    help_window.geometry('800x350')
    text = Text(help_window, width=600, height=350, wrap=WORD)
    text.pack(side=TOP)
    with open("help.txt", encoding="utf-8") as f:
        user_guide = f.read()
    text.insert(END, user_guide)
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
        
        self.texts_diff = list()
        self.texts_edge = list()
        self.maximal_dimension = maximal_dimension

    def __get_data_from_gui(self):
        f = list()
        for i in range(self.maximal_dimension):
            expr = self.texts_diff[i].get("1.0", "end-1c")
            if expr:
                f.append(parse_expr(expr, transformations=TRANSFORMATIONS, evaluate=True))
        R = list()
        for i in range(self.maximal_dimension):
            expr = self.texts_edge[i].get("1.0", "end-1c")
            if expr:
                R.append(parse_expr(expr, transformations=TRANSFORMATIONS, evaluate=True))
        self.n_ = len(f)
        a = eval(self.init_a.get("1.0", "end-1c"))
        b = eval(self.init_b.get("1.0", "end-1c"))
        t_star = eval(self.init_time.get("1.0", "end-1c"))
        p0 = eval(self.init_value.get("1.0", "end-1c"))
        if self.n_ == 1:
            p0 = [p0]
        else:
            p0 = list(p0)
        inner_method = self.inner_method.get()
        external_method = self.external_method.get()
        return {'n': self.n_, 'a': a, 'b': b, 't*': t_star, 'p0': p0, 'f': f, 'R': R,
                'inner method': inner_method, 'external method': external_method}

    def exec(self):
        self.window.update()
        self.button_draw["state"] = DISABLED
        self.button_functional["state"] = DISABLED
        self.window.update()
        data = self.__get_data_from_gui()
        self.solver = Solver(data)
        self.t_res, self.x_res = self.solver.solve(self.window, self.progressbar)
        print(f"t_res = {self.t_res}")
        print(f"x_res = {self.x_res}")
        self.button_draw["state"] = NORMAL
        self.button_functional["state"] = NORMAL
        self.window.update()

    def clear_fields(self):
        for text_value in *self.texts_diff, *self.texts_edge:
            text_value.delete("1.0", "end-1c")
        self.init_a.delete("1.0", "end-1c")
        self.init_b.delete("1.0", "end-1c")
        self.init_time.delete("1.0", "end-1c")
        self.init_value.delete("1.0", "end-1c")

    def open_from_file(self):
        self.clear_fields()
        filetypes = (('text files', '*.txt'), ('All files', '*.*'))
        f = askopenfile(
            title='Open a file',
            initialdir='./problems',
            filetypes=filetypes)
        s = f.readline()
        i = 0
        while s != '\n':
            self.texts_diff[i].insert("1.0", s.strip())
            i += 1
            s = f.readline()
        i = 0
        s = f.readline()
        while s != '\n':
            self.texts_edge[i].insert("1.0", s.strip())
            i += 1
            s = f.readline()
        self.init_a.insert("1.0", f.readline().strip())
        self.init_b.insert("1.0", f.readline().strip())
        self.init_value.insert("1.0", f.readline().strip())
        self.init_time.insert("1.0", f.readline().strip())

    def save_problem(self):
        f = asksaveasfile(
            initialfile='untitled.txt',
            defaultextension=".txt",
            filetypes=[("All Files", "*.*"), ("Text Documents", "*.txt")]
        )
        for i in range(self.maximal_dimension):
            ode_rhs = self.texts_diff[i].get("1.0", "end-1c")
            if ode_rhs == "":
                break
            f.write(ode_rhs)
            f.write("\n")
        f.write("\n")
        for i in range(self.maximal_dimension):
            boundary_lhs = self.texts_edge[i].get("1.0", "end-1c")
            if boundary_lhs == "":
                break
            f.write(boundary_lhs)
            f.write("\n")
        f.write("\n")
        f.write(self.init_a.get("1.0", "end-1c"))
        f.write("\n")
        f.write(self.init_b.get("1.0", "end-1c"))
        f.write("\n")
        f.write(self.init_value.get("1.0", "end-1c"))
        f.write("\n")
        f.write(self.init_time.get("1.0", "end-1c"))

    def save_solution(self):
        solution_dict = {"t": self.t_res}
        for i in range(self.n_):
            solution_dict[f"x_{i + 1}"] = self.x_res[i]
        df = DataFrame(solution_dict)
        f = asksaveasfile(
            initialfile='untitled.csv',
            defaultextension=".csv",
            filetypes=[("All Files", "*.*"), ("Text Documents", "*.txt"), ("Table documents", "*.csv")]
        )
        df.to_csv(path_or_buf=f)

    def calc_functional(self):
        expr = self.integral_txt.get("1.0", "end-1c")
        if not expr:
            return
        fun = parse_expr(expr, transformations=TRANSFORMATIONS, evaluate=True)
        functional = self.solver.calculate_functional(fun, self.t_res, self.x_res)
        self.integral_value.config(text=str(functional))

    def __frame_ode(self, frame: Frame):
        get_frame(frame=frame, text="Введите систему дифференциальных уравнений:", side=TOP)
        for i in range(self.maximal_dimension):
            fr = Frame(frame)
            fr.pack(side=TOP)
            self.texts_diff.append(Text(fr, width=30, height=1, font=("Arial", 15)))
            Label(fr, text=f"dx_{i + 1}/dt =", font=("Inter", 15)).pack(side=LEFT, expand=True, pady=5, padx=40)
            self.texts_diff[i].pack(side=LEFT, fill="x", expand=True)

    def __frame_boundary_values(self, frame: Frame):
        frame = get_frame(frame=frame, text="Введите систему краевых условий:", side=TOP)
        for i in range(self.maximal_dimension):
            fr = Frame(frame)
            fr.pack(side=TOP)
            self.texts_edge.append(Text(fr, width=30, height=1, font=("Arial", 15)))
            self.texts_edge[i].pack(side=LEFT, fill="x", expand=True)
            Label(fr, text="= 0", font=("Inter", 15)).pack(side=LEFT, expand=True, pady=5, padx=40)

    def __draw(self, dr1, dr2):
        conv = {f"x_{i + 1}": self.x_res[i] for i in range(self.n_)}
        conv["t"] = self.t_res
        plot = Tk()
        plot.title("График")
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

    def __frame_axes(self, frame: Frame):
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
        self.__frame_ode(frame_left)
        self.init_a, self.init_b = frame_time_interval(frame_left)
        self.init_time = frame_time_star(frame_left)
        self.inner_method = frame_inner_method(frame_left)
        self.external_method = frame_external_method(frame_left)
        frame_right = Frame(self.window)
        frame_right.pack(side=RIGHT)
        self.__frame_boundary_values(frame_right)
        self.init_value = frame_initial_value(frame_right)
        return frame_right

    def execute(self):
        frame = self.create_interface()
        self.progressbar = Progressbar(frame, orient="horizontal", length=200, value=0)
        Button(frame, text="Решить задачу", font=("Arial", 15), command=self.exec).pack(side=TOP, fill="y", pady=5)
        Label(frame, text="Алгоритм выполняется...", font=("Inter", 15)).pack(side=TOP, expand=True, pady=5, padx=40)
        self.progressbar.pack(side=TOP)

        self.__frame_axes(frame)

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
