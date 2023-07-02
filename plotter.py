import matplotlib.pyplot as plt


class Plotter:
    def __init__(self, t, x):
        self.t_ = t
        self.x_ = x

    def plot_x(self):
        plt.xlabel('t')
        plt.ylabel('x(t)')
        plt.title('x(t) = (x0(t), x1(t), ..., x_n-1(t))')
        colors = {0: 'red', 1: 'blue', 2: 'green', 3: 'orange'}
        plt.plot(self.x_[0], self.x_[1])
        # plt.grid()
        # for i in range(len(self.x_)):
        #     plt.plot(self.t_, self.x_[i], color=colors[i], label=f'x{i}(t)')
        legend = plt.legend(loc='upper left', shadow=True, fontsize='x-large')
        legend.get_frame()
        plt.show()
