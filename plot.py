import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import interact

class FunctionPlotter:
    def f(self, func, x):
        return eval(func)

    def plot_function(self, func):
        x = np.linspace(-10, 10, 1000)
        y = self.f(func, x)

        plt.figure(figsize=(8, 6))
        plt.plot(x, y)
        plt.title("Plot of " + func)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.grid(True)
        plt.show()

fp = FunctionPlotter()

interact(fp.plot_function, func="x**2")