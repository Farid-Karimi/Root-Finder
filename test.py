import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import matplotlib
matplotlib.use('TkAgg')


# Get the function from the user
func_str = input("Enter the function to plot (e.g., 'x**2'): ")
func = lambda x: eval(func_str)


# Create the figure and axis
fig, ax = plt.subplots(figsize=(8, 6))
plt.subplots_adjust(bottom=0.25)
plt.style.use('Solarize_Light2')

# Plot the function
x = np.linspace(-100, 100, 1000)
line, = ax.plot(x, func(x))

# Show the grid
ax.grid(True)


# Enable interactive features
fig.canvas.mpl_connect('scroll_event', lambda event: plt.axis('auto' if abs(event.delta) > 0 else 'tight'))
fig.canvas.mpl_connect('key_press_event', lambda event: [fig.canvas.key_press_event(event) for fig in plt.get_figures()])


plt.show()

