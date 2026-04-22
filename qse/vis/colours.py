import numpy as np

colors = np.sqrt(np.linspace(0.1, 0.9, 10))

qse_palette = {
    "colors": colors,
    "rads": np.linspace(12, 1, 10),
}

qse_green = (0.1, colors[5], 0.5)
qse_red = (colors[6], 0.1, 0.5)
