import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure(figsize=(5.8, 2.4))
n = 20
font_sizes = np.linspace(180, 1, n)[::-1]
colors = np.sqrt(np.linspace(0.1, 0.9, n))[::-1]
alphas = np.linspace(0, 1, n)

for font_size, c, alpha in zip(font_sizes, colors, alphas):
    plt.text(x=-0.06, y=0.1, s="QSE", size=font_size, color=(0.1, c, 0.5), alpha=alpha)

plt.axis("off")
plt.savefig("logo.svg", format="svg", transparent=True)
# plt.show()
