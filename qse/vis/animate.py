import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

viridis = mpl.colormaps["viridis"]


def animate_qbits(qbits, quantity):

    # Set up the figure and axis
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_aspect("equal")
    ax.axis("off")  # Hide axes

    for i in range(qbits.nqbits - 1):
        for j in range(i + 1, qbits.nqbits):
            # if laplacian[i,j] == -1:
            ax.plot(
                [qbits[i].x, qbits[j].x],
                [qbits[i].y, qbits[j].y],
                zorder=-1,
                c="k",
                alpha=0.3,
            )

    circles = []
    for q in qbits:
        circle = plt.Circle((q.x, q.y), 0.1, color=viridis(quantity[0, q.index]))
        ax.add_patch(circle)
        circles.append(circle)

    # Animation update function
    def update(frame):
        for q in qbits:
            circles[q.index].set_color(viridis(quantity[frame, q.index]))
        return circles

    # Create the animation
    ani = FuncAnimation(
        fig,
        update,
        frames=len(quantity),
        interval=200,  # 50ms between frames
        blit=True,
    )
    return ani
