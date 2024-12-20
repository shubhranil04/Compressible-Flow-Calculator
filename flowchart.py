import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# Function to create a fancy box with text inside
def fancy_box(ax, xy, text):
    ax.add_patch(FancyBboxPatch(xy, width=2.5, height=0.5, boxstyle="round,pad=0.1", edgecolor="black", facecolor="lightblue"))
    ax.text(xy[0] + 1.25, xy[1] + 0.25, text, ha="center", va="center", fontsize=10)

# Set up the figure
fig, ax = plt.subplots(figsize=(10, 10))
ax.set_xlim(-3, 7)
ax.set_ylim(-8, 3)
ax.axis('off')

# Creating boxes for each step
fancy_box(ax, (0, 2), "Start Program")
fancy_box(ax, (0, 1), "User Input")
fancy_box(ax, (-4, 0), "1. Isentropic Flow Relations")
fancy_box(ax, (0, 0), "2. Normal Shock Relations")
fancy_box(ax, (4, 0), "3. Oblique Shock Relations")

fancy_box(ax, (2, -2), "Check for Normal Shock")
fancy_box(ax, (2, -4), "Normal Shock Calculations")

fancy_box(ax, (-2, -6), "Check for Oblique Shock")
fancy_box(ax, (-2, -8), "Oblique Shock Calculations")

fancy_box(ax, (2, -6), "Handle Errors or Warnings")
fancy_box(ax, (2, -8), "Display Results")

# Add arrows to show flow
ax.annotate("", xy=(-2, 1.5), xytext=(-2, 2.5), arrowprops=dict(arrowstyle="->"))
ax.annotate("", xy=(-2, -0.5), xytext=(-2, 0.5), arrowprops=dict(arrowstyle="->"))
ax.annotate("", xy=(-2, -2.5), xytext=(-2, -1.5), arrowprops=dict(arrowstyle="->"))
ax.annotate("", xy=(-2, -4.5), xytext=(-2, -3.5), arrowprops=dict(arrowstyle="->"))

ax.annotate("", xy=(2, -2.5), xytext=(-2, -2.5), arrowprops=dict(arrowstyle="->"))
ax.annotate("", xy=(2, -4.5), xytext=(2, -3.5), arrowprops=dict(arrowstyle="->"))

ax.annotate("", xy=(-2, -6.5), xytext=(-2, -5.5), arrowprops=dict(arrowstyle="->"))
ax.annotate("", xy=(-2, -8.5), xytext=(-2, -7.5), arrowprops=dict(arrowstyle="->"))

ax.annotate("", xy=(2, -6.5), xytext=(2, -5.5), arrowprops=dict(arrowstyle="->"))
ax.annotate("", xy=(2, -8.5), xytext=(2, -7.5), arrowprops=dict(arrowstyle="->"))

# Display the figure
plt.show()

