import numpy as np
from matplotlib import pyplot as plt


def plot_scores(best_scores, worst_scores):
    # fig, ax1 = plt.subplots()
    #
    # # Plotting best fitness scores
    # color = 'tab:blue'
    # ax1.set_xlabel('Generation')
    # ax1.set_ylabel('Best Fitness Score', color=color)
    # ln1 = ax1.plot(self.best_scores, label='Best Fitness Score', color=color)
    #
    # # Plotting worst fitness scores
    # ax2 = ax1.twinx()
    # color = 'tab:red'
    # ax2.set_ylabel('Worst Fitness Score', color=color)
    # ln2 = ax2.plot(self.worst_scores, label='Worst Fitness Score', color=color)
    #
    # # Adding a single scatter plot for the last scores
    # last_gen = len(self.best_scores) - 1
    # ln3 = ax1.scatter(last_gen, self.best_scores[-1], color='purple',
    #                   label=f'Last Best Score: {self.best_scores[-1]:.2f}\nLast Worst Score: {self.worst_scores[-1]:.2f}')
    # ax2.scatter(last_gen, self.worst_scores[-1],
    #             color='purple')
    #
    # # Combining all legends into one
    # lns = ln1 + ln2 + [ln3]
    # labels = [l.get_label() for l in lns]
    # ax1.legend(lns, labels, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
    #
    # fig.tight_layout(rect=[0, 0, 1, 0.95])  # Adjust layout to leave space for the title
    # fig.suptitle('Fitness Scores Through Generations', y=0.98)
    #
    # plt.show()

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Plot for Best Fitness Score
    axes[0].plot(best_scores, label='Best Fitness Score', color='blue')
    axes[0].scatter(len(best_scores) - 1, best_scores[-1], color='blue',
                    label=f'Last Score: {best_scores[-1]:.2f}')
    axes[0].set_xlabel('Generation')
    axes[0].set_ylabel('Best Fitness Score')
    axes[0].set_title('Best Fitness Score Through Generations')
    axes[0].legend()

    # Plot for Worst Fitness Score
    axes[1].plot(worst_scores, label='Worst Fitness Score', color='red')
    axes[1].scatter(len(worst_scores) - 1, worst_scores[-1], color='red',
                    label=f'Last Score: {worst_scores[-1]:.2f}')
    axes[1].set_xlabel('Generation')
    axes[1].set_ylabel('Worst Fitness Score')
    axes[1].set_title('Worst Fitness Score Through Generations')
    axes[1].legend()

    plt.tight_layout()
    plt.show()


def create_radar_chart(target_values, final_values):
    # Labels and data
    labels = ['Molecular\nWeight', 'Isoelectric\nPoint', 'Hydrophobicity', 'Charge',
              'Aliphatic\nIndex', 'Instability\nIndex', 'Boman', 'Hydrophobic\nMoment']

    # Normalize the data
    max_values = [max(tv, fv) for tv, fv in zip(target_values, final_values)]
    normalized_target_values = [tv / mv for tv, mv in zip(target_values, max_values)]
    normalized_final_values = [fv / mv for fv, mv in zip(final_values, max_values)]

    # Compute the angle of each axis
    num_vars = len(labels)
    angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()

    # Complete the loop for the plot
    normalized_target_values += normalized_target_values[:1]
    normalized_final_values += normalized_final_values[:1]
    angles += angles[:1]

    # Plot
    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))

    ax.plot(angles, normalized_target_values, color='blue', linewidth=2, label='Target Values')
    ax.plot(angles, normalized_final_values, color='red', linewidth=2, label='Final Values')

    # Add labels to the plot
    ax.set_yticklabels([])
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels)

    # Increase radial distance of the labels
    ax.tick_params(axis='x', pad=20)  # Adjust the padding between labels and plot

    # Add value annotations on one specified axis

    start, num = calculate_axis_crossing_points(normalized_target_values, normalized_final_values)

    specified_axis = 0  # Index of the specified axis (e.g., 0 for 'Molecular Weight')
    for grid_value in np.linspace(start, 0.8, num):  # Values at 0.0, 0.2, 0.4, 0.6, 0.8, and 1.0
        ax.text(angles[specified_axis], grid_value, f'{grid_value:.1f}',
                horizontalalignment='center', verticalalignment='center', fontsize=10, color='black',
                bbox=dict(facecolor='white', edgecolor='none', pad=0.2))

    # Customize gridlines
    ax.grid(color='gray', linestyle='-', linewidth=0.5)

    # Customize the legend
    plt.legend(loc='upper right', bbox_to_anchor=(1.1, 1.1))

    # Show the plot
    plt.show()


def calculate_axis_crossing_points(normalized_target_values, normalizes_final_values):
    min_value = min(min(normalized_target_values), min(normalizes_final_values))

    if min_value >= 0.0:
        return 0.0, 5

    # Step 1: Scale the min_value by 10
    scaled_value = min_value * 10

    # Step 2: Floor the scaled value
    floored_value = round(scaled_value)

    # Step 3: Adjust to the nearest even number if needed
    if floored_value % 2 != 0:
        floored_value -= 1

    # Step 4: Scale back to the original decimal place
    start_value = floored_value / 10.0

    step_size = 0.2
    target_value = 1.0

    # Calculate the difference
    difference = target_value - start_value

    # Calculate the number of steps
    steps = difference / step_size

    if start_value < -0.4:
        steps = 5

    return start_value, int(steps) + 1
