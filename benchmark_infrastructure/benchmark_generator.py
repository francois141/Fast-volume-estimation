import numpy as np
from configparser import ConfigParser
import argparse
import glob
from datetime import datetime
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# parser = argparse.ArgumentParser(description='Process input data')
# parser.add_argument('--configuration', type=str, required=True)

# args = parser.parse_args()

config = ConfigParser()
config.read('configurations/computer_config_constantin.ini')
# config.read('configurations/computer_config_francois_home.ini')

# data = ConfigParser()
# data.read(args.configuration)


title = "{} @ {} <br>L1d cache: {} <br>L1i cache: {} <br>L2 cache: {} <br>L3 cache: {} <br>Compiler:  {} Flags: {}".format(
    config.get('cpu', 'type'),
    config.get('cpu', 'frequence'),
    config.get('cache', 'l1d'),
    config.get('cache', 'l1i'),
    config.get('cache', 'l2'),
    config.get('cache', 'l3'),
    config.get('compiler', 'version'),
    config.get('compiler', 'flags'),
)


def draw_values():
    with open('./files_to_plot/sizes.txt') as f:
        sizes = f.readlines()

    files = glob.glob('./files_to_plot/cycles_*.txt')

    fig = make_subplots()
    fig.update_layout(barmode='group',
                  title=dict(
                      text=title,
                      x=0.05,
                      y=0.95,
                      xanchor='left',
                      yanchor='top',
                      font=dict(
                          size=20,
                          color='#000000'
                      )
                  ))
    fig.update_layout(font=dict(size=20))
    fig.update_xaxes(title_text="Input size M (2M vectors of size M)")
    fig.update_yaxes(title_text="Cycles", secondary_y=False)
    fig.update_layout(legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01))
    # ticks
    # fig.update_layout(xaxis=dict(tickmode="linear", tick0=0, dtick=1))
    fig.update_layout(yaxis=dict(tickmode="linear", tick0=0, dtick=1e10))
    fig.update_layout(margin=dict(t=200))
    def plot_values(x, y, label):
        fig.add_trace(go.Scatter(x=x, y=y, name=label))

    for file in files:
        with open(file, 'r') as f:
            cycles = f.readlines()

        cycles = np.array(cycles).astype('float64')
        plot_values(sizes, cycles, file)
    fig.show()

def save():
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    current_date = now.strftime("%Y-%m-%d")
    # print("Current Time =", current_time)
    # print("Current Date =", current_date)


draw_values()
# save()
