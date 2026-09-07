import matplotlib
matplotlib.use("Agg")
from solar_dynamo import Parameters
from solar_dynamo.simulation import run
from solar_dynamo.visualization import make_figure
def test_exact_grid_plot():
    result=run(Parameters(),None,2,1);fig=make_figure(result);assert fig.axes
