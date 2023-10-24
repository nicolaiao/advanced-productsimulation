import BeamModels_with_TODO as models
import SolverAlgorithms_with_TODO as algs
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

num_nodes = 10
#model = models.SimplySupportedBeamModel(num_nodes)
model = models.CantileverWithEndMoment(num_nodes)

#algs.solveLinearSteps(model,load_steps=0.01,max_steps=100)
#algs.solveNonlinLoadControl(model,load_steps=0.01,max_steps=100)
algs.solveArchLength(model,archLength=0.01,max_steps=100)

num_steps = len(model.load_history)

for iStep in range(num_steps):
   print("LoadFactor= {:12.3e}".format(model.load_history[iStep]))
   print("dispVec={:}".format(iStep))
   print(model.disp_history[iStep])

def animate_model(model, step_inc):
    num_steps = len(model.load_history)

    # Create a figure
    fig, (ax, ax_shape) = plt.subplots(nrows=1, ncols=2, figsize=(20, 10))

    def update(iStep):
        model.plotDispState(iStep, fig=fig, ax=ax, ax_shape=ax_shape)  # Pass the figure and axes

    # Create an animation object
    ani = FuncAnimation(fig, update, frames=range(0, num_steps, step_inc), repeat=False)

    plt.show()

# Create matplotlib plots
step_inc = (num_steps // 10)
animate_model(model, step_inc)
#for iStep in range(0,len(model.load_history), step_inc):
    #model.plotDispState(iStep)
    # Create an animation object
    #animate_model(model, step_inc)


# write vtu-files for animation in ParaView
for iStep in range(num_steps):
    model.vtu_print_state(iStep)
print("End")
