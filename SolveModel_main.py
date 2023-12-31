import BeamModels_with_TODO as models
import SolverAlgorithms_with_TODO as algs
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pdb

def animate_model(model, step_inc):
    num_steps = len(model.load_history)

    # Create a figure
    fig, (ax, ax_shape) = plt.subplots(nrows=1, ncols=2, figsize=(20, 10))

    def update(iStep):
        model.plotDispState(iStep, fig=fig, ax=ax, ax_shape=ax_shape)  # Pass the figure and axes

    # Create an animation object
    ani = FuncAnimation(fig, update, frames=range(0, num_steps, step_inc), repeat=False)

    plt.show()

num_nodes = 21 # NEED TO BE ODD NUMBERS OF NODES (To get the force in the midle of the beam)
#model = models.SimplySupportedBeamModel(num_nodes)
#model = models.CantileverWithEndMoment(num_nodes)

#algs.solveLinearSteps(model,load_steps=0.01,max_steps=100)
#algs.solveNonlinLoadControl(model,load_steps=0.01,max_steps=100)
#model = models.DeepArchModel(num_nodes)
models = [models.SimplySupportedBeamModel(num_nodes), 
          models.CantileverWithEndMoment(num_nodes),
          models.DeepArchModel(num_nodes)
        ]
max_steps=1000
for model in models:
    #algs.solveLinearSteps(model,load_steps=0.01,max_steps=100)
    #algs.solveNonlinLoadControl(model,load_steps=0.01,max_steps=100)
    algs.solveArchLength(model,archLength=5.0,max_steps=max_steps)
    num_steps = len(model.load_history)
    for iStep in range(num_steps):
        print("LoadFactor= {:12.3e}".format(model.load_history[iStep]))
        #print("dispVec={:}".format(iStep))
        #print(model.disp_history[iStep])

    # Create matplotlib plots
    step_inc = max((num_steps // 10), 1)
    print(len(model.load_history), step_inc)
    animate_model(model, step_inc)
    #for iStep in range(0,len(model.load_history), step_inc):
        #model.plotDispState(iStep)
        # Create an animation object
        #animate_model(model, step_inc)


    # write vtu-files for animation in ParaView
    for iStep in range(num_steps):
        model.vtu_print_state(iStep)
    print("End")
