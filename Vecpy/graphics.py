import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np

#General
def quiver(X,Y,U,V,title,units):
    plt.figure()
    plt.title(title)
    Q = plt.quiver(X,Y,U,V,np.sqrt(U**2+V**2),angles='xy',cmap='coolwarm')
    qk = plt.quiverkey(Q, 0.9, 0.9, 1, r'1 '+units, labelpos='E',coordinates='figure')
    plt.show()

#Frame
def frame_quiver(frame):
    X,Y = frame.return_grid()
    U,V = frame.return_all_velocities()
    title = 'Quiver of frame: '+ frame.properties.frame
    units = frame.properties.velocity_units
    quiver(X,Y,U,V,title,units)
    

def frame_tke_graph():
    plt.tricontourf(X,Y,tke)

#Run
def run_quiver_animation(run):
    fig, ax = plt.subplots(figsize=(21, 18))
    if run.check_same_grid_run():
        X,Y = run.run_grid()
        frames = run.frames()
        shape = (X.shape[0],X.shape[1],len(frames))
        U = np.zeros(shape)
        V = np.zeros(shape)
        for ind in range(len(frames)):                
            x,y,u,v = run.fields[frames[ind]].create_mesh_grid()
            U[:,:,ind] = u[::-1]
            V[:,:,ind] = v[::-1]
    
    else:
        print('Not all run in same grid')
        return
    
    Q = ax.quiver(X,Y,U[...,0],V[...,0],np.sqrt(U[...,0]**2+V[...,0]**2),angles='xy',cmap='coolwarm')
    
    def update_quiver(i,Q,U,V):
        Q.set_UVC(U[...,i],V[...,i])
        return Q,
    
    anim = animation.FuncAnimation(fig, update_quiver,frames=shape[2],interval=200,fargs=(Q,U,V), blit=False,repeat=True)
    fig.tight_layout()
    plt.show()
    return anim
    

def mean_run_quiver(run):
    X,Y = run.run_grid()
    U,V = run.run_mean_velocities()
    title = 'Quiver of mean velocity over all frames'
    units = run.fields[run.frames()[0]].properties.velocity_units
    quiver(X,Y,U,V,title,units)

