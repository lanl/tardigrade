"""
===============================================================================
|                               process_data.py                               |
===============================================================================
| A collection of scripts to read data from exodus files using the yt library.|
| Note: A good way to generate the functions is to:                           |
|       from functools import partial                                         |
|       fxns = [partial(f,*args,**keywords) for args,keywords in ...]         |
|                                                                             |
|       This makes independent functions which can be given different args    |
|       and keywords.                                                         |
===============================================================================
"""

#Import modules
import yt
import numpy as np
import functools

def get_dof_coordinate_data(data_set, dof, meshname = 'connect1'):
    """
    =================================
    |    get_dof_coordinate_data    |
    =================================

    Get the degree of freedom data and the
    x,y,z coordinates at which that data is defined.

    note that dof can either be a string or a list of strings.
    """

    if type(dof)==str:
        dof_data = [[float(d) for d in data] for data in data_set.all_data()[meshname, dof]]
    elif type(dof)==list:
        _dof_data = [[[float(d) for d in data] for data in data_set.all_data()[meshname, _dof]] for _dof in dof]
        dof_data = [zip(*v) for v in zip(*_dof_data)]

    coord_x  = data_set.all_data()[meshname, 'vertex_x']
    coord_y  = data_set.all_data()[meshname, 'vertex_y']
    coord_z  = data_set.all_data()[meshname, 'vertex_z']

    return [[(d,float(_x),float(_y),float(_z)) for d,_x,_y,_z in zip(data,x,y,z)] for data,x,y,z in zip(dof_data, coord_x, coord_y, coord_z)]

def evaluate_functions_at_coordinates(list_of_functions,coordinates):
    """
    ===========================================
    |    evaluate_functions_at_coordinates    |
    ===========================================

    Evaluate functions at the given coordinate points.
    Functions should be of the form v = f(x,y,z)
    """

    return [tuple([f(x,y,z) for f in list_of_functions]) for x,y,z in coordinates]

def evaluate_manufactured_solution_result_at_step(filename,dof,list_of_functions,step=-1,meshname='connect1',rtol=1e-5):
    """
    =======================================================
    |    evaluate_manufactured_solution_result_at_step    |
    =======================================================

    Evaluate the manufactured solution result and 
    return a true-false statement if the solution 
    has converged at a given step. Defaults to the 
    last step of the simulation.
    """

    data_set                = yt.load(filename,step=-1)
    simulation_results      = get_dof_coordinate_data(data_set, dof, meshname = meshname);
    flat_simulation_results = [item for sublist in simulation_results for item in sublist]
    manufactured_solution   = evaluate_functions_at_coordinates(list_of_functions,[sr[1:] for sr in flat_simulation_results])

    if(len(flat_simulation_results) != len(manufactured_solution)):
        print("Error: there aren't as many simulation results as manufactured solutions")
        return False
    result = all([np.allclose(a,b,rtol=rtol) for a,b in zip([r[0] for r in flat_simulation_results],manufactured_solution)])

    if(result==False):
        print("Result failed. Computing maximum differences...\n")
        diffs = np.array([np.array(a)-np.array(b) for a,b in zip([r[0] for r in flat_simulation_results],manufactured_solution)])
        print("maximum abs differences: {0}".format(np.max(np.abs(diffs),axis=0)))
        

    return result
    

### UTILITY FUNCTIONS ###

def const_fxn(x,y,z,v=0.):
    return v

def linear_fxn(x,y,z,a=0.,b=0.,c=0.,d=0.):
    return a*x + b*y + c*z + d

def generate_linear_functions(n,bounds=[1.0,-1.0],seed=123):
    """
    ===================================
    |    generate_linear_functions    |
    ===================================

    Generate linear functions with random coefficients
    """

    np.random.seed(seed)
    coefs = [(bounds[1] - bounds[0])*np.random.rand(4) + bounds[0] for _ in range(n)]
    strings = ["{0}*x+{1}*y+{2}*z+{3}".format(coef[0],coef[1],coef[2],coef[3]) for coef in coefs]
    return [functools.partial(linear_fxn,a=coef[0],b=coef[1],c=coef[2],d=coef[3]) for coef in coefs],coefs,strings

def generate_random_phis(stretch_scale_bounds = [0.5,2.0], theta_bounds = [-0.5*np.pi,0.5*np.pi],seed=123):
    """
    ==============================
    |    generate_random_phis    |
    ==============================

    Generate random values of phi that will be physical.
    """

    np.random.seed(seed)
    S = np.diag([(b-a)*np.random.rand(3) + a])*np.eye(3)

    thetas = (theta_bounds[1]-theta_bounds[0])*np.random.rand(3) + theta_bounds[0]

    Rx = np.array([[                 1,                 0,                  0],\
                   [                 0, np.cos(thetas[0]), -np.sin(thetas[0])],\
                   [                 0, np.sin(thetas[0]),  np.cos(thetas[0])]])

    Ry = np.array([[ np.cos(thetas[1]),                 0,  np.sin(thetas[1])],\
                   [                 0,                 1,                  0],\
                   [-np.sin(thetas[1]),                 0,  np.cos(thetas[1])]])

    Rz = np.array([[ np.cos(thetas[2]),-np.sin(thetas[2]),                  0],\
                   [ np.sin(thetas[2]), np.cos(thetas[2]),                  0],\
                   [                 0,                 0,                  1]])
    
    phi = Rx*Ry*Rz*S

def rotate_matrix(A,thetas):
    """
    =======================
    |    rotate_matrix    |
    =======================

    Rotate the given matrix by the 
    provided angles. The order with which 
    these are applied are x rotation, then y, then z.

    thetas = [theta_x, theta_y, theta_z]
    """

    Rx = np.array([[                 1,                 0,                  0],\
                   [                 0, np.cos(thetas[0]), -np.sin(thetas[0])],\
                   [                 0, np.sin(thetas[0]),  np.cos(thetas[0])]])

    Ry = np.array([[ np.cos(thetas[1]),                 0,  np.sin(thetas[1])],\
                   [                 0,                 1,                  0],\
                   [-np.sin(thetas[1]),                 0,  np.cos(thetas[1])]])

    Rz = np.array([[ np.cos(thetas[2]),-np.sin(thetas[2]),                  0],\
                   [ np.sin(thetas[2]), np.cos(thetas[2]),                  0],\
                   [                 0,                 0,                  1]])
    
    return np.dot(Rz,np.dot(Ry,np.dot(Rx,A)))

def form_linear_phi_eqns(stretch_scale_bounds = [0.5,2.0],theta_bounds = [-.5*np.pi,0.5*np.pi],seed=123):
    """
    ==============================
    |    form_linear_phi_eqns    |
    ==============================

    Form equations that result in a linear stretch of phi.
    """

    np.random.seed(seed)
    terms = [np.diag((stretch_scale_bounds[1]-stretch_scale_bounds[0])*np.random.rand(3) + stretch_scale_bounds[0])*np.eye(3) for _ in range(4)]

    thetas = (theta_bounds[1] - theta_bounds[0])*np.random.rand(3) + theta_bounds[0]

    #Compute the total rotation matrix
    R = rotate_matrix(np.eye(3),thetas)

    #Compute the rotation terms
    rotated_terms = [np.dot(R,t) for t in terms]

    #Compute the coefficients
    coefs = zip(*[(rt[0,0],rt[1,1],rt[2,2],rt[1,2],rt[0,2],rt[0,1],rt[2,1],rt[2,0],rt[1,0])\
                  for rt in rotated_terms])

    strings = ["{0}*x+{1}*y+{2}*z+{3}".format(coef[0],coef[1],coef[2],coef[3]) for coef in coefs]
    return [functools.partial(linear_fxn,a=coef[0],b=coef[1],c=coef[2],d=coef[3]) for coef in coefs],coefs,strings

