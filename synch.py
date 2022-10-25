import xgi 
import numpy as np


def EM():
    p_ij = 1
    p_ijk = 1
    e_i = 1
    
    
def p_ij():
    pass
    
def p_ijk():
    pass

def e_i():
    pass

def rho_j():    
    pass
    
def rho_jk():
    pass

def rho_ei():
    pass


def theta_timeseries(H, k2, k3, w, theta, timesteps=10000, dt=0.002):
    """This function calculates the order parameter for the Kuramoto model on hypergraphs.
    This solves the Kuramoto model ODE on hypergraphs with edges of sizes 2 and 3
    using the Euler Method. It returns an order parameter which is a measure of synchrony.

    Parameters
    ----------
    H : Hypergraph object
        The hypergraph on which you run the Kuramoto model
    k2 : float
        The coupling strength for links
    k3 : float
        The coupling strength for triangles
    w : numpy array of real values
        The natural frequency of the nodes.
    theta : numpy array of real values
        The initial phase distribution of nodes.
    timesteps : int greater than 1, default: 10000
        The number of timesteps for Euler Method.
    dt : float greater than 0, default: 0.002
        The size of timesteps for Euler Method.

    Returns
    -------
    r_time : numpy array of floats
        timeseries for Kuramoto model order parameter

    References
    ----------
    "Synchronization of phase oscillators on complex hypergraphs"
    by Sabina Adhikari, Juan G. Restrepo and Per Sebastian Skardal
    https://doi.org/10.48550/arXiv.2208.00909

    Examples
    --------
    >>> import numpy as np
    >>> import xgi
    >>> n = 100
    >>> H = xgi.random_hypergraph(n, [0.05, 0.001], seed=None)
    >>> w = 2*np.ones(n)
    >>> theta = np.linspace(0, 2*np.pi, n)
    >>> R = compute_kuramoto_order_parameter(H, k2 = 2, k3 = 3, w = w, theta = theta)
    """
    thetas = [theta]
    
    H_int = xgi.convert_labels_to_integers(H, "label")

    links = H_int.edges.filterby("size", 2).members()
    triangles = H_int.edges.filterby("size", 3).members()
    n = H_int.num_nodes

    r_time = np.zeros(timesteps)

    for t in range(timesteps):

        r1 = np.zeros(n, dtype=complex)
        r2 = np.zeros(n, dtype=complex)

        for i, j in links:

            r1[i] += np.exp(1j * theta[j])
            r1[j] += np.exp(1j * theta[i])

        for i, j, k in triangles:

            r2[i] += np.exp(2j * theta[j] - 1j * theta[k]) + np.exp(
                2j * theta[k] - 1j * theta[j]
            )
            r2[j] += np.exp(2j * theta[i] - 1j * theta[k]) + np.exp(
                2j * theta[k] - 1j * theta[i]
            )
            r2[k] += np.exp(2j * theta[i] - 1j * theta[j]) + np.exp(
                2j * theta[j] - 1j * theta[i]
            )

        d_theta = (
            w
            + k2 * np.multiply(r1, np.exp(-1j * theta)).imag
            + k3 * np.multiply(r2, np.exp(-1j * theta)).imag
        )
        theta_new = theta + d_theta * dt
        theta = theta_new
        
        
        thetas.append(theta)
        
    return np.array(thetas)



#%% GET THETAS


n = 100
H = xgi.random_hypergraph(n, [0.05, 0.001], seed=None) #Make random hypergraph

w = 2*np.ones(n) #Natural frequency of oscillators
theta_i = np.linspace(0, 2*np.pi, n) #Initial phases


'''
If the phase of oscillator i is thetai(t), then 
thetas is a 2d array [ [theta1(0),theta2(0),...], 
                      [theta1(dt),theta2(dt),...] ,
                      [theta1(2dt),theta2(2dt),...],...]
'''
thetas = theta_timeseries(H, k2 = 2, k3 = 3, w = w, theta = theta_i)

thetas = np.mod(thetas,2*np.pi) #Wrap theta around


#%% CREATE 'INCOMPLETENESS' IN SPACE


'''
Replace fraction p_missing of node data with 0s
Nodes are chosen randomly to be replaced

e.g if random process with chooses only 2nd oscillator to be removed, then
thetas --> [ [theta1(0),0,...], 
             [theta1(dt),0,...] ,
             [theta1(2dt),0,...],...]
'''

spacial_incompleteness = True
p_missing = 0.1
missing_nodes = np.random.randint(0,n-1,size=round(p_missing*n))
if spacial_incompleteness:
    thetas[:,missing_nodes]=0

#%% DISCRETISE IF NEEDED


n_bins = 3
thetas[2*np.pi/n_bins > thetas ] == 1
thetas[4*np.pi/n_bins > thetas >= 2*np.pi/n_bins] == 2
thetas[thetas>= 4*np.pi/n_bins] == 3



#%% Calculate P for decreases in theta_i 
N = len(H.edges)

P_ij = np.zeros((N,N))
P_ijk = np.zeros((N,N,N))
        
for i,theta_i in enumerate(thetas.T):
    for j,theta_j in enumerate(np.delete(thetas.T, i, 0)):
        gt_ind = np.where(theta_i>theta_j)[0]
        gt_ind = gt_ind[gt_ind<len(theta_i)-1]
        if len(gt_ind) == 0:
            break
        else:
            P_ij[i][j] = len(np.where(theta_i[gt_ind+1]<theta_i[gt_ind])[0])/len(gt_ind)

        
for i,theta_i in enumerate(thetas.T):
    for j,theta_j in enumerate(np.delete(thetas.T, i, 0)):
        for k,theta_k in enumerate(np.delete(thetas.T,[i,j],0)):
            gt_ind = np.where( (theta_i>theta_j) & (theta_i>theta_k) )[0]
            gt_ind = gt_ind[gt_ind<len(theta_i)-1]
            if len(gt_ind) == 0:
                break
            else:
                P_ijk[i][j][k] = len(np.where(theta_i[gt_ind+1]<theta_i[gt_ind])[0])/len(gt_ind)
    


