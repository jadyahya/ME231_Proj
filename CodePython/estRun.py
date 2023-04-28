import numpy as np
import scipy as sp
import time
#NO OTHER IMPORTS ALLOWED (However, you're allowed to import e.g. scipy.linalg)

def estRun(times, dt, internalStateIn, steeringAngle, pedalSpeed, measurement):
    # In this function you implement your estimator. The function arguments
    # are:
    #  time: current time in [s] 
    #  dt: current time step [s]
    #  internalStateIn: the estimator internal state, definition up to you. 
    #  steeringAngle: the steering angle of the bike, gamma, [rad] 
    #  pedalSpeed: the rotational speed of the pedal, omega, [rad/s] 
    #  measurement: the position measurement valid at the current time step
    #
    # Note: the measurement is a 2D vector, of x-y position measurement.
    #  The measurement sensor may fail to return data, in which case the
    #  measurement is given as NaN (not a number).
    #
    # The function has four outputs:
    #  x: your current best estimate for the bicycle's x-position
    #  y: your current best estimate for the bicycle's y-position
    #  theta: your current best estimate for the bicycle's rotation theta
    #  internalState: the estimator's internal state, in a format that can be understood by the next call to this function

    # Example code only, you'll want to heavily modify this.
    # this internal state needs to correspond to your init function:
    # x = internalStateIn[0]
    # y = internalStateIn[1]
    # theta = internalStateIn[2]

    # x = x + pedalSpeed
    # y = y + pedalSpeed
    xm0 = internalStateIn[0]
    ym0 = internalStateIn[1]
    thetam0 = internalStateIn[2]
    r = internalStateIn[3]
    b = internalStateIn[4]
    Pm0 = internalStateIn[5]



    # Prediction Step:
    xp1 = xm0+ dt*(5*r*pedalSpeed*np.cos(thetam0))
    yp1 = ym0+ dt*(5*r*pedalSpeed*np.sin(thetam0))
    thetap1 = thetam0 +dt*(5*r*pedalSpeed*np.tan(steeringAngle))/b
    bp1 = b
    rp1 = r
    xp = np.array([[xp1],[yp1],[thetap1],[rp1],[bp1]])
    
    # Matrices of Concern:
    A = np.eye(5)+dt*np.array([[0,0,-5*r*pedalSpeed*np.sin(thetam0),np.cos(thetam0)*5*pedalSpeed,0],[0,0,5*r*pedalSpeed*np.cos(thetam0),np.sin(thetam0)*5*pedalSpeed,0],[0,0,0,np.tan(steeringAngle)*5*pedalSpeed/b,-r*np.tan(steeringAngle)*5*pedalSpeed/b**2],[0,0,0,0,0],[0,0,0,0,0]])
    L = dt*np.array([[5*r*pedalSpeed*np.cos(thetam0),0, 5*r*np.cos(thetam0)],[5*r*pedalSpeed*np.sin(thetam0),0,5*r*np.sin(thetam0)],[5*r*pedalSpeed*np.tan(steeringAngle)/b,(5*1*r*pedalSpeed*(np.tan(steeringAngle)**2+1))/b,5*r*np.tan(steeringAngle)/b],[0,0,0],[0,0,0]])
    H = np.array([[1,0,-b/2*np.sin(thetap1),0,0.5*np.cos(thetap1)],[0,1,b/2*np.cos(thetap1),0,0.5*np.sin(thetap1)]])
    M = np.array([[1,0],[0,1]])

    sigvv= np.array([[0.05,0,0],[0,0.01,0],[0,0,4]])
    sigww = np.array([[40,0],[0,40]])
    # sigww = np.array([[1.09,0],[0,2.99]])
    Pp1 = A@Pm0@A.T+L@sigvv@L.T
    Kk1 = Pp1@H.T@sp.linalg.inv(H@Pp1@H.T+M@sigww@M.T)
    
    if not (np.isnan(measurement[0]) or np.isnan(measurement[1])):
        xtot = xp+Kk1@(np.array([[measurement[0]],[measurement[1]]])-np.array([[xp1+0.5*0.8*np.cos(thetap1)],[yp1+0.5*0.8*np.sin(thetap1)]]))
        Pm1 = (np.eye(5)-Kk1@H)@Pp1
        x = xtot[0][0]
        y = xtot[1][0]
        theta = xtot[2][0]
        r = xtot[3][0]
        b = xtot[4][0]
        print(r)
        print(b)
    else:
        # x = xm0
        # y = ym0
        # theta = thetam0
        # Pm1 = Pm0
        x = xp1
        y = yp1
        theta=thetap1
        r = rp1
        b = bp1
        Pm1 = Pp1




    #### OUTPUTS ####
    # Update the internal state (will be passed as an argument to the function
    # at next run), must obviously be compatible with the format of
    # internalStateIn:
    internalStateOut = [x,
                     y,
                     theta,
                     r,
                     b,
                     Pm1
                     ]

    # DO NOT MODIFY THE OUTPUT FORMAT:
    return x, y, theta, internalStateOut


