""" Builds Lorenz Attractors 
Author: Patrick Harrison
"""

import os
from tqdm import tqdm
import numpy as np
import random
import xml.etree.ElementTree as et

from Lorenz import LorenzAttractor

FILENAME = "LorenzData.xml"

def randomize_params(states=(0, 10), sigma=(0, 20), beta=(0, 10), rho=(20, 30)):
    """ returns parameters for a lorenz attractor """
    
    initialState = (random.uniform(states[0],states[1]),
                    random.uniform(states[0],states[1]),
                    random.uniform(states[0],states[1])
                   )
    
    sigma = random.uniform(sigma[0], sigma[1])
    beta = random.uniform(beta[0], beta[1])
    rho = random.uniform(rho[0], rho[1])
    
    return initialState, sigma, beta, rho

def create_file():
    """ Creates a new xml file for storing lorenz data """
    root = et.Element("LorenzData")
    tree = et.ElementTree(root)

    with open (FILENAME, "wb") as files :
        tree.write(files)
    
    return None

def element_writer(initialState, sigma, beta, rho):
    """ Creates a xml subelement from an attractor with 
        subelements being the XYZ states
    """
    t = np.arange(0.0, 30.0, 0.001)

    LA = LorenzAttractor(t, initialState, sigma, beta, rho)   
    X, Y, Z = LA.get_states()

    attributes = {"initialState": str(LA.get_initialState()), 
                  "sigma": str(LA.get_sigma()), 
                  "beta": str(LA.get_beta()), 
                  "rho": str(LA.get_rho())
                 }
    Attractor = et.Element("Attractor", attrib=attributes)

    x = et.SubElement(Attractor, "X")
    y = et.SubElement(Attractor, "Y")
    z = et.SubElement(Attractor, "Z")
    

    x.text = ','.join([str(elem) for elem in X])
    y.text = ','.join([str(elem) for elem in Y])
    z.text = ','.join([str(elem) for elem in Z])
    
    return Attractor
    

def append_data(initialState, sigma, beta, rho):
    """ Calculates the data and writes it to xml """

    tree = et.parse(FILENAME)
    root = tree.getroot()
    a = element_writer(initialState, sigma, beta, rho)
    root.append(a)

    # write to file
    with open (FILENAME, "wb") as files :
        tree.write(files)
        files.close()

    del(tree)
    del(root)
    del(a)

if __name__ == "__main__":
   
    print("starting")
    # check if file exists
    if os.path.isfile(FILENAME):
        print("File exists")
    else:
        print("Creating new file:", FILENAME)
        create_file()

    print("Writing Data")
    for attractor in tqdm(range(1024)):
        initialState, sigma, beta, rho = randomize_params(states=(0, 10), 
                                                          sigma=(0, 20), 
                                                          beta=(0, 10), 
                                                          rho=(20, 30)
                                                         )
        print("initialState:", initialState, 
              "sigma:", sigma, 
              "beta", beta, 
              "rho", rho,
              end='\r'
             )

        append_data(initialState, sigma, beta, rho)

    print("--- done ---")