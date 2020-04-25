from epydemic import *
import networkx
import pandas as pd

class SIR(CompartmentedModel):
    SUSCEPTIBLE = 'S'
    INFECTED = 'I'
    REMOVED = 'R'

    P_INFECTED = 'pInfected'
    P_INFECT = 'pInfect'
    P_REMOVE = 'pRemove'

    SI = 'SI'

    def build(self, params):
        pInfected = params[self.P_INFECTED]
        pInfect = params[self.P_INFECT]
        pRemove = params[self.P_REMOVE]
        
        self.addCompartment(self.INFECTED, pInfected)
        self.addCompartment(self.REMOVED, 0.0)
        self.addCompartment(self.SUSCEPTIBLE, 1 - pInfected)

        self.trackNodesInCompartment(self.INFECTED)
        self.trackEdgesBetweenCompartments(self.SUSCEPTIBLE, self.INFECTED, name=self.SI)

        self.addEventPerElement(self.SI, pInfect, self.infect)
        self.addEventPerElement(self.INFECTED, pRemove, self.remove)

    def infect( self, t, e ):
        (n, m) = e
        self.changeCompartment(n, self.INFECTED)
        self.markOccupied(e, t)

    def remove( self, t, n ):
        self.changeCompartment(n, self.REMOVED)




def main():

    param = dict()
    param[SIR.P_INFECT] = 0.1
    param[SIR.P_REMOVE] = 0.5
    param[SIR.P_INFECTED] = 0.01

    N = 10000                 # order (number of nodes) of the network
    kmean = 5                 # mean node degree
    phi = (kmean + 0.0) / N   # probability of attachment between two nodes chosen at random

    # create the network
    g = networkx.erdos_renyi_graph(N, phi)

    # create a model and a dynamics to run it
    m = SIR()                      # the model (process) to simulate
    e = StochasticDynamics(m, g)   # use stochastic (Gillespie) dynamics

    # set the parameters we want and run the simulation
    rc = e.set(param).run()

    pass

if __name__ == "__main__":
    main()