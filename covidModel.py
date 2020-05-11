from epydemic import *
import networkx
import pandas as pd
import numpy as np
import matplotlib
import random
import matplotlib.pyplot as plt

from covidCharacteristics import get_random_days_till_death


death_curve_for_run = [0] * 365

death_day_average_for_run = {}

class SIR(CompartmentedModel):
    SUSCEPTIBLE = 'S'
    INFECTED = 'I'
    REMOVED = 'R'
    DEAD = 'D'

    P_INFECTED = 'pInfected'
    P_INFECT = 'pInfect'
    P_DEAD = 'pDead'
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

    def infect(self, t, e):
        (n, m) = e
        day = int(t * 10)
        # if 17 < day < 20:
        if 17 <= day < 20:
            if random.random() > 0.66:
                pass
                return
        if 20 <= day < 50:
            if random.random() > 0.2:
                pass
                return

        self.changeCompartment(n, self.INFECTED)
        self.markOccupied(e, t)

        d = get_random_days_till_death()
        if d > 0:
            death_curve_for_run[day + d] = death_curve_for_run[day + d] + 1
            if (day + d) in death_day_average_for_run:
                x = death_day_average_for_run[day+d][0] * death_day_average_for_run[day+d][1]
                death_day_average_for_run[day + d][1] = death_day_average_for_run[day + d][1] + 1
                death_day_average_for_run[day + d][0] = (x + d) / death_day_average_for_run[day+d][1]

            else:
                death_day_average_for_run[day+d] = [d,1]

    def remove(self, t, n):
        self.changeCompartment(n, self.REMOVED)


def main():
    param = dict()
    param[SIR.P_INFECT] = 0.3
    param[SIR.P_REMOVE] = 0.05
    param[SIR.P_INFECTED] = 0.01

    N = 50000  # order (number of nodes) of the network

    # days_till_death_distribution = [0] * 54
    #
    # for n in range(N):
    #     d = get_random_days_till_death()
    #     if d > 0:
    #         days_till_death_distribution[d] = days_till_death_distribution[d] + 1
    #
    # plt.subplot(121)
    # plt.plot(days_till_death_distribution)
    # plt.show()


    kmean = 5  # mean node degree
    phi = (kmean + 0.0) / N  # probability of attachment between two nodes chosen at random

    # create the network
    g = networkx.erdos_renyi_graph(N, phi)

    # create a model and a dynamics to run it
    m = SIR()  # the model (process) to simulate
    e = StochasticDynamics(m, g)  # use stochastic (Gillespie) dynamics

    # set the parameters we want and run the simulation
    rc = e.set(param).run()

    fig, ax = plt.subplots()
    ax.set_xlabel('Day')
    ax.set_ylabel('Deaths per day; average days-to-death')
    ax.set_title('Deaths per day: lockdown days 20-50; lesser mitigation day 17')
    major_ticks = np.arange(0, 121, 10)
    minor_ticks = np.arange(0, 121, 5)
    major_yticks = np.arange(0, 131, 10)
    minor_yticks = np.arange(0, 131, 5)
    ax.set_xticks(major_ticks)
    ax.set_xticks(minor_ticks, minor=True)
    ax.set_yticks(major_yticks)
    ax.set_yticks(minor_yticks, minor=True)
    ax.set_ylim([0, 130])
    plt.grid(which='both')
    plt.plot(death_curve_for_run[:120])
    death_day_average_array = [0]*365
    for entry in death_day_average_for_run:
        death_day_average_array[entry] = death_day_average_for_run[entry][0]
    plt.plot(death_day_average_array[:120])
    plt.show()

    pass


if __name__ == "__main__":
    main()