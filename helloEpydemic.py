import epydemic
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


class MonitoredSIR(epydemic.SIR):
    INTERVAL = 'interval'
    PROGRESS = 'progress'

    def setUp(self, params):
        """Schedule the monitoring event.

        :param params: the simulation parameters"""
        super(MonitoredSIR, self).setUp(params)

        # add a monitoring event to fill-in the evolution of the process
        self._series = []
        self.postRepeatingEvent(0, params[self.INTERVAL], None, self.monitor)

    def monitor(self, t, e):
        """Record the sizes of each compartment.

        :param t: the simulation time
        :param e: the element (ignored)"""
        s = dict()
        for k in [epydemic.SIR.SUSCEPTIBLE, epydemic.SIR.INFECTED, epydemic.SIR.REMOVED]:
            s[k] = len(self.compartment(k))
        self._series.append((t, s))

    def results(self):
        """Add the time series to the experiment's results.

        :returns: a results dict including the monitored time series"""
        rc = super(MonitoredSIR, self).results()

        rc[self.PROGRESS] = self._series
        return rc


def main():
    # use an ER network as the substrate
    N = 1000
    kmean = 3
    phi = kmean / N
    g = nx.erdos_renyi_graph(N, phi)
    options = {
        'node_color': 'black',
        'node_size': 5,
        'width': 10,
    }

    plt.subplot(121)
    nx.draw(g, **options)
    # plt.subplot(122)
    # nx.draw_shell(g, nlist=[range(5, 10), range(5)], with_labels=True, font_weight='bold')

    plt.show()

    # set the parameters the same as above
    params = dict()
    params[epydemic.SIR.P_INFECT] = 0.02  # infection probability
    params[epydemic.SIR.P_REMOVE] = 0.002  # recovery probability
    params[epydemic.SIR.P_INFECTED] = 0.01  # initial fraction infected

    # capture every 10s
    params[MonitoredSIR.INTERVAL] = 10

    e = epydemic.StochasticDynamics(MonitoredSIR(), g=g)
    e.process().setMaximumTime(1000)
    rc = e.set(params).run()

    plt.plot([x[1]['S'] for x in rc['results']['progress']])
    plt.plot([x[1]['I'] for x in rc['results']['progress']])
    plt.plot([x[1]['R'] for x in rc['results']['progress']])

    plt.show()
    pass


if __name__ == "__main__":
    main()
