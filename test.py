import multiprocessing

from watson.watson import Watson

Watson('./KIC9944201/', './KIC9944201/').vetting("KIC 9944201", 0.721526, 131.71, 91, 21.730, 'all', 0.01,
                                                 ra=288.49090, dec=46.80500,
                                       cadence='long', cpus=multiprocessing.cpu_count() // 2, clean=False)