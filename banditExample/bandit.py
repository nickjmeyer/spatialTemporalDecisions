"""Bandit example to show why Thompson Sampling is necessary"""

import numpy as np

import matplotlib.pyplot as plt

import multiprocessing as mp

import cPickle as pickle

import traceback

import pandas

class Bandit(object):
    def __init__(self, handle_means = (0, 0)):
        assert isinstance(handle_means, tuple)

        assert len(handle_means) == 2

        self.handle_means_ = handle_means

        self.values_ = None

    def set_values(self):
        self.values_ = [np.random.normal(loc = self.handle_means_[i],
                                         scale = 1.0)
                        for i in range(2)]


    def pull_handle(self, arm_number):
        ## pull handle and return the value
        assert isinstance(arm_number, int)
        assert (arm_number == 0 or arm_number == 1)

        return self.values_[arm_number]


    def regret(self, value):
        return max(self.values_) - value


class KstepAlternatingAgent(object):
    """An agent that explores by alternating for k steps then picks the best"""

    def __init__(self, ksteps):
        assert isinstance(ksteps, int)

        self.ksteps_ = ksteps

        self.nsteps_ = 0

        self.history_ = [[], []]

        self.offset_ = np.random.binomial(1, 0.5)


    def choose_handle(self):
        if self.nsteps_ < self.ksteps_:
            ## if still eploring, return alternating handle
            return (self.nsteps_ + self.offset_) % 2
        else:
            ## if done exploring, return the best or choose randomly if tied
            if np.mean(self.history_[0]) > np.mean(self.history_[1]):
                return 0
            elif np.mean(self.history_[0]) < np.mean(self.history_[1]):
                return 1
            else:
                return np.random.binomial(1,0.5)


    def feedback(self, value):
        ## record history
        if self.nsteps_ < self.ksteps_:
            self.history_[(self.nsteps_ + self.offset_) % 2].append(value)


    def turn_clock(self):
        ## increment steps
        self.nsteps_ += 1



class TsExploreAgent(object):
    """An agent that uses Thompson Sampling to balance exploration and
    exploitation

    """

    def __init__(self):

        self.nsteps_ = 0

        self.history_ = [[], []]

        self.last_handle_ = None


    def choose_handle(self):
        beliefs = self.draw_beliefs()

        if beliefs[0] > beliefs[1]:
            ## if belief for 0 is larger than 1, return 0
            self.last_handle_ = 0
            return 0
        elif beliefs[0] < beliefs[1]:
            ## if belief for 0 is smaller than 1, return 1
            self.last_handle_ = 1
            return 1
        else:
            ## if beliefs are tied, choose randomly
            self.last_handle_ = np.random.binomial(1, 0.5)
            return self.last_handle_


    def draw_beliefs(self):
        beliefs = []
        for i in range(2):
            ## calculate mean and standard error for each handle and
            ## draw a belief
            handle_mean = None
            handle_se = None

            n_obs = len(self.history_[i])

            if n_obs == 0:
                ## use prior
                handle_mean = 0.0
                handle_se = 1.0
            else:
                ## use posterior
                handle_mean = sum(self.history_[i]) / (1. / 100. + n_obs)
                handle_se = np.sqrt(1. / (1. / 100. + n_obs))


            assert handle_mean is not None
            assert handle_se is not None

            ## draw a belief
            beliefs.append(np.random.normal(loc = handle_mean,
                                            scale = handle_se))

        return beliefs


    def feedback(self, value):
        ## record history
        self.history_[self.last_handle_].append(value)


    def turn_clock(self):
        ## increment steps
        self.nsteps_ += 1



class Runner(object):
    def __init__(self, bandit, agent, max_steps):
        self.bandit_ = bandit
        self.agent_ = agent

        self.pull_history_ = []
        self.value_history_ = []
        self.regret_history_ = []

        self.nsteps_ = 0
        self.max_steps_ = max_steps


    def turn_clock(self):
        ## bandit selects values for each arm
        self.bandit_.set_values()

        ## get handle from the agent
        handle = self.agent_.choose_handle()
        self.pull_history_.append(handle)

        ## pull the handle and get the value
        value = self.bandit_.pull_handle(handle)
        self.value_history_.append(value)
        self.regret_history_.append(self.bandit_.regret(value))

        ## provide feedback to the agent
        self.agent_.feedback(value)

        ## increment clock
        self.agent_.turn_clock()


    def run(self):
        ## run for max_steps
        for steps in range(self.max_steps_):
            self.turn_clock()

        ## return histories
        return self.pull_history_, self.value_history_, self.regret_history_


def create_agent(kind, *args):
    if kind == "TsExploreAgent":
        return TsExploreAgent(*args)
    elif kind == "KstepAlternatingAgent":
        return KstepAlternatingAgent(*args)
    else:
        raise ValueError("Don't know how to create agent of type: "
                         + str(kind))


def run_experiment_for_agent(a_name, a_args, num_reps,
                             handle_means, final_t):
    res = {"pull": [],
           "value": [],
           "regret": [],
           "time": [],
           "rep": [],
           "name": []}

    a_res_pull = []
    a_res_value = []
    for rep in range(num_reps):
        np.random.seed(rep)

        agent = create_agent(*a_args)
        bandit = Bandit(handle_means)
        runner = Runner(bandit, agent, final_t)

        pull_history, value_history, regret_history = runner.run()

        assert len(pull_history) == final_t
        assert len(value_history) == final_t

        res["pull"].extend(pull_history)
        res["value"].extend(value_history)
        res["regret"].extend(regret_history)
        res["time"].extend(range(1, final_t + 1))

        res["rep"].extend([rep] * final_t)
        res["name"].extend([a_name] * final_t)


    return pandas.DataFrame.from_dict(res)


def wrapper(args):
    return run_experiment_for_agent(*args)


def run_experiment():
    num_reps = 1000

    final_t = 100

    handle_means = (1, 0)

    ## create arguments for map
    arg_list = []
    for i in [1, 5, 10]:
        arg_list.append(("%d-Step Alternating" % (2 * i),
                         ("KstepAlternatingAgent", 2 * i),
                         num_reps, handle_means, final_t))

    arg_list.append(("Thompson Sampling", ("TsExploreAgent",),
                     num_reps, handle_means, final_t))


    ## run sims
    df = pandas.concat(map(wrapper, arg_list))

    ## save data
    df.to_csv("all_res.csv", index = False)


if __name__ == "__main__":
    run_experiment()
