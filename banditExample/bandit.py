"""Bandit example to show why Thompson Sampling is necessary"""

import numpy as np

import matplotlib.pyplot as plt

import multiprocessing as mp

import cPickle as pickle

import traceback

class Bandit(object):
    def __init__(self, handle_means = (0, 0), handle_sds = (1, 1)):
        assert isinstance(handle_means, tuple)
        assert isinstance(handle_sds, tuple)

        assert len(handle_means) == 2
        assert len(handle_sds) == 2

        self.handle_sds_ = handle_sds
        self.handle_means_ = handle_means


    def pull_handle(self, arm_number):
        ## pull handle and return the value
        assert isinstance(arm_number, int)
        assert (arm_number == 0 or arm_number == 1)

        return np.random.normal(loc = self.handle_means_[arm_number],
                                scale = self.handle_sds_[arm_number])


class KstepAlternatingAgent(object):
    """An agent that explores by alternating for k steps then picks the best"""

    def __init__(self, ksteps):
        assert isinstance(ksteps, int)

        self.ksteps_ = ksteps

        self.nsteps_ = 0

        self.history_ = [[], []]


    def choose_handle(self):
        if self.nsteps_ < self.ksteps_:
            ## if still eploring, return alternating handle
            return self.nsteps_ % 2
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
        self.history_[self.nsteps_ % 2].append(value)


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

        self.default_mean_ = 0
        self.default_se_ = 1


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

            if len(self.history_[i]) == 0:
                ## use defaults when no data are available
                handle_mean = self.default_mean_
                handle_se = self.default_se_
            elif len(self.history_[i]) == 1:
                ## use default se when only one observation is available
                handle_mean = self.history_[i][0]
                handle_se = self.default_se_
            else:
                ## use sample mean and sample se
                handle_mean = np.mean(self.history_[i])
                handle_se = (np.std(self.history_[i])
                             / np.sqrt(len(self.history_[i])))

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

        self.nsteps_ = 0
        self.max_steps_ = max_steps


    def turn_clock(self):
        ## get handle from the agent
        handle = self.agent_.choose_handle()
        self.pull_history_.append(handle)

        ## pull the handle and get the value
        value = self.bandit_.pull_handle(handle)
        self.value_history_.append(value)

        ## provide feedback to the agent
        self.agent_.feedback(value)

        ## increment clock
        self.agent_.turn_clock()


    def run(self):
        ## run for max_steps
        for steps in range(self.max_steps_):
            self.turn_clock()

        ## return histories
        return self.pull_history_, self.value_history_


def create_agent(kind, *args):
    if kind == "TsExploreAgent":
        return TsExploreAgent(*args)
    elif kind == "KstepAlternatingAgent":
        return KstepAlternatingAgent(*args)
    else:
        raise ValueError("Don't know how to create agent of type: "
                         + str(kind))


def run_experiment_for_agent(a_args, num_reps,
                             handle_means, handle_sds,
                             final_t):
    a_res_pull = []
    a_res_value = []
    for rep in range(num_reps):
        agent = create_agent(*a_args)
        bandit = Bandit(handle_means, handle_sds)
        runner = Runner(bandit, agent, final_t)

        pull_history, value_history = runner.run()

        a_res_pull.append(np.mean([pull == 0 for pull in pull_history]))

        a_res_value.append(np.mean(value_history))

    return a_res_pull, a_res_value


def wrapper(args):
    try :
        return run_experiment_for_agent(*args)
    except Exception as e:
        print traceback.format_exc()
        print e
        raise e



def run_experiment():
    num_reps = 100

    final_t = 100

    handle_means = (2, 0)
    handle_sds = (5, 1)

    ## create arguments for map
    arg_list = []
    for i in range(50):
        arg_list.append((("KstepAlternatingAgent", 2 * (i + 1)),
                         num_reps, handle_means, handle_sds, final_t))

    arg_list.append((("TsExploreAgent",),
                         num_reps, handle_means, handle_sds, final_t))


    ## run sims
    pool = mp.Pool()

    all_res = pool.map(wrapper, arg_list)


    ## save data
    with open("all_res.p", "w") as f:
        f.write(pickle.dumps(all_res))


if __name__ == "__main__":
    run_experiment()
