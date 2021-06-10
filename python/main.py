#!/usr/bin/env python3
import contextlib
import time
import sys
import os
from collections import defaultdict
import pulp as pl
import pyhull.convex_hull as hull
import logging
import math
import json
import itertools
from typing import *

logging.basicConfig(level=logging.DEBUG)


def die(msg):
    print(msg, file=sys.stderr)
    exit(1)


def det(p, q, r):
    sum1 = q[0] * r[1] + p[0] * q[1] + r[0] * p[1]
    sum2 = q[0] * p[1] + r[0] * q[1] + p[0] * r[1]
    return sum1 - sum2


def is_left_turn(p, q, r):
    assert p != q and q != r and p != r
    return det(p, q, r) > 0


class Timer(object):
    def __init__(self):
        self.start_time = None
        self.elapsed = 0
        self.count = 0
        self.last = None
        self.snap_elapsed = None

    def start(self):
        if self.start_time is not None:
            raise ValueError
        self.start_time = time.time()
        return self.start_time

    def stop(self):
        if self.start_time is None:
            raise ValueError
        elapsed = time.time() - self.start_time
        self.last = elapsed
        self.elapsed += elapsed
        self.count += 1
        self.start_time = None

    def snap(self):
        if self.start_time is not None:
            raise ValueError
        self.snap_elapsed = self.elapsed

    def since_snap(self):
        if self.start_time is not None:
            raise ValueError
        return self.elapsed - self.snap_elapsed

    def __enter__(self):
        self.start()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.stop()

    def total(self):
        if self.start_time is not None:
            raise ValueError
        return self.elapsed

    def stats(self):
        if self.start_time is not None:
            raise ValueError
        return {
            "count": self.count,
            "total": self.elapsed
        }

    def is_running(self):
        return self.start_time is not None

    def stop_if_running(self):
        if self.start_time is not None:
            self.stop()


class ParetoSet(object):
    @staticmethod
    def combine(weighted, combine_timer: Timer = contextlib.suppress(), hull_timer: Timer = contextlib.suppress()):
        with combine_timer:
            points: Dict[float, float] = dict()

            for distribution in weighted:
                for combination in itertools.product(*[[(weight, point) for point in pareto_set.points]
                                                       for weight, pareto_set in distribution]):
                    probability, expectation = 0.0, 0.0
                    for weight, point in combination:
                        probability += weight * point[0]
                        expectation += weight * point[1]

                    if expectation < points.get(probability, math.inf):
                        points[probability] = expectation
            point_list = [(p, e) for p, e in points.items()]

        with hull_timer:
            return ParetoSet(point_list)

    def __init__(self, points: Collection[Tuple[float, float]]):
        if not points:
            raise ValueError

        for point in points:
            if len(point) != 2:
                raise ValueError
        self.points: List[Tuple[float, float]]

        rightmost_point = max(points, key=lambda x: (x[0], -x[1]))
        if rightmost_point[0] == 0.0:
            self.points = [rightmost_point]
        elif rightmost_point == (1.0, 0.0):
            self.points = [(0.0, rightmost_point[1]), rightmost_point]
        else:
            pruned_points = list(filter(lambda x: x[1] <= rightmost_point[1], points))
            if len(pruned_points) <= 2:
                self.points = pruned_points
            else:
                top_left_point = 0.0, rightmost_point[1] + 1.0
                top_right_point = rightmost_point[1], rightmost_point[1] + 1.0
                pruned_points.extend([top_left_point, top_right_point])
                self.points = [pruned_points[v[0]] for v in hull.ConvexHull(pruned_points).vertices
                               if v[0] < len(pruned_points) - 2]

    def __len__(self):
        return len(self.points)

    def __str__(self):
        return "[" + ", ".join(map(lambda x: f"{x[0]:f}@{x[1]:f}",
                                   sorted(self.points, key=lambda x: x[0]))) + "]"

    def best_expectation(self, p):
        previous = None
        for point in sorted(self.points, key=lambda x: x[0], reverse=True):
            if point[0] <= p:
                if previous is None:
                    return math.inf
                p_left, e_left = point
                p_right, e_right = previous
                assert p_left < p_right, f"{p_left:f} {p_right:f}"
                assert p_left <= p <= p_right, f"{p_left:f} {p:f} {p_right:f}"
                return e_left + (e_right - e_left) * (p - p_left) / (p_right - p_left)
            else:
                previous = point
        raise AssertionError


if __name__ == '__main__':
    if len(sys.argv) < 6:
        die(f"Usage: {sys.argv[0]} <modes: vi,lp> <prism file> <goal label> <cost file> <threshold> <output>")

    modes, filename, goal_label, cost_file, threshold, output_file = sys.argv[1:]
    modes = modes.split(",")

    parse_timer = Timer()
    with parse_timer:
        if filename.endswith("."):
            filename = filename[:-1]
        if threshold.startswith("range"):
            start, end, count = threshold.split(",")[1:4]
            start, end, count = float(start), float(end), int(count)
            thresholds = list(start + i * (end - start) / (count - 1) for i in range(count))
            assert len(thresholds) == count
        else:
            thresholds = [float(threshold)]
        logging.debug("Considering thresholds: %s", ", ".join(f"{t:.3f}" for t in thresholds))

        for suffix in ["tra", "lab"]:
            if not os.path.isfile(filename + "." + suffix):
                die(f"Could not find file {os.path.basename(filename)}.{suffix}")
        if cost_file != "reach" and not os.path.isfile(cost_file):
            die(f"Could not find file {cost_file}")
        for mode in modes:
            if mode not in ["vi", "lp"]:
                die("Mode must be vi or lp")

        goal_states = set()
        absorbing_states = set()
        with open(file=filename + ".lab", mode="rt") as lab:
            first = lab.readline().strip()
            initial_index, goal_index = None, None
            for label in first.split(" "):
                index, label_name = label.split("=")
                label_name = label_name[1:-1]
                if label_name == "init":
                    if initial_index is not None:
                        die("Multiple initial labels")
                    initial_index = int(index)
                elif label_name == goal_label:
                    if goal_index is not None:
                        die("Multiple goal labels")
                    goal_index = int(index)
            if initial_index is None:
                die("No initial state label")

            initial_state = None
            for row in lab:
                state_row = row.strip().split(" ")
                state = state_row[0][:-1]
                for index in state_row[1:]:
                    if int(index) == initial_index:
                        if initial_state is not None:
                            die("Multiple initial states")
                        initial_state = int(state)
                    if int(index) == goal_index:
                        goal_states.add(int(state))

        if not goal_index and not goal_states:
            logging.info("No goal states found from goal label \"%s\"", goal_label)

        with open(file=filename + ".tra", mode="rt") as tra:
            metadata_line = tra.readline()
            metadata = metadata_line.split(" ")
            if not len(metadata) == 3:
                die(f"Invalid metadata line: {metadata_line} (is this an MDP?)")
            states = int(metadata[0])
            mdp = [defaultdict(dict) for _ in range(states)]

            for row in tra.readlines():
                if not row:
                    continue
                transition_data = row.split(" ")
                if len(transition_data) <= 3 or len(transition_data) > 5:
                    die(f"Invalid transition line: {row} (is this an MDP?)")
                source, action, dest, prob = transition_data[:4]
                source, dest, prob = int(source), int(dest), float(prob)
                if source in goal_states or not prob:
                    continue
                mdp[source][int(action)][dest] = prob

        for state, actions in enumerate(mdp):
            if state not in goal_states and all(list(dist.keys()) == [state] for dist in actions.values()):
                absorbing_states.add(state)

        if absorbing_states:
            logging.warning("Adding %d absorbing states to goal states", len(absorbing_states))
            goal_states.update(absorbing_states)

        if not goal_states:
            die("No goal states")
        if initial_state in goal_states:
            die("Initial state in goal states")

        non_goal_states = set(state for state in range(states) if state not in goal_states)

        if cost_file != "reach":
            state_costs = dict()

            with open(file=filename + ".srew", mode="rt") as rew:
                rew.readline()
                for row in rew.readlines():
                    if not row:
                        continue
                    reward_data = row.split(" ")
                    if len(reward_data) != 2:
                        die(f"Invalid state reward line: {row}")
                    state, reward = int(reward_data[0]), int(reward_data[1])
                    if reward < 0:
                        die(f"Negative reward in {state}")
                    if not reward:
                        continue
                    state_costs[state] = reward

            if not state_costs:
                die("No costs")
        else:
            state_costs = {state: 1 for state in non_goal_states}

        fixed_goal_states = 0
        for state in goal_states:
            if state_costs.get(state):
                del state_costs[state]
                fixed_goal_states += 1
        if fixed_goal_states:
            logging.warning("Removed rewards from %d goal states", fixed_goal_states)
        for state in non_goal_states:
            if state not in state_costs:
                die(f"Non-goal state {state} has no cost")

    maximal_cost = max(state_costs.values())

    non_probabilistic_states, deterministic_states, non_probabilistic_actions = 0, 0, 0
    for state, state_actions in enumerate(mdp):
        if len(state_actions) <= 1:
            deterministic_states += 1
            if not state_actions:
                state_actions[0][state] = 1.0

        probabilistic_transition = False
        for action, transitions in state_actions.items():
            if len(transitions) <= 1:
                non_probabilistic_actions += 1
                if not transitions:
                    transitions[state] = 1.0
            else:
                probabilistic_transition = True
            if math.fabs(sum(transitions.values()) - 1.0) > 0.001:
                die(f"Transition probabilities of {state} do not sum up to 1: {sum(transitions.values()):.3f}")
        if not probabilistic_transition:
            non_probabilistic_states += 1

    logging.debug("Found %d goal states, %d cost assignments (maximum: %d), %d / %d non-prob. states / actions, "
                  "%d deterministic states", len(goal_states), len(state_costs), maximal_cost,
                  non_probabilistic_states, non_probabilistic_actions, deterministic_states)

    output = {
        "name": os.path.basename(filename),
        "commandline": sys.argv,
        "data": {
            "states": states,
            "initial-state": initial_state,
            "non-prob-states": non_probabilistic_states,
            "non-prob-actions": non_probabilistic_actions,
            "det-states": deterministic_states
        },
        "thresholds": [f"{threshold:.3f}" for threshold in thresholds],
        "result": {},
        "time": {
            "parse": parse_timer.stats()
        }
    }

    solvers = pl.listSolvers(onlyAvailable=True)
    logging.debug("Available solvers: %s", ", ".join(solvers))
    if 'CPLEX_PY' in solvers:
        solver = pl.CPLEX_PY(msg=False)
    elif 'GUROBI' in solvers:
        solver = pl.GUROBI(msg=False)
    elif 'CPLEX_CMD' in solvers:
        solver = pl.CPLEX_CMD(msg=False)
    else:
        solver = pl.getSolver(solvers[0], msg=False)
    logging.debug("Using solver %s", solver.name)

    solve_timer = Timer()
    with solve_timer:
        logging.debug("Computing optimal SSP vi LP")
        prob = pl.LpProblem("SSP", pl.LpMaximize)
        state_value = []
        for state in range(states):
            state_value.append(pl.LpVariable(f"x_{state}", lowBound=0))

        prob += pl.lpSum(state_value), "Minimize expected cost"

        for state in range(states):
            if state in goal_states:
                prob += state_value[state] == 0
                continue
            for action, successors in mdp[state].items():
                prob += state_value[state] <= state_costs.get(state, 0) + \
                        pl.lpSum(state_value[successor] * probability
                                 for successor, probability in successors.items())

        prob.solve(solver)
        if prob.status != pl.LpStatusOptimal:
            die("Infeasible base problem")
        optimal_ssp = pl.value(state_value[initial_state])
        logging.info("Optimal SSP: %f", optimal_ssp)

        state_expected_cost = list(map(pl.value, state_value))
        for state, cost in enumerate(state_expected_cost):
            if cost is None or math.isnan(cost) or cost < 0:
                die(f"Invalid expected cost {cost} in state {state}")
            if cost < state_costs.get(state, 0):
                die(f"Inconsistent expected costs in state {state}")

    output["result"]["ssp"] = optimal_ssp
    output["time"]["ssp"] = solve_timer.stats()

    thresholds.sort(reverse=True)
    if "vi" in modes:
        best_cvar: List[float] = [math.inf] * len(thresholds)
        best_var: List[Optional[int]] = [None] * len(thresholds)
        threshold_timers = [Timer() for _ in thresholds]

        cvar_timer = Timer()
        hull_timer = Timer()
        combine_timer = Timer()
        for timer in threshold_timers:
            timer.start()

        with cvar_timer:
            current_sets: List[ParetoSet] = [
                ParetoSet([(1.0, 0.0) if state in goal_states else (0.0, state_expected_cost[state])])
                for state in range(states)]
            step = 1
            while step <= best_cvar[-1]:
                logging.debug("VI iteration %d; %d points, %d non-trivial sets", step,
                              sum(map(len, current_sets)), sum(1 for pareto in current_sets if len(pareto.points) > 2))

                hull_timer.snap()
                combine_timer.snap()
                current_sets = [ParetoSet.combine(itertools.chain([(probability, current_sets[successor])
                                                                   for successor, probability in transitions.items()]
                                                                  for action, transitions in mdp[state].items()),
                                                  hull_timer=hull_timer, combine_timer=combine_timer)
                                for state in range(states)]
                logging.debug("Step %d: %f combine, %f hull", step, combine_timer.since_snap(), hull_timer.since_snap())

                for index, threshold in enumerate(thresholds):
                    cvar = step + current_sets[initial_state].best_expectation(1.0 - threshold) / threshold

                    if cvar < optimal_ssp:
                        die(f"CVaR {cvar} better than SSP {optimal_ssp}")

                    if cvar < best_cvar[index]:
                        logging.info("Improving best CVaR from %f to %f for %.3f", best_cvar[index], cvar, threshold)
                        best_cvar[index] = cvar
                        best_var[index] = step

                    if step > best_cvar[index]:
                        threshold_timers[index].stop_if_running()
                step += 1

            for timer in threshold_timers:
                timer.stop_if_running()

        output["time"]["vi"] = {
            "total": cvar_timer.stats(),
            "hull": hull_timer.stats(),
            "combine": combine_timer.stats(),
            "thresholds": {
                f"{threshold:.3f}": threshold_timers[index].stats() for index, threshold in enumerate(thresholds)
            }
        }
        output["result"]["vi"] = {
            f"{threshold:.3f}": {
                "cvar": best_cvar[index],
                "var": best_var[index]
            } for index, threshold in enumerate(thresholds)
        }

    if "lp" in modes:
        best_cvar: List[float] = [math.inf] * len(thresholds)
        best_var: List[Optional[int]] = [None] * len(thresholds)
        threshold_timers = [Timer() for _ in thresholds]

        cvar_timer = Timer()
        cvar_building_timer = Timer()
        cvar_solving_timer = Timer()

        with cvar_timer:
            for index, threshold in enumerate(thresholds):
                logging.info("Solving threshold %.3f", threshold)
                threshold_timers[index].start()
                var_guess, solution = 1, None

                while var_guess <= best_cvar[index]:
                    with cvar_building_timer:
                        logging.debug("Solving for guess %d (current best: %f)", var_guess, best_cvar[index])

                        cvar_prob = pl.LpProblem("SSP_CVaR", pl.LpMinimize)
                        cvar_state_prob = [dict() for _ in range(var_guess)]
                        cvar_action_prob = [dict() for _ in range(var_guess)]
                        cvar_cost_transitions = defaultdict(list)
                        cvar_goal_prob = [pl.LpVariable(f"goal_{i}") for i in range(var_guess + 1)]

                        for i in range(var_guess):
                            for state in non_goal_states:
                                cvar_state_prob[i][state] = pl.LpVariable(f"x_{i}_{state}", lowBound=0)
                                for action in mdp[state].keys():
                                    cvar_action_prob[i][(state, action)] = \
                                        pl.LpVariable(f"x_{i}_{state}_{action}", lowBound=0)

                        cvar_prob += cvar_goal_prob[0] == 0
                        for state in non_goal_states:
                            cvar_prob += cvar_state_prob[0][state] == (1 if state == initial_state else 0)

                        for i in range(var_guess):
                            for state in non_goal_states:
                                cvar_prob += cvar_state_prob[i][state] == pl.lpSum(cvar_action_prob[i][(state, action)]
                                                                                   for action in mdp[state].keys())

                        cvar_transition_prob = [defaultdict(list) for i in range(var_guess)]
                        cvar_goal_transition_prob = [[] for i in range(var_guess + 1)]

                        for i in range(var_guess):
                            for state in non_goal_states:
                                state_cost = state_costs[state]
                                for action, transitions in mdp[state].items():
                                    action_prob = cvar_action_prob[i][(state, action)]
                                    for successor, probability in transitions.items():
                                        if i + state_cost < var_guess:
                                            if successor in goal_states:
                                                target = cvar_goal_transition_prob[i + state_cost]
                                            else:
                                                target = cvar_transition_prob[i + state_cost][successor]
                                        else:
                                            target = cvar_cost_transitions[i + state_cost - var_guess +
                                                                           state_expected_cost[successor]]
                                        target.append(probability * action_prob)

                        for i, cost_transitions in enumerate(cvar_transition_prob):
                            for state, state_transitions in cost_transitions.items():
                                cvar_prob += cvar_state_prob[i][state] == pl.lpSum(state_transitions)
                        for i in range(var_guess):
                            cvar_prob += cvar_goal_prob[i + 1] == \
                                         cvar_goal_prob[i] + pl.lpSum(cvar_goal_transition_prob[i + 1])

                        cvar_prob += cvar_goal_prob[var_guess - 1] <= 1 - threshold
                        cvar_prob += cvar_goal_prob[var_guess] + cvar_cost_transitions[0] >= 1 - threshold

                        cvar_prob += pl.lpSum(cost * pl.lpSum(variables)
                                              for cost, variables in cvar_cost_transitions.items()), "Minimize CVaR"

                    logging.debug("Building took %f, running solver", cvar_building_timer.last)
                    with cvar_solving_timer:
                        cvar_prob.solve(solver)
                    if cvar_prob.status == pl.LpStatusOptimal:
                        value = pl.value(cvar_prob.objective)
                        cvar = var_guess + value / threshold
                        logging.debug("Got CVaR %f for VaR %d after %f", cvar, var_guess, cvar_solving_timer.last)

                        # Debugging
                        #
                        # for v in cvar_prob.variables():
                        #     if not v.varValue:
                        #         continue
                        #     print(v.name, "=", v.varValue)

                        if cvar < optimal_ssp:
                            die(f"CVaR {cvar} better than SSP {optimal_ssp}")

                        if cvar < best_cvar[index]:
                            logging.info("Improving best CVaR from %f to %f for %.3f",
                                         best_cvar[index], cvar, threshold)
                            best_cvar[index] = cvar
                            best_var[index] = var_guess

                            # solution = {state: [dict() for _ in range(var_guess)] for state in non_goal_states
                            #             if len(mdp[state]) > 1}
                            # for state in non_goal_states:
                            #     if len(mdp[state]) == 1:
                            #         continue
                            #     for i in range(var_guess):
                            #         state_value = pl.value(cvar_state_prob[i][state])
                            #         if state_value < 1e-10:
                            #             solution[state][i] = None
                            #             continue
                            #         for action in mdp[state].keys():
                            #             action_value = pl.value(cvar_action_prob[i][(state, action)])
                            #             if action_value < 1e-10:
                            #                 continue
                            #             solution[state][i][action] = action_value / state_value
                    else:
                        logging.debug("Got got non-solution after %f: %s",
                                      cvar_solving_timer.last, pl.LpStatus[cvar_prob.status])

                    var_guess += 1

                threshold_timers[index].stop()

                # compressed_strategy = defaultdict(dict)
                #
                # for state, state_strategy in solution.items():
                #     step = 0
                #     while step < len(state_strategy):
                #         current_strategy = state_strategy[step]
                #         while current_strategy is None and step + 1 < len(state_strategy):
                #             step += 1
                #             current_strategy = state_strategy[step]
                #
                #         start = step
                #         while step < len(state_strategy):
                #             next_strategy = state_strategy[step]
                #             if next_strategy is not None:
                #                 if next_strategy.keys() != current_strategy.keys():
                #                     break
                #                 if any(math.fabs(next_strategy[action] - current_strategy[action]) > 1.0001
                #                        for action in current_strategy.keys()):
                #                     break
                #             step += 1
                #         if start == 0 and not current_strategy or current_strategy is None:
                #             continue
                #         compressed_strategy[state][start] = current_strategy
                #
                # for state in sorted(compressed_strategy.keys()):
                #     print(state)
                #     for step, strategy in compressed_strategy[state].items():
                #         print(f"  {step}: {strategy}")

        output["result"]["lp"] = {
            f"{threshold:.3f}": {
                "cvar": best_cvar[index],
                "var": best_var[index]
            } for index, threshold in enumerate(thresholds)
        }
        output["time"]["lp"] = {
            "total": cvar_timer.stats(),
            "building": cvar_building_timer.stats(),
            "solving": cvar_solving_timer.stats(),
            "thresholds": {
                f"{threshold:.3f}": threshold_timers[index].stats() for index, threshold in enumerate(thresholds)
            }
        }

    with open(output_file, mode="wt") as f:
        json.dump(output, f)
