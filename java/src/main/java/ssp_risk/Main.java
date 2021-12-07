package ssp_risk;

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.base.Preconditions.checkState;

import com.google.common.base.Stopwatch;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.SetMultimap;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonPrimitive;
import com.google.gson.JsonSerializer;
import gurobi.GRB;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Duration;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LongSummaryStatistics;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import javax.annotation.Nullable;

public final class Main {
  private static final Logger log = Logger.getLogger("cvar");

  private Main() {}

  public static void main(String[] args) throws IOException, GRBException {
    if (args.length != 7) {
      System.out.println("Usage: <tra file> <label file> <goal label> <cost file> <modes: vi,lp> <threshold> <output>");
      System.exit(1);
    }

    String transitionFile = args[0];
    String labelsFile = args[1];
    String goalLabel = args[2];
    String costs = args[3];
    Set<String> modes = Set.of(args[4].split(","));
    Set<Double> thresholdSet = new TreeSet<>();
    String[] thresholdSpecs = args[5].split("\\|");
    for (String threshold : thresholdSpecs) {
      if (threshold.startsWith("range:")) {
        String[] rangeData = threshold.substring("range:".length()).split(",");
        double rangeStart = Double.parseDouble(rangeData[0]);
        double rangeEnd = Double.parseDouble(rangeData[1]);
        int rangeCount = Integer.parseInt(rangeData[2]);
        checkArgument(rangeCount >= 2);
        for (int i = 0; i < rangeCount; i++) {
          thresholdSet.add(rangeStart + (rangeEnd - rangeStart) / (rangeCount - 1) * i);
        }
      } else if (threshold.startsWith("logrange:")) {
        String[] rangeData = threshold.substring("logrange:".length()).split(",");
        double rangeStart = Double.parseDouble(rangeData[0]);
        double rangeEnd = Double.parseDouble(rangeData[1]);
        int rangeCount = Integer.parseInt(rangeData[2]);
        checkArgument(rangeCount >= 2);
        for (int i = 0; i < rangeCount; i++) {
          thresholdSet.add(StrictMath.pow(10, rangeStart + (rangeEnd - rangeStart) / (rangeCount - 1) * i));
        }
      } else {
        thresholdSet.add(Double.parseDouble(threshold));
      }
    }
    double[] thresholds = thresholdSet.stream().mapToDouble(Double::doubleValue).sorted().toArray();
    String outputFile = args[6];

    Result.InputData inputData = new Result.InputData(Path.of(transitionFile).getFileName().toString(), args);

    SetMultimap<String, Integer> labels = Parser.parseLabels(Path.of(labelsFile));

    MDP mdp = Parser.parseMdp(Path.of(transitionFile), labels);

    Set<Integer> goalLabelStates = labels.get(goalLabel);
    Set<Integer> absorbingStates = Parser.getAbsorbingStates(mdp);
    Set<Integer> goalStates;
    if (goalLabelStates.isEmpty() && goalLabel.equals("absorbing")) {
      goalStates = ImmutableSet.copyOf(absorbingStates);
    } else {
      goalStates = ImmutableSet.copyOf(goalLabelStates);
    }

    checkArgument(!goalStates.isEmpty(), "Did not find any goal states!");
    checkArgument(!goalStates.contains(mdp.initialState), "Initial state is goal state");

    int states = mdp.size();
    int[] stateCosts = new int[states];
    if (costs.equals("reach")) {
      for (int state = 0; state < states; state++) {
        stateCosts[state] = goalStates.contains(state) ? 0 : 1;
      }
    } else {
      Map<Integer, Integer> costMap = Parser.parseStateRewards(Path.of(costs));
      for (int state = 0; state < states; state++) {
        stateCosts[state] = goalStates.contains(state) ? 0 : costMap.getOrDefault(state, 0);
      }
    }

    Result.ModelData modelData = new Result.ModelData(states, mdp.actionCount(), mdp.transitionCount(), mdp.initialState,
        absorbingStates.size());
    Result.ProblemParameters problemParameters = new Result.ProblemParameters(goalStates.size(), thresholds);

    log.log(Level.INFO, "Problem has {0} states, {1} goal states (label {3}), {4} absorbing, solving for thresholds {5}",
        new Object[] {states, goalStates.size(), goalLabelStates.size(), goalLabel, absorbingStates.size(),
            Arrays.stream(thresholds).mapToObj("%.3f"::formatted).collect(Collectors.toList())});


    GRBEnv env = new GRBEnv();
    env.set(GRB.IntParam.LogToConsole, 0);

    log.log(Level.FINE, "Solving SSP LP");

    Stopwatch sspTimer = Stopwatch.createStarted();

    double[] stateSSP = new double[states];
    Duration sspConstructionTime;
    //noinspection UnnecessaryCodeBlock
    {
      GRBModel sspModel = new GRBModel(env);
      GRBVar[] sspStateVars = new GRBVar[states];
      for (int state = 0; state < states; state++) {
        sspStateVars[state] = sspModel.addVar(0, Double.POSITIVE_INFINITY, -1.0, GRB.CONTINUOUS, "x_" + state);
      }
      for (int state = 0; state < mdp.transitions.length; state++) {
        if (goalStates.contains(state)) {
          sspModel.addConstr(sspStateVars[state], '=', 0.0, "f_%d".formatted(state));
        } else {
          Distribution[] actions = mdp.transitions[state];
          for (int action = 0; action < actions.length; action++) {
            Distribution distribution = actions[action];
            if (distribution == null) {
              continue;
            }
            GRBLinExpr sum = new GRBLinExpr();
            sum.addConstant(stateCosts[state]);
            distribution.forEach((s, p) -> sum.addTerm(p, sspStateVars[s]));
            sspModel.addConstr(sspStateVars[state], '<', sum, "f_%d_%d".formatted(state, action));
          }
        }
      }
      sspConstructionTime = sspTimer.elapsed();

      sspModel.optimize();
      int status = sspModel.get(GRB.IntAttr.Status);
      if (status != GRB.Status.OPTIMAL) {
        log.log(Level.SEVERE, "Could not solve SSP LP, status: {0}", new Object[] {status});
        System.exit(1);
      }
      for (int state = 0; state < states; state++) {
        stateSSP[state] = sspStateVars[state].get(GRB.DoubleAttr.X);
      }
    }
    Duration sspTime = sspTimer.stop().elapsed();
    double optimalSsp = stateSSP[mdp.initialState];
    log.log(Level.INFO, "Got optimal SSP {0} after {1}", new Object[] {optimalSsp, sspTimer});

    @Nullable
    Result.ValueIterationSolution viSolution = null;
    if (modes.contains("vi")) {
      //noinspection ErrorNotRethrown
      try {
        log.log(Level.INFO, "Solving with VI");

        double[] bestCVaR = new double[thresholds.length];
        Arrays.fill(bestCVaR, Double.POSITIVE_INFINITY);
        int[] bestVaR = new int[thresholds.length];
        Arrays.fill(bestVaR, Integer.MAX_VALUE);

        Stopwatch[] thresholdStopwatches = new Stopwatch[thresholds.length];
        Arrays.setAll(thresholdStopwatches, i -> Stopwatch.createStarted());
        Stopwatch viTimer = Stopwatch.createStarted();

        ParetoSet[] sets = new ParetoSet[states];
        for (int state = 0; state < states; state++) {
          sets[state] = goalStates.contains(state)
              ? ParetoSet.of(1.0, 0.0)
              : ParetoSet.of(0.0, stateSSP[state]);
        }

        int step = 1;
        while (step <= Math.ceil(bestCVaR[0])) {
          Stopwatch iterationStopwatch = Stopwatch.createStarted();
          ParetoSet[] currentSets = sets;
          ParetoSet[] nextSets = new ParetoSet[states];
          mdp.stateStream().parallel()
              .forEach(state -> {
                if (goalStates.contains(state)) {
                  nextSets[state] = currentSets[state];
                } else {
                  ParetoSet combination = ParetoSet.combine(currentSets, mdp.transitions[state]);
                  // assert ParetoSet.checkCombination(combination.set(), currentSets, mdp.transitions[state]);
                  nextSets[state] = combination;
                }
              });
          sets = nextSets;
          if (log.isLoggable(Level.FINE)) {
            LongSummaryStatistics sizeStats = Arrays.stream(sets)
                .filter(p -> !p.isSimple())
                .mapToLong(ParetoSet::size)
                .summaryStatistics();
            log.log(Level.FINE, ("VI iteration %d took %s; %d points, %d non-trivial sets with %.2f avg, %d max points").formatted(
                step, iterationStopwatch, Arrays.stream(sets).mapToLong(ParetoSet::size).sum(),
                sizeStats.getCount(),
                sizeStats.getCount() == 0 ? 0 : sizeStats.getAverage(),
                sizeStats.getCount() == 0 ? 0 : sizeStats.getMax()));
          }

          Map<Integer, Double> improvedThresholds = new HashMap<>();
          for (int tIndex = 0; tIndex < thresholds.length; tIndex++) {
            double threshold = thresholds[tIndex];
            double expectation = sets[mdp.initialState].bestExpectation(1.0 - threshold);
            if (Double.isNaN(expectation)) {
              continue;
            }
            double cvar = step + expectation / threshold;
            checkState(Util.doublesLessOrEqual(optimalSsp, cvar));
            if (cvar < bestCVaR[tIndex]) {
              improvedThresholds.put(tIndex, bestCVaR[tIndex]);
              bestCVaR[tIndex] = cvar;
              bestVaR[tIndex] = step;
            }
            if (step >= Math.ceil(bestCVaR[tIndex])) {
              if (thresholdStopwatches[tIndex].isRunning()) {
                thresholdStopwatches[tIndex].stop();
              }
            }
          }
          if (!improvedThresholds.isEmpty() && log.isLoggable(Level.INFO)) {
            String updates = improvedThresholds.entrySet().stream()
                .sorted(Comparator.comparingDouble(e -> -thresholds[e.getKey()]))
                .map(entry -> "%.3f: %.3f -> %.3f".formatted(thresholds[entry.getKey()], entry.getValue(), bestCVaR[entry.getKey()]))
                .collect(Collectors.joining(", "));
            log.log(Level.INFO, "Guess {0} improved CVaR: {1}", new Object[] {step, updates});
          }

          step += 1;
        }

        viTimer.stop();

        Map<Double, Result.ThresholdResult> resultMap = IntStream.range(0, thresholds.length).boxed().collect(Collectors.toUnmodifiableMap(
            tIndex -> thresholds[tIndex],
            tIndex -> new Result.ThresholdResult(bestVaR[tIndex], bestCVaR[tIndex], thresholdStopwatches[tIndex].elapsed())
        ));
        viSolution = new Result.ValueIterationSolution(viTimer.elapsed(), resultMap);
      } catch (OutOfMemoryError e) {
        log.log(Level.WARNING, "Out of memory during VI", e);
      }
      //noinspection CallToSystemGC
      Runtime.getRuntime().gc();
    }

    @Nullable
    Result.LinearProgrammingSolution lpSolution = null;
    if (modes.contains("lp")) {
      //noinspection ErrorNotRethrown
      try {
        double[] bestCVaR = new double[thresholds.length];
        Arrays.fill(bestCVaR, Double.POSITIVE_INFINITY);
        int[] bestVaR = new int[thresholds.length];
        Arrays.fill(bestVaR, Integer.MAX_VALUE);

        Stopwatch[] thresholdStopwatches = new Stopwatch[thresholds.length];
        Arrays.setAll(thresholdStopwatches, i -> Stopwatch.createUnstarted());
        Stopwatch lpBuildTimer = Stopwatch.createUnstarted();
        Stopwatch lpSolveTimer = Stopwatch.createUnstarted();
        Stopwatch lpTimer = Stopwatch.createStarted();

        Stopwatch lpReach = Stopwatch.createStarted();
        int[] minimalReachabilityStep = computeMinimalReach(mdp, goalStates, thresholds);
        lpReach.stop();

        for (int tIndex = 0; tIndex < thresholds.length; tIndex++) {
          thresholdStopwatches[tIndex].start();
          double threshold = thresholds[tIndex];
          log.log(Level.INFO, "Solving for threshold {0} with LP, minimal guess is {1}",
              new Object[] {threshold, minimalReachabilityStep[tIndex]});

          for (int varGuess = minimalReachabilityStep[tIndex]; varGuess <= bestCVaR[tIndex]; varGuess += 1) {
            Stopwatch iterationTimer = Stopwatch.createStarted();
            lpBuildTimer.start();
            GRBModel cvarModel = new GRBModel(env);
            buildLinearProgrammingModel(cvarModel, threshold, varGuess, mdp, goalStates, stateSSP, stateCosts);
            if (Double.isFinite(bestCVaR[tIndex])) {
              double cutoff = (bestCVaR[tIndex] - bestVaR[tIndex]) * threshold + 1.0e-8;
              cvarModel.set(GRB.DoubleParam.Cutoff, cutoff);
            }
            // cvarModel.set(GRB.DoubleParam.BestObjStop, optimalSsp);
            lpBuildTimer.stop();

            log.log(Level.FINER, "Model for guess {0} (threshold {1}) built, calling solver",
                new Object[] {varGuess, threshold});
            lpSolveTimer.start();
            cvarModel.optimize();
            lpSolveTimer.stop();
            iterationTimer.stop();
            log.log(Level.FINE, "LP iteration took {0}", new Object[] {iterationTimer});

            int status = cvarModel.get(GRB.IntAttr.Status);
            checkState(status != GRB.Status.USER_OBJ_LIMIT);
            if (status == GRB.Status.OPTIMAL) {
              double value = cvarModel.get(GRB.DoubleAttr.ObjVal);
              double cvar = varGuess + value / threshold;
              checkState(cvar >= optimalSsp, cvar);
              if (cvar < bestCVaR[tIndex]) {
                log.log(Level.INFO, "Improving best CVaR from {0} to {1} for {2}, guess {3}",
                    new Object[] {bestCVaR[tIndex], cvar, threshold, varGuess});
                bestCVaR[tIndex] = cvar;
                bestVaR[tIndex] = varGuess;
              }
            } else {
              log.log(Level.FINE, "No solution for {0}, guess {1}, (status: {2})", new Object[] {threshold, varGuess, status});
            }
          }
          thresholdStopwatches[tIndex].stop();
        }

        lpTimer.stop();

        Map<Double, Result.ThresholdResult> resultMap = IntStream.range(0, thresholds.length).boxed().collect(Collectors.toUnmodifiableMap(
            tIndex -> thresholds[tIndex],
            tIndex -> new Result.ThresholdResult(bestVaR[tIndex], bestCVaR[tIndex], thresholdStopwatches[tIndex].elapsed())
        ));
        lpSolution = new Result.LinearProgrammingSolution(lpTimer.elapsed(), lpReach.elapsed(), lpBuildTimer.elapsed(),
            lpSolveTimer.elapsed(), resultMap);
      } catch (OutOfMemoryError e) {
        log.log(Level.WARNING, "Out of memory during LP", e);
      } catch (GRBException e) {
        log.log(Level.WARNING, "Gurobi error during LP", e);
      }
      //noinspection CallToSystemGC
      Runtime.getRuntime().gc();
    }

    if (viSolution != null && lpSolution != null) {
      for (double threshold : thresholds) {
        Result.ThresholdResult viResult = viSolution.results().get(threshold);
        Result.ThresholdResult lpResult = lpSolution.results().get(threshold);
        if (!Util.doublesEqual(viResult.cvar(), lpResult.cvar())) {
          log.log(Level.INFO, "Difference for threshold %f, LP: %.6f, VI: %.6f, delta: %g"
              .formatted(threshold, lpResult.cvar(), viResult.cvar(), Math.abs(lpResult.cvar() - viResult.cvar())));
        }
      }
    }

    Result.SspSolution sspSolution = new Result.SspSolution(sspConstructionTime, sspTime, optimalSsp);
    Result result = new Result(inputData, modelData, problemParameters, sspSolution, viSolution, lpSolution);
    GsonBuilder gsonBuilder = new GsonBuilder();
    if (outputFile.equals("-")) {
      gsonBuilder.setPrettyPrinting();
    }
    gsonBuilder.registerTypeAdapter(Duration.class, (JsonSerializer<Duration>) (src, typeOfSrc, context) ->
        new JsonPrimitive(src.toSeconds() + (src.getNano() / (double) TimeUnit.SECONDS.toNanos(1))));
    Gson gson = gsonBuilder.create();
    try (Writer writer = outputFile.equals("-") ? new OutputStreamWriter(System.out) : Files.newBufferedWriter(Path.of(outputFile))) {
      gson.toJson(result, writer);
    }
  }

  private static void buildLinearProgrammingModel(GRBModel model, double threshold, int varGuess, MDP mdp, Set<Integer> goalStates,
      double[] stateSSP, int[] stateCosts)
      throws GRBException {
    int states = mdp.size();

    GRBVar[][] stepStateProbability = new GRBVar[varGuess + 1][states];
    GRBVar[][][] stepActionProbability = new GRBVar[varGuess + 1][states][];

    for (int step = 0; step < varGuess; step++) {
      for (int state = 0; state < states; state++) {
        if (!goalStates.contains(state)) {
          stepStateProbability[step][state] = model.addVar(0.0, 1.0, 0.0, GRB.CONTINUOUS,
              "x_%d_%d".formatted(step, state));

          Distribution[] stateActions = mdp.transitions[state];
          GRBVar[] actionVars = new GRBVar[stateActions.length];
          GRBLinExpr actionSum = new GRBLinExpr();
          for (int action = 0; action < stateActions.length; action++) {
            if (stateActions[action] != null) {
              actionVars[action] = model.addVar(0.0, 1.0, 0.0, GRB.CONTINUOUS,
                  "a_%d_%d_%d".formatted(step, state, action));
              actionSum.addTerm(1.0, actionVars[action]);
            }
          }
          model.addConstr(stepStateProbability[step][state], '=', actionSum, "x_a_%d_%d".formatted(step, state));
          stepActionProbability[step][state] = actionVars;
        }
      }
    }

    GRBVar[] flowInsideGoal = new GRBVar[varGuess + 1];
    flowInsideGoal[0] = model.addVar(0.0, 0.0, 0.0, GRB.CONTINUOUS, "g_0");
    for (int i = 1; i < flowInsideGoal.length; i++) {
      flowInsideGoal[i] = model.addVar(0.0, 1.0, 0.0, GRB.CONTINUOUS, "g_%d".formatted(i));
    }

    for (int state = 0; state < states; state++) {
      if (!goalStates.contains(state)) {
        model.addConstr(stepStateProbability[0][state], '=', state == mdp.initialState ? 1.0 : 0.0, "init_%d".formatted(state));
      }
    }

    GRBLinExpr[] stepIncomingGoalFlow = new GRBLinExpr[varGuess + 1];
    Arrays.setAll(stepIncomingGoalFlow, i -> new GRBLinExpr());
    GRBLinExpr[][] incomingStateFlow = new GRBLinExpr[varGuess][states];
    for (int step = 1; step < varGuess; step++) {
      Arrays.setAll(incomingStateFlow[step], i -> new GRBLinExpr());
    }
    Map<Double, GRBLinExpr> flowRemainingCost = new HashMap<>();

    for (int step = 0; step < varGuess; step++) {
      for (int state = 0; state < states; state++) {
        if (goalStates.contains(state)) {
          continue;
        }
        int cost = stateCosts[state];

        Distribution[] stateActions = mdp.transitions[state];
        GRBVar[] actionVars = stepActionProbability[step][state];
        for (int actionIndex = 0; actionIndex < stateActions.length; actionIndex++) {
          Distribution action = stateActions[actionIndex];
          if (action == null) {
            continue;
          }

          GRBVar actionVar = actionVars[actionIndex];
          int _step = step;

          action.forEach((successor, probability) -> {
            GRBLinExpr expr;
            if (_step + cost < varGuess) {
              if (goalStates.contains(successor)) {
                expr = stepIncomingGoalFlow[_step + cost];
              } else {
                expr = incomingStateFlow[_step + cost][successor];
                assert expr != null;
              }
            } else {
              double remainingCost = _step + cost - varGuess + stateSSP[successor];
              assert (remainingCost == 0) == goalStates.contains(successor);
              expr = flowRemainingCost.computeIfAbsent(remainingCost, k -> new GRBLinExpr());
            }
            expr.addTerm(probability, actionVar);
          });
        }
      }
    }

    for (int step = 1; step < varGuess; step++) {
      GRBLinExpr[] stepIncomingFlows = incomingStateFlow[step];
      for (int state = 0; state < states; state++) {
        if (!goalStates.contains(state)) {
          if (stepIncomingFlows[state].size() > 0) {
            model.addConstr(stepStateProbability[step][state], '=', stepIncomingFlows[state], "t_%d_%d".formatted(step, state));
          } else {
            model.addConstr(stepStateProbability[step][state], '=', 0, "t_%d_%d".formatted(step, state));
          }
        }
      }
    }
    for (int step = 1; step <= varGuess; step++) {
      GRBLinExpr flowInGoalExpression = new GRBLinExpr(stepIncomingGoalFlow[step]);
      flowInGoalExpression.addTerm(1.0, flowInsideGoal[step - 1]);
      model.addConstr(flowInsideGoal[step], '=', flowInGoalExpression, "goal_%d".formatted(step));
    }

    model.addConstr(flowInsideGoal[varGuess - 1], '<', 1 - threshold, "cvar_lower");
    GRBLinExpr flowInGoalAtGuess = new GRBLinExpr(flowRemainingCost.getOrDefault(0.0, new GRBLinExpr()));
    flowInGoalAtGuess.addTerm(1.0, flowInsideGoal[varGuess]);
    model.addConstr(flowInGoalAtGuess, '>', 1 - threshold, "cvar_upper");

    for (Map.Entry<Double, GRBLinExpr> entry : flowRemainingCost.entrySet()) {
      double remainingCost = entry.getKey();
      if (remainingCost == 0.0) {
        continue;
      }
      long bits = Double.doubleToLongBits(remainingCost);
      GRBVar flowWithRemainingCost = model.addVar(0.0, 1.0, remainingCost, GRB.CONTINUOUS, "c_%d".formatted(bits));
      model.addConstr(flowWithRemainingCost, '=', entry.getValue(), "c_%d_c".formatted(bits));
    }
  }

  private static int[] computeMinimalReach(MDP mdp, Set<Integer> goalStates, double[] thresholds) {
    assert IntStream.range(0, thresholds.length - 1).allMatch(i -> thresholds[i] <= thresholds[i + 1]);
    int states = mdp.size();

    double[] stateValues = new double[states];
    double[] nextStateValues = new double[states];

    for (int state : goalStates) {
      stateValues[state] = 1.0;
      nextStateValues[state] = 1.0;
    }
    int[] minimalStep = new int[thresholds.length];
    Arrays.fill(minimalStep, Integer.MAX_VALUE);
    int reachStep = 1;
    while (minimalStep[0] == Integer.MAX_VALUE) {
      double[] currentValues = stateValues;
      double[] nextValues = nextStateValues;

      for (int state = 0; state < states; state++) {
        if (!goalStates.contains(state)) {
          double maximum = 0.0;
          Distribution[] stateActions = mdp.transitions[state];
          for (Distribution action : stateActions) {
            if (action != null) {
              double value = action.sum(currentValues);
              if (value > maximum) {
                maximum = value;
              }
            }
          }
          nextValues[state] = maximum;
        }
      }
      for (int tIndex = 0; tIndex < thresholds.length; tIndex++) {
        if (minimalStep[tIndex] == Integer.MAX_VALUE
            && nextValues[mdp.initialState] >= 1 - thresholds[tIndex]) {
          minimalStep[tIndex] = reachStep;
        }
      }

      stateValues = nextValues;
      nextStateValues = currentValues;
      reachStep += 1;
    }
    return minimalStep;
  }
}
