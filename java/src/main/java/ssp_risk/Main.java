package ssp_risk;

import static com.google.common.base.Preconditions.checkArgument;
import static picocli.CommandLine.Command;
import static picocli.CommandLine.Option;

import com.google.common.collect.ImmutableSet;
import com.google.common.collect.SetMultimap;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonPrimitive;
import com.google.gson.JsonSerializer;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Duration;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.concurrent.Callable;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import javax.annotation.Nullable;
import picocli.CommandLine;

@Command(name = "solve", mixinStandardHelpOptions = true, version = "1.0")
public final class Main implements Callable<Integer> {
  private static final Logger log = Logger.getLogger("cvar");
  private final String[] args;

  @Option(names = {"-m", "--model"}, description = "Model transition file (.tra)", required = true)
  private Path transitionFile;
  @Option(names = {"-l", "--labels"}, description = "Model labels file (.lab)")
  @Nullable
  private Path labelsFile;
  @Option(names = {"-s", "--srew"}, description = "Model state rewards file (.srew)")
  @Nullable
  private Path stateRewardsFile;
  @Option(names = {"--trew"}, description = "Model transition rewards file (.trew)")
  @Nullable
  private Path transitionRewardsFile;
  @Option(names = {"-g", "--goal"}, description = "Goal label")
  @Nullable
  private String goalLabel;
  @Option(names = {"-t", "--thresholds"}, description = "CVaR thresholds", split = "\\|")
  private List<String> thresholdSpecs;
  @Option(names = {"--method"}, description = "Solution methods", split = ",")
  private Set<Method> methods;
  @Option(names = {"-o", "--output"}, description = "Output (- for stdout)")
  private String outputFile;
  @Option(names = {"--offset-rewards"}, description = "Add a constant to make all rewards positive", type = Boolean.class)
  private boolean offsetRewards;

  public Main(String[] args) {
    this.args = args;
  }

  @Override
  public Integer call() throws Exception {
    SortedSet<Double> thresholdSet = new TreeSet<>();
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
    Result.InputData inputData = new Result.InputData(transitionFile.getFileName().toString(), args);
    if (!Files.exists(transitionFile)) {
      Path transitionSibling = transitionFile.resolveSibling(transitionFile.getFileName() + ".tra");
      if (Files.exists(transitionSibling)) {
        if (labelsFile == null) {
          labelsFile = transitionFile.resolveSibling(transitionFile.getFileName() + ".lab");
        }
        if (stateRewardsFile == null) {
          Path stateRewards = transitionFile.resolveSibling(transitionFile.getFileName() + ".srew");
          if (Files.exists(stateRewards)) {
            stateRewardsFile = stateRewards;
          }
        }
        if (transitionRewardsFile == null) {
          Path transitionRewards = transitionFile.resolveSibling(transitionFile.getFileName() + ".trew");
          if (Files.exists(transitionRewards)) {
            transitionRewardsFile = transitionRewards;
          }
        }
        transitionFile = transitionSibling;
      }
    }

    SetMultimap<String, Integer> labels = Parser.parseLabels(labelsFile);
    MDP mdp = Parser.parseMdp(transitionFile, labels);

    Set<Integer> goalLabelStates = labels.get(goalLabel);
    Set<Integer> absorbingStates = Parser.getAbsorbingStates(mdp);
    Set<Integer> goalStates;
    if (goalLabelStates.isEmpty() && Objects.equals(goalLabel, "absorbing")) {
      goalStates = ImmutableSet.copyOf(absorbingStates);
    } else {
      goalStates = ImmutableSet.copyOf(goalLabelStates);
    }

    checkArgument(!goalStates.isEmpty(), "Did not find any goal states!");
    checkArgument(!goalStates.contains(mdp.initialState), "Initial state is goal state");

    int states = mdp.size();
    int[][] costs = new int[states][];
    if (stateRewardsFile == null && transitionRewardsFile == null) {
      for (int state = 0; state < states; state++) {
        int[] stateActionCosts = new int[mdp.transitions[state].length];
        Arrays.fill(stateActionCosts, goalStates.contains(state) ? 0 : 1);
        costs[state] = stateActionCosts;
      }
    } else {
      Map<Integer, Integer> stateCosts = stateRewardsFile == null
          ? Map.of() : Parser.parseStateRewards(stateRewardsFile);
      Map<Integer, Map<Integer, Integer>> actionCosts = transitionRewardsFile == null
          ? Map.of() : Parser.parseActionRewards(transitionRewardsFile);

      int offset = 0;
      if (offsetRewards) {
        int minimalCost = Integer.MAX_VALUE;
        for (int state = 0; state < states; state++) {
          if (goalStates.contains(state)) {
            continue;
          }
          int stateCost = stateCosts.getOrDefault(state, 0);
          var stateActionCostsMap = actionCosts.getOrDefault(state, Map.of());
          for (int action = 0; action < mdp.transitions[state].length; action++) {
            int cost = stateCost + stateActionCostsMap.getOrDefault(action, 0);
            if (cost < minimalCost) {
              minimalCost = cost;
            }
          }
        }
        offset = -minimalCost + 1;
      }
      if (offset != 0) {
        log.log(Level.INFO, "Offsetting rewards by {0}", offset);
      }

      for (int state = 0; state < states; state++) {
        int[] stateActionCosts = new int[mdp.transitions[state].length];
        if (!goalStates.contains(state)) {
          int stateCost = stateCosts.getOrDefault(state, 0);
          var stateActionCostsMap = actionCosts.getOrDefault(state, Map.of());
          for (int action = 0; action < mdp.transitions[state].length; action++) {
            stateActionCosts[action] = offset + stateCost + stateActionCostsMap.getOrDefault(action, 0);
          }
        }
        costs[state] = stateActionCosts;
      }
    }

    var modelData = new Result.ModelData(states, mdp.actionCount(), mdp.transitionCount(), mdp.initialState, absorbingStates.size());
    Result.ProblemParameters problemParameters = new Result.ProblemParameters(goalStates.size(), thresholds);
    int maximalCost = Arrays.stream(costs).flatMapToInt(Arrays::stream).max().orElseThrow();
    log.log(Level.INFO, "Problem has {0} states, {1} goal states (label {3}), {4} absorbing, maximal cost {5}, solving for thresholds {6}",
        new Object[] {states, goalStates.size(), goalLabelStates.size(), goalLabel, absorbingStates.size(),
            maximalCost, Arrays.stream(thresholds).mapToObj("%.3f"::formatted).collect(Collectors.toList())});

    var solution = new Solver(mdp, costs, methods, goalStates, thresholds).solve();
    Result result = new Result(inputData, modelData, problemParameters,
        solution.sspSolution(), solution.viSolution(), solution.lpSolution());

    GsonBuilder gsonBuilder = new GsonBuilder();
    if (outputFile == null || outputFile.equals("-")) {
      gsonBuilder.setPrettyPrinting();
    }
    gsonBuilder.registerTypeAdapter(Duration.class, (JsonSerializer<Duration>) (src, typeOfSrc, context) ->
        new JsonPrimitive(src.toSeconds() + (src.getNano() / (double) TimeUnit.SECONDS.toNanos(1))));
    Gson gson = gsonBuilder.create();
    try (Writer writer = (outputFile == null || outputFile.equals("-"))
        ? new OutputStreamWriter(System.out) : Files.newBufferedWriter(Path.of(outputFile))) {
      gson.toJson(result, writer);
    }
    return 0;
  }


  public static void main(String[] args) {
    System.exit(new CommandLine(new Main(args)).execute(args));
  }

  public enum Method {
    VI, LP
  }
}
