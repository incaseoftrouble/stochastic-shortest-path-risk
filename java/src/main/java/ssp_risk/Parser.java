package ssp_risk;

import static com.google.common.base.Preconditions.checkArgument;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSetMultimap;
import com.google.common.collect.Iterators;
import com.google.common.collect.SetMultimap;
import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.Collectors;

public final class Parser {
  private Parser() {}

  public static SetMultimap<String, Integer> parseLabels(Path labelPath) throws IOException {
    try (BufferedReader labelReader = Files.newBufferedReader(labelPath)) {
      String labelLine = labelReader.readLine().strip();
      Map<Integer, String> labelIndices = new HashMap<>();

      for (String labelEntry : labelLine.split(" ")) {
        String[] labelData = labelEntry.split("=");
        int index = Integer.parseInt(labelData[0]);
        String labelName = labelData[1].substring(1, labelData[1].length() - 1);
        checkArgument(!labelIndices.containsKey(index));
        labelIndices.put(index, labelName);
      }

      ImmutableSetMultimap.Builder<String, Integer> stateLabels = ImmutableSetMultimap.builder();
      labelReader.lines().map(String::strip).forEach(line -> {
        String[] lineData = line.split(" ");
        int state = Integer.parseInt(lineData[0].substring(0, lineData[0].length() - 1));
        for (int i = 1; i < lineData.length; i++) {
          int labelIndex = Integer.parseInt(lineData[i]);
          String label = labelIndices.get(labelIndex);
          stateLabels.put(label, state);
        }
      });
      return stateLabels.build();
    }
  }

  public static MDP parseMdp(Path transitionMatrixPath, SetMultimap<String, Integer> stateLabels) throws IOException {
    int initialState = Iterators.getOnlyElement(stateLabels.get("init").iterator());

    List<Map<Integer, SortedMap<Integer, Double>>> transitions;
    try (BufferedReader transitionReader = Files.newBufferedReader(transitionMatrixPath)) {
      String[] metadata = transitionReader.readLine().split(" ");
      checkArgument(metadata.length == 3);
      int states = Integer.parseInt(metadata[0]);

      transitions = new ArrayList<>(states);
      for (int i = 0; i < states; i++) {
        transitions.add(new HashMap<>());
      }
      transitionReader.lines().map(String::strip).filter(string -> !string.isEmpty()).map(string -> string.split(" ")).forEach(row -> {
        checkArgument(4 <= row.length && row.length <= 5);
        int source = Integer.parseInt(row[0]);
        int action = Integer.parseInt(row[1]);
        int dest = Integer.parseInt(row[2]);
        double probability = Double.parseDouble(row[3]);

        Double oldProbability = transitions.get(source).computeIfAbsent(action, k -> new TreeMap<>()).put(dest, probability);
        checkArgument(oldProbability == null);
      });
    }

    Distribution[][] transitionArray = new Distribution[transitions.size()][];
    ListIterator<Map<Integer, SortedMap<Integer, Double>>> iterator = transitions.listIterator();
    while (iterator.hasNext()) {
      int state = iterator.nextIndex();
      Map<Integer, SortedMap<Integer, Double>> stateActions = iterator.next();

      int maximalAction = stateActions.keySet().stream().mapToInt(Integer::intValue).max().orElse(0);
      transitionArray[state] = new Distribution[maximalAction + 1];
      for (Map.Entry<Integer, SortedMap<Integer, Double>> entry : stateActions.entrySet()) {
        int action = entry.getKey();
        SortedMap<Integer, Double> actionTransitions = entry.getValue();
        if (!actionTransitions.isEmpty()) {
          transitionArray[state][action] = Distribution.of(actionTransitions);
        }
      }
    }

    return new MDP(initialState, transitionArray);
  }

  public static Map<Integer, Integer> parseStateRewards(Path stateRewardsPath) throws IOException {
    ImmutableMap.Builder<Integer, Integer> rewards = ImmutableMap.builder();
    try (BufferedReader rewardReader = Files.newBufferedReader(stateRewardsPath)) {
      rewardReader.readLine();
      rewardReader.lines().map(String::strip).filter(string -> !string.isEmpty()).map(string -> string.split(" ")).forEach(row -> {
        checkArgument(row.length == 2);
        int state = Integer.parseInt(row[0]);
        int reward = Integer.parseInt(row[1]);
        checkArgument(reward >= 0);
        rewards.put(state, reward);
      });
    }
    return rewards.build();
  }

  public static Set<Integer> getAbsorbingStates(MDP mdp) {
    return mdp.stateStream().filter(state -> Arrays.stream(mdp.transitions[state])
            .allMatch(distribution -> distribution == null || Util.doublesEqual(distribution.value(state), 1.0)))
        .boxed()
        .collect(Collectors.toSet());
  }
}
