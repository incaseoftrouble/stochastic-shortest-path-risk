package ssp_risk;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.Arrays;
import java.util.Objects;
import java.util.Set;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public final class MDP {
  public final int initialState;
  public final Distribution[][] transitions;
  private Set<Integer> reachable;

  @SuppressWarnings("AssignmentOrReturnOfFieldWithMutableType")
  public MDP(int initialState, Distribution[][] transitions, Set<Integer> reachable) {
    assert reachable.contains(initialState);
    checkArgument(initialState < transitions.length);
    this.reachable = reachable;
    this.initialState = initialState;
    this.transitions = transitions;

  }

  public IntStream stateStream() {
    return reachable.size() == transitions.length ? IntStream.range(0, transitions.length) : reachable.stream().mapToInt(Integer::intValue);
  }

  public int size() {
    return transitions.length;
  }

  public long actionCount() {
    return Arrays.stream(transitions).flatMap(Arrays::stream).filter(Objects::nonNull).count();
  }

  public long transitionCount() {
    return Arrays.stream(transitions).flatMap(Arrays::stream).filter(Objects::nonNull).mapToInt(Distribution::size).asLongStream().sum();
  }

  public Stream<Distribution> stateActionStream(int state) {
    return Arrays.stream(transitions[state]).filter(Objects::nonNull);
  }
}
