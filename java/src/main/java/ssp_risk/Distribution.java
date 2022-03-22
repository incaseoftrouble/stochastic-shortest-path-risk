package ssp_risk;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.Arrays;
import java.util.Map;
import java.util.SortedMap;
import java.util.stream.IntStream;

public record Distribution(int[] support, double[] probability) {
  public static Distribution of(SortedMap<Integer, Double> map) {
    int[] support = new int[map.size()];
    double[] probability = new double[support.length];
    int index = 0;
    for (Map.Entry<Integer, Double> entry : map.entrySet()) {
      support[index] = entry.getKey();
      probability[index] = entry.getValue();
      index += 1;
    }
    assert Util.doublesEqual(Arrays.stream(probability).sum(), 1.0) : "Invalid distribution " + Arrays.stream(probability).sum();
    return new Distribution(support, probability);
  }

  public static Distribution of(Map<Integer, Double> map) {
    if (map instanceof SortedMap<Integer, Double>) {
      return of((SortedMap<Integer, Double>) map);
    }

    int[] support = new int[map.size()];
    double[] probability = new double[support.length];
    int index = 0;
    for (Integer key : map.keySet()) {
      support[index] = key;
    }
    Arrays.sort(support);
    for (int i = 0; i < support.length; i++) {
      probability[i] = map.get(support[i]);
    }
    checkArgument(Util.doublesEqual(Arrays.stream(probability).sum(), 1.0));
    return new Distribution(support, probability);
  }

  public double value(int key) {
    for (int i = 0; i < support.length; i++) {
      if (support[i] == key) {
        return probability[i];
      }
    }
    return 0.0;
  }

  public double sum(double[] values) {
    return IntStream.range(0, support.length).mapToDouble(i -> values[support[i]] * probability[i]).sum();
  }

  public int support(int key) {
    return support[key];
  }

  public double probability(int key) {
    return probability[key];
  }

  public void forEach(DistributionAction action) {
    for (int i = 0; i < support.length; i++) {
      action.accept(support[i], probability[i]);
    }
  }

  public int size() {
    return support.length;
  }

  @Override
  public String toString() {
    StringBuilder builder = new StringBuilder(support.length * 20);
    builder.append('{');
    for (int i = 0; i < support.length; i++) {
      builder.append(support[i]).append(": ").append(probability[i]);
      if (i < support.length - 1) {
        builder.append(", ");
      }
    }
    builder.append('}');
    return builder.toString();
  }
}
