package ssp_risk;

import static com.google.common.base.Preconditions.checkArgument;

import com.google.common.collect.Iterables;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Map;
import java.util.NavigableMap;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.IntStream;

public record ParetoSet(double[] p, double[] e) {
  // Building the Minkowski sum

  public static ParetoSet combine(ParetoSet[] stateSets, Distribution[] distributions) {
    if (distributions.length == 1 && distributions[0].size() == 1) {
      return stateSets[distributions[0].support(0)];
    }

    NavigableMap<Double, Double> map = new TreeMap<>();
    for (Distribution distribution : distributions) {
      int successors = distribution.size();
      if (successors == 1) {
        ParetoSet set = stateSets[distribution.support(0)];
        for (int i = 0; i < set.p.length; i++) {
          insert(map, set.p[i], set.e[i]);
        }
        continue;
      }

      ParetoSet[] sets = new ParetoSet[successors];
      double[] maxP = new double[successors];
      double[] minE = new double[successors];
      double maxPSum = 0.0;
      double minESum = 0.0;
      for (int index = 0; index < successors; index++) {
        ParetoSet set = stateSets[distribution.support(index)];
        sets[index] = set;
        double w = distribution.probability(index);
        maxPSum += set.maxP() * w;
        maxP[index] = maxPSum;
        minESum += set.minE() * w;
        minE[index] = minESum;
      }

      int[] indices = new int[successors];
      double[] pVector = new double[successors];
      double[] eVector = new double[successors];

      int pos = successors - 1;
      while (pos < successors && updateVector(sets, distribution, pos, indices, pVector, eVector, map, maxP, minE)) {
        assert checkConsistency(stateSets, distribution, indices, pVector, eVector);
        insert(map, pVector[0], eVector[0]);
        pos = increaseCartesianIndices(0, indices, sets);
      }
    }
    return of(map);
  }

  private static void insert(NavigableMap<Double, Double> map, double p, double e) {
    Map.Entry<Double, Double> entry = map.ceilingEntry(p);
    if (entry == null || e <= entry.getValue()) {
      map.put(p, e);
    }
  }

  private static int increaseCartesianIndices(int start, int[] indices, ParetoSet[] sets) {
    for (int i = 0; i < start; i++) {
      indices[i] = 0;
    }

    int inc = start;
    while (inc < indices.length) {
      indices[inc] += 1;
      int size = sets[inc].size();
      assert indices[inc] <= size;
      if (indices[inc] == size) {
        indices[inc] = 0;
        inc += 1;
      } else {
        break;
      }
    }
    return inc;
  }

  private static boolean checkConsistency(ParetoSet[] sets, Distribution distribution,
      int[] setIndices, double[] pVector, double[] eVector) {
    double checkP = 0.0;
    double checkE = 0.0;
    for (int s = 0; s < distribution.size(); s++) {
      int successor = distribution.support(s);
      double weight = distribution.probability(s);
      ParetoSet set = sets[successor];

      int setIndex = setIndices[s];
      checkP += set.p[setIndex] * weight;
      checkE += set.e[setIndex] * weight;
    }
    assert Util.doublesEqual(checkP, pVector[0]);
    assert Util.doublesEqual(checkE, eVector[0]);
    return true;
  }

  private static boolean updateVector(ParetoSet[] sets, Distribution distribution, int pos, int[] indices,
      double[] pVector, double[] eVector, NavigableMap<Double, Double> map, double[] maxP, double[] minE) {
    int s = pos;
    while (s >= 0) {
      double w = distribution.probability(s);
      ParetoSet set = sets[s];
      int i = indices[s];
      if (s == pVector.length - 1) {
        pVector[s] = set.p[i] * w;
        eVector[s] = set.e[i] * w;
      } else {
        pVector[s] = pVector[s + 1] + set.p[i] * w;
        eVector[s] = eVector[s + 1] + set.e[i] * w;
      }

      s = pruneBranch(s, indices, pVector, eVector, maxP, minE, map, sets);
      if (s == distribution.size()) {
        return false;
      }
    }
    return true;
  }

  private static int pruneBranch(int s, int[] indices, double[] pVector, double[] eVector,
      double[] maxP, double[] minE, NavigableMap<Double, Double> map, ParetoSet[] sets) {
    if (s <= 2) {
      return s - 1;
    }
    double branchMaxP = maxP[s - 1] + pVector[s];
    double branchMinE = minE[s - 1] + eVector[s];
    Map.Entry<Double, Double> bestEntry = map.ceilingEntry(branchMaxP);
    if (bestEntry != null && branchMinE > bestEntry.getValue()) {
      return increaseCartesianIndices(s, indices, sets);
    }
    return s - 1;
  }

  // Reference implementation for Minkowski sum

  public static boolean checkCombination(ParetoSet set, ParetoSet[] stateSets, Distribution[] distributions) {
    ParetoSet simple = combineReference(stateSets, distributions);
    assert set.similar(simple) : "\nArg: " + set + "\nRef: " + simple;
    return true;
  }

  public static ParetoSet combineReference(ParetoSet[] stateSets, Distribution[] distributions) {
    if (distributions.length == 1 && distributions[0].size() == 1) {
      return stateSets[distributions[0].support()[0]];
    }

    NavigableMap<Double, Double> map = new TreeMap<>();
    for (Distribution distribution : distributions) {
      int successors = distribution.size();
      int[] support = distribution.support();
      double[] probability = distribution.probability();

      if (successors == 1) {
        ParetoSet set = stateSets[support[0]];
        for (int i = 0; i < set.size(); i++) {
          map.merge(set.p[i], set.e[i], Double::min);
        }
        continue;
      }

      ParetoSet[] sets = new ParetoSet[successors];
      for (int s = 0; s < successors; s++) {
        sets[s] = stateSets[support[s]];
      }

      int[] indices = new int[successors];

      int pos;
      do {
        double p = 0.0;
        double e = 0.0;
        for (int s = 0; s < successors; s++) {
          double w = probability[s];
          ParetoSet set = sets[s];
          int i = indices[s];

          p += set.p[i] * w;
          e += set.e[i] * w;
        }

        map.merge(p, e, Double::min);
        pos = increaseCartesianIndices(0, indices, sets);
      } while (pos < successors);
    }

    return ParetoSet.of(map);
  }

  // Convex hull

  private static double orientation(double px, double py, double qx, double qy, double rx, double ry) {
    return (qx - px) * (ry - py) - (rx - px) * (qy - py);
  }

  public static ParetoSet of(NavigableMap<Double, Double> points) {
    checkArgument(points.size() >= 1);
    if (points.size() == 1) {
      Map.Entry<Double, Double> point = Iterables.getOnlyElement(points.entrySet());
      return of(point.getKey(), point.getValue());
    }

    double smallestExpectation = points.values().stream().mapToDouble(Double::doubleValue).min().orElseThrow();
    Map.Entry<Double, Double> bottomRightEntry = points.entrySet().stream()
        .filter(entry -> Util.doublesEqual(entry.getValue(), smallestExpectation))
        .max(Map.Entry.comparingByKey()).orElseThrow();
    SortedMap<Double, Double> largerPoints = points.tailMap(bottomRightEntry.getKey(), true);
    Iterator<Map.Entry<Double, Double>> iterator = largerPoints.entrySet().iterator();

    int sizeGuess = Math.max(points.size() / 8, 4);

    double[] hullP = new double[sizeGuess];
    hullP[0] = bottomRightEntry.getKey();
    double[] hulLE = new double[sizeGuess];
    hulLE[0] = bottomRightEntry.getValue();
    int top = 0;

    while (iterator.hasNext()) {
      Map.Entry<Double, Double> next = iterator.next();
      if (next.getKey() <= hullP[top]) {
        continue;
      }
      while (top >= 1 && (hulLE[top] >= next.getValue()
          || Util.lessZeroLoose(orientation(hullP[top - 1], hulLE[top - 1], hullP[top], hulLE[top],
          next.getKey(), next.getValue())))) {
        top -= 1;
      }
      assert hullP[top] < next.getKey() && hulLE[top] <= next.getValue();

      top += 1;
      if (top == hullP.length) {
        hullP = Arrays.copyOf(hullP, hullP.length * 2);
        hulLE = Arrays.copyOf(hulLE, hulLE.length * 2);
      }
      hullP[top] = next.getKey();
      hulLE[top] = next.getValue();
    }

    if (hullP.length == top + 1) {
      return new ParetoSet(hullP, hulLE);
    }
    return new ParetoSet(Arrays.copyOf(hullP, top + 1), Arrays.copyOf(hulLE, top + 1));
  }

  public static ParetoSet of(double p, double e) {
    return new ParetoSet(new double[] {p}, new double[] {e});
  }

  // Pareto set

  @SuppressWarnings("AssignmentOrReturnOfFieldWithMutableType")
  public ParetoSet(double[] p, double[] e) {
    this.p = p;
    this.e = e;

    assert p.length == e.length;
    assert IntStream.range(0, p.length - 1).allMatch(i -> orientation(p[i], e[i], p[i + 1], e[i + 1], p[i + 2], e[i + 2]) > 0.0) : this;
    assert IntStream.range(0, p.length).allMatch(i -> p[i] < p[i + 1]) : this;
    assert IntStream.range(0, p.length).allMatch(i -> e[i + 1] > e[i]) : this;
  }

  public int size() {
    return p.length;
  }

  public double maxP() {
    return p[p.length - 1];
  }

  private double minE() {
    return e[0];
  }

  public double bestExpectation(double p) {
    if (this.p.length == 1) {
      return p <= this.p[0] ? e[0] : Double.NaN;
    }

    int index = Arrays.binarySearch(this.p, p);
    if (index >= 0) {
      return e[index];
    }
    index = -(index + 1);
    assert index <= this.p.length;
    if (index == this.p.length) {
      assert p > this.p[this.p.length - 1];
      return Double.NaN;
    }
    if (index == 0) {
      return e[0];
    }

    double pLeft = this.p[index - 1];
    double eLeft = e[index - 1];
    double pRight = this.p[index];
    double eRight = e[index];

    assert pLeft < p && p < pRight;
    assert eLeft <= eRight;
    return eLeft + (eRight - eLeft) * (p - pLeft) / (pRight - pLeft);
  }

  @Override
  public String toString() {
    StringBuilder builder = new StringBuilder(20 * p.length);
    builder.append('{');
    for (int i = 0; i < p.length; i++) {
      builder.append(p[i]).append(':').append(e[i]);
      if (i < p.length - 1) {
        builder.append(", ");
      }
    }
    builder.append('}');
    return builder.toString();
  }

  public boolean similar(ParetoSet o) {
    for (int i = 0; i < p.length; i++) {
      if (!o.containsPoint(this.p[i], this.e[i])) {
        return false;
      }
    }
    for (int i = 0; i < o.p.length; i++) {
      if (!containsPoint(o.p[i], o.e[i])) {
        return false;
      }
    }
    return true;
  }

  private boolean containsPoint(double p, double e) {
    int i = Arrays.binarySearch(this.p, p);
    if (i >= 0) {
      return Util.doublesEqual(this.e[i], e);
    }
    int j = -(i + 1);
    return (j > 0 && Util.doublesEqual(this.p[j - 1], p) && Util.doublesEqual(this.e[j - 1], e))
        || (j < this.p.length && Util.doublesEqual(this.p[j], p) && Util.doublesEqual(this.e[j], e));
  }

  public boolean isSimple() {
    //noinspection FloatingPointEquality
    return p.length == 1 && (p[0] == 0.0 || p[0] == 1.0);
  }
}
