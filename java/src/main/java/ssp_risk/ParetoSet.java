package ssp_risk;

import static com.google.common.base.Preconditions.checkArgument;

import com.google.common.collect.Iterables;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Map;
import java.util.NavigableMap;
import java.util.TreeMap;
import java.util.stream.IntStream;
import javax.annotation.Nullable;

public record ParetoSet(double[] p, double[] e) {
  // TODO Navigable map for primitives :-(

  public static ParetoSet combine(ParetoSet[] stateSets, Distribution[] distributions) {
    // Building the Minkowski sum

    if (distributions.length == 1 && distributions[0].size() == 1) {
      // Simple case
      return stateSets[distributions[0].support(0)];
    }

    // Track an approximation of all points that are part of the hull
    NavigableMap<Double, Double> map = new TreeMap<>();

    for (Distribution distribution : distributions) {
      int successors = distribution.size();
      int[] support = distribution.support();
      if (successors == 1) {
        ParetoSet set = stateSets[support[0]];
        for (int i = 0; i < set.p.length; i++) {
          insert(map, set.p[i], set.e[i]);
        }
        continue;
      }

      double[] probability = distribution.probability();
      ParetoSet[] sets = new ParetoSet[successors];

      // maxP + minE give the maximal probability / minimal expectation that can be achieved with all successors up to index i
      // This is used to prune the sum computation
      double[] maxP = new double[successors];
      double[] minE = new double[successors];
      double maxPSum = 0.0;
      double minESum = 0.0;
      for (int index = 0; index < successors; index++) {
        ParetoSet set = stateSets[support[index]];
        sets[index] = set;
        double w = probability[index];
        maxPSum += set.maxP() * w;
        maxP[index] = maxPSum;
        minESum += set.minE() * w;
        minE[index] = minESum;
      }

      int[] indices = new int[successors];
      // pVector + eVector are partial sums of all successors >= i
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
    // Insert if this point can be a part of the hull
    @Nullable
    Map.Entry<Double, Double> entry = map.ceilingEntry(p);
    if (entry == null || e <= entry.getValue()) {
      map.put(p, e);
    }
  }

  private static int increaseCartesianIndices(int start, int[] indices, ParetoSet[] sets) {
    Arrays.fill(indices, 0, 0, start);

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
    assert Util.doublesEqual(checkP, pVector[0]) : "Expected %f got %f".formatted(checkP, pVector[0]);
    assert Util.doublesEqual(checkE, eVector[0]) : "Expected %f got %f".formatted(checkE, eVector[0]);
    return true;
  }

  private static boolean updateVector(ParetoSet[] sets, Distribution distribution, int pos, int[] indices,
      double[] pVector, double[] eVector, NavigableMap<Double, Double> map, double[] maxP, double[] minE) {
    // TODO Reliability of this could be improved by adding Kahan sum or similar

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
    // Branch pruning: If the current partial sum + what the next successors possibly can yield is worse than what is already in the sum
    // skip the entire branch
    if (s <= 2) {
      return s - 1;
    }
    double branchMaxP = maxP[s - 1] + pVector[s];
    double branchMinE = minE[s - 1] + eVector[s];
    Map.Entry<Double, Double> bestEntry = map.ceilingEntry(branchMaxP);
    if (bestEntry != null && branchMinE > bestEntry.getValue()) {
      // Prune this branch
      return increaseCartesianIndices(s, indices, sets);
    }
    return s - 1;
  }

  // Reference implementation for Minkowski sum

  public static boolean checkCombination(ParetoSet set, ParetoSet[] stateSets, Distribution[] distributions) {
    ParetoSet simple = combineReference(stateSets, distributions);
    assert set.similar(simple) : "\nImplementation: " + set + "\nReference:      " + simple;
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
    int size = points.size();
    checkArgument(size >= 1);
    if (size == 1) {
      Map.Entry<Double, Double> point = Iterables.getOnlyElement(points.entrySet());
      return of(point.getKey(), point.getValue());
    }

    // Build the pareto frontier in one pass over the sorted points

    Iterator<Map.Entry<Double, Double>> iterator = points.entrySet().iterator();
    Map.Entry<Double, Double> first = iterator.next();
    int sizeGuess = Math.max(size / 8, 4);
    double[] hullP = new double[sizeGuess];
    hullP[0] = first.getKey();
    double[] hullE = new double[sizeGuess];
    hullE[0] = first.getValue();
    double smallestExpectation = hullE[0];
    int top = 0;
    while (iterator.hasNext()) {
      Map.Entry<Double, Double> entry = iterator.next();
      double key = entry.getKey();
      assert hullP[top] < key;
      double value = entry.getValue();

      if (Util.lessOrEqual(value, smallestExpectation)) {
        if (value <= smallestExpectation) {
          smallestExpectation = value;
        }
        // Found a new lowest point, completely reset
        top = 0;
        hullP[0] = key;
        hullE[0] = value;
      } else {
        while (top > 0 && (value <= hullE[top]
            || Util.doubleLessOrEqualZero(orientation(hullP[top - 1], hullE[top - 1], hullP[top], hullE[top], key, value)))) {
          top -= 1;
        }
        assert (hullP[top] < key && hullE[top] <= value);

        if (Util.doublesEqual(key, hullP[top])) {
          hullP[top] = Math.max(key, hullP[top]);
          hullE[top] = Math.min(value, hullE[top]);
        } else {
          top += 1;
          if (top == hullP.length) {
            hullP = Arrays.copyOf(hullP, Math.min(hullP.length * 2, size));
            hullE = Arrays.copyOf(hullE, Math.min(hullE.length * 2, size));
          }
          hullP[top] = key;
          hullE[top] = value;
        }
      }
    }

    if (hullP.length == top + 1) {
      return new ParetoSet(hullP, hullE);
    }
    return new ParetoSet(Arrays.copyOf(hullP, top + 1), Arrays.copyOf(hullE, top + 1));
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
    assert IntStream.range(0, p.length - 2).allMatch(i -> orientation(p[i], e[i], p[i + 1], e[i + 1], p[i + 2], e[i + 2]) > 0.0) : this;
    assert IntStream.range(0, p.length - 1).allMatch(i -> p[i] < p[i + 1]) : this;
    assert IntStream.range(0, p.length - 1).allMatch(i -> e[i + 1] > e[i]) : this;
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
      if (Util.doublesEqual(p, this.p[this.p.length - 1])) {
        return this.e[this.p.length - 1];
      }
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
    return Util.lessOrEqual(bestExpectation(p), e);
  }

  public boolean isSimple() {
    //noinspection FloatingPointEquality
    return p.length == 1 && (p[0] == 0.0 || p[0] == 1.0);
  }
}
