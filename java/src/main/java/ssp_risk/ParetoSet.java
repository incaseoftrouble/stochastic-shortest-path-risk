package ssp_risk;

import java.util.Arrays;
import java.util.Comparator;
import java.util.Objects;
import java.util.PriorityQueue;
import java.util.stream.IntStream;

public record ParetoSet(double[] p, double[] e) {
  private static boolean isConvex(double[] x, double[] y) {
    return isConvex(x, y, x.length);
  }

  private static boolean isConvex(double[] x, double[] y, int limit) {
    return IntStream.range(0, limit - 2).allMatch(i ->
        Util.doubleLessOrEqualZero(-orientation(x[i], y[i], x[i + 1], y[i + 1], x[i + 2], y[i + 2])));
  }

  private static boolean isIncreasing(double[] x) {
    return isIncreasing(x, x.length);
  }

  private static boolean isIncreasing(double[] x, int limit) {
    return IntStream.range(0, limit - 1).allMatch(i -> Util.doublesLessOrEqual(x[i], x[i + 1]));
  }

  public static ParetoSet combine(ParetoSet[] stateSets, Distribution[] distributions) {
    int actions = distributions.length;
    if (actions == 1 && distributions[0].size() == 1) {
      return stateSets[distributions[0].support(0)];
    }

    // For each action, the points in the Minkowski sum
    int[] sizes = new int[actions];
    double[][] ps = new double[actions][];
    double[][] es = new double[actions][];

    for (int action = 0; action < actions; action++) {
      // Building the Minkowski for each action

      Distribution distribution = distributions[action];
      int successors = distribution.size();
      int[] support = distribution.support();

      if (successors == 1) {
        ParetoSet set = stateSets[support[0]];
        ps[action] = set.p;
        es[action] = set.e;
        continue;
      }

      // Build line segments in the Minkowski sum
      // First point is the sum of all lowest points of each successor set
      // Then always add the segment which increases the angle by the least amount
      // This implementation uses Java streams for a stable summation algorithm

      double[] probability = distribution.probability();
      ParetoSet[] sets = new ParetoSet[successors];
      Arrays.setAll(sets, i -> stateSets[support[i]]);

      // Segment indices for each of the considered successors
      int[] positions = new int[successors];
      PriorityQueue<Integer> priorities = new PriorityQueue<>(successors, Comparator.<Integer>comparingDouble(i -> {
        // Sort by slope as a proxy for angle
        int s = positions[i];
        ParetoSet set = sets[i];
        return (set.e[s + 1] - set.e[s]) / (set.p[s + 1] - set.p[s]);
      }).thenComparingDouble(i -> {
        // Then sort by length of the segment to be added
        int s = positions[i];
        ParetoSet set = sets[i];
        return set.p[s + 1] * probability[i];
      }));
      IntStream.range(0, successors).filter(i -> sets[i].size() > 1).forEach(priorities::add);

      // Summed segments
      int segments = Arrays.stream(sets).mapToInt(ParetoSet::size).sum() - sets.length + 1;
      double[] p = new double[segments];
      p[0] = IntStream.range(0, successors).mapToDouble(i -> sets[i].p[0] * probability[i]).sum();
      double[] e = new double[segments];
      e[0] = IntStream.range(0, successors).mapToDouble(i -> sets[i].e[0] * probability[i]).sum();

      int pos = 0;
      for (int r = 1; r < segments; r++) {
        int index = Objects.requireNonNull(priorities.poll());
        positions[index] += 1;
        if (positions[index] < sets[index].size() - 1) {
          // We could instead do something like pSum += next - current or similar, however this would accumulate errors
          priorities.add(index);
        }

        double pSum = IntStream.range(0, successors).mapToDouble(i -> sets[i].p[positions[i]] * probability[i]).sum();
        double eSum = IntStream.range(0, successors).mapToDouble(i -> sets[i].e[positions[i]] * probability[i]).sum();
        assert Util.doublesLessOrEqual(p[pos], pSum);
        if (pSum <= p[pos]) {
          // This could happen due to float imprecision in the sum
          e[pos] = Math.min(e[pos], eSum);
        } else {
          if (pos < 2 || orientation(p[pos - 1], e[pos - 1], p[pos], e[pos], pSum, eSum) > 0.0) {
            // Discard co-linear points
            pos += 1;
          }
          p[pos] = pSum;
          e[pos] = eSum;
        }
      }
      assert priorities.isEmpty();

      sizes[action] = pos + 1;
      ps[action] = p;
      es[action] = e;

      assert isConvex(p, e, pos + 1) : "Not convex:\n%s\n%s".formatted(Arrays.toString(p), Arrays.toString(e));
      assert isIncreasing(p, pos + 1) && isIncreasing(e, pos + 1)
          : "Consistency error:\n%s\n%s".formatted(Arrays.toString(p), Arrays.toString(e));
    }

    int[] indices = new int[actions];
    PriorityQueue<Integer> queue = new PriorityQueue<>(actions,
        Comparator.<Integer>comparingDouble(a -> ps[a][indices[a]]).thenComparingDouble(a -> es[a][indices[a]]));
    IntStream.range(0, actions).forEach(queue::add);

    assert !queue.isEmpty();
    int first = queue.poll();
    indices[first] += 1;
    if (indices[first] < sizes[first]) {
      queue.add(first);
    }

    int size = Arrays.stream(sizes).sum();
    int sizeGuess = Math.max(size / 8, 4);
    double[] hullP = new double[sizeGuess];
    hullP[0] = ps[first][0];
    double[] hullE = new double[sizeGuess];
    hullE[0] = es[first][0];
    double smallestExpectation = hullE[0];
    int top = 0;
    while (!queue.isEmpty()) {
      int a = queue.poll();
      double p = ps[a][indices[a]];
      double e = es[a][indices[a]];
      assert hullP[top] <= p;

      //noinspection FloatingPointEquality
      if (hullP[top] == p) {
        // This can happen if multiple successors have, e.g. pSum = 0.0
        // Since the queue sorts accordingly, we can skip these cases
        assert hullE[top] <= e;
      } else if (e <= smallestExpectation) {
        smallestExpectation = e;
        // Found a new lowest point, completely reset
        top = 0;
        hullP[0] = p;
        hullE[0] = e;
      } else {
        while (top > 0 && (e <= hullE[top] || orientation(hullP[top - 1], hullE[top - 1], hullP[top], hullE[top], p, e) <= 0.0)) {
          top -= 1;
        }
        assert (hullP[top] < p && hullE[top] <= e);
        top += 1;
        if (top == hullP.length) {
          hullP = Arrays.copyOf(hullP, Math.min(hullP.length * 2, size));
          hullE = Arrays.copyOf(hullE, Math.min(hullE.length * 2, size));
        }
        hullP[top] = p;
        hullE[top] = e;
      }

      indices[a] += 1;
      if (indices[a] < sizes[a]) {
        queue.add(a);
      }
    }

    if (top + 1 == hullP.length) {
      return new ParetoSet(hullP, hullE);
    }
    return new ParetoSet(Arrays.copyOf(hullP, top + 1), Arrays.copyOf(hullE, top + 1));
  }

  // Convex hull

  private static double orientation(double px, double py, double qx, double qy, double rx, double ry) {
    return (qx - px) * (ry - py) - (rx - px) * (qy - py);
  }

  public static ParetoSet of(double p, double e) {
    return new ParetoSet(new double[] {p}, new double[] {e});
  }

  // Pareto set

  @SuppressWarnings("AssignmentOrReturnOfFieldWithMutableType")
  public ParetoSet(double[] p, double[] e) {
    assert p.length == e.length;

    this.p = p;
    this.e = e;

    assert isConvex(p, e) : this;
    assert isIncreasing(p) : this;
    assert isIncreasing(e) : this;
  }

  public int size() {
    return p.length;
  }

  public double bestExpectation(double p) {
    if (this.p.length == 1) {
      return Util.doublesLessOrEqual(p, this.p[0]) ? e[0] : Double.NaN;
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
    return Util.doublesLessOrEqual(bestExpectation(p), e);
  }

  public boolean isSimple() {
    //noinspection FloatingPointEquality
    return p.length == 1 && (p[0] == 0.0 || p[0] == 1.0);
  }
}
